#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
#include <algorithm>
using namespace std;


int argmin(double *array,int len){
	int minloc = std::distance(array, std::min_element(array,array+len));
	return minloc;
	
	
}

double error_spread(double value[2]) {
	//construct a trivial random generator engine from a time-based seed
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> distribution(value[0],value[1]);
	return distribution(generator);
}

double DHI_calc(double MHI) {
 	// Jing Wang, Koribalski, et al 2016 
	// HI Mass - Diameter relationship 
	// https://arxiv.org/abs/1605.01489

	double in_slope[2] = {0.506,0.003} ;
	double in_constant[2] = {-3.293,0.009} ;
	double scatr = 0.06 ;

	double constant = error_spread(in_constant);
	double slope = error_spread(in_slope);
	
	return pow(10.0,slope*log10(MHI)+constant);

};

double Mstar_calc(double MHI) {
	// Bradford et al 2015, Fig 5 
	// HI Mass - Stellar Mass Relationship 
	// https://arxiv.org/abs/1505.04819
	double Mgas = MHI*1.4;
	double in_slope[2] ;
	double in_constant[2] ; 
	double scatr[2] ;
	
	if (Mgas  < 9.2832) {
		in_slope[0] = 1.052 ;
		in_slope[1] = 0.058 ;
		in_constant[0] = 0.461 ;
		in_constant[1] = 0.011 ;
		scatr[0] = 0.285 ;
		scatr[1] = 0.019 ;
	} else {
		in_slope[0] = 0.461 ;
		in_slope[1] = 0.011 ;
		in_constant[0] = 5.329 ;
		in_constant[1] = 0.112 ;
		scatr[0] = 0.221 ;
		scatr[1] = 0.006 ;
	}

	double constant = error_spread(in_constant);
	double slope = error_spread(in_slope);
	double scatter = error_spread(scatr);

	double Mstar = slope*log10(Mgas)+constant;
	double spread[2] = {Mstar,scatter};
	Mstar = error_spread(spread);
	
	return Mstar;
	
}	

double BTFR(double Mbar) {
	// Bradford et al 2015, Fig 6
	// Baryonic Tully-Fisher relationship
	// https://arxiv.org/abs/1505.04819
	
	double in_slope[2] = {0.277,0.004};
	double in_constant[2] = {-0.672,0.041};
	double scatr[2] = {0.075,0.002};

	double constant = error_spread(in_constant);
	double slope = error_spread(in_slope);
	double scatter = error_spread(scatr);

	double logv = log10(Mbar) * slope + constant;
	double spread[2] = {logv,scatter};
	logv = error_spread(spread);

	return pow(10.0,logv);
}

double Ropt_calc(double vflat) {
	// Saintonge et al 2007
	// Optical Radius (r83) - Vflat Relationship
	// https://arxiv.org/abs/0710.0760
	
	double in_slope[2] = {0.56,0.04};
	double in_constant[2] = {-0.36,0.08};
	double scatr = 0.16;

	double constant = error_spread(in_constant);
        double slope = error_spread(in_slope);

	double Ropt = constant + slope * log10(vflat);
	double spread[2] = {Ropt,scatr};
	Ropt = error_spread(spread);

	return pow(10.0,Ropt);

}

double V0_func(double x){
	// best fitting parameters for V0
	// Polyex turnover velocity
	// based on Catinella et al 2007
	// https://arxiv.org/abs/astro-ph/0512051
	// Best-fitting parameters from SciPy
	double a = 2.25799719;
	double b = 0.40588194;
	double c = 4.94486961;
	double d = -1.88611724;
	return a*exp(-x*b -c)+d*x;
}

double a_func(double x){
	// best fitting parameters for a
	// Polyex outer slope
	// based on Catinella et al 2007
	// https://arxiv.org/abs/astro-ph/0512051
	// Best-fitting parameters from SciPy
	double a = 0.00766497;
	double b = 0.17990300;
	return a*x+b;
}

double rt_func(double x){
	// best fitting parameters for rPE
	// Polyex turnover radius
	// based on Catinella et al 2007
	// https://arxiv.org/abs/astro-ph/0512051
	// Best-fitting parameters from SciPy
	double a = -2.07573321;
	double b = -2.97620420;
	return a*log10(-x)-b;
}

double Mag_calc(double vrot, double Ropt, double RHI, double mstar){
	// Find Mag, and slope based on
	// Catinella et al 2007 (Mag)
    	// https://arxiv.org/abs/astro-ph/0512051
    	// Dutton et al 2018 (Slope)
    	// https://arxiv.org/abs/1807.10518
	
	// Create range of magnitudes
	// and set parameters for all mags
	int Mag_length = 25000;
	int guess_a = 440;
	double Mag[Mag_length];
	double vt[Mag_length];
	double vt_0[Mag_length];
	double rt[Mag_length];
	double vrot_compare[Mag_length];
	double a[Mag_length];

	for (int i=0;i <= Mag_length ; i++){
		Mag[i] = -25.0 + 0.001*i;
		vt_0[i] = V0_func(Mag[i]);
		a[i] = a_func(Mag[i]);
		rt[i] = rt_func(Mag[i]);
	}
	
	// Set slope from NIHAO 17
	double slope_sparc = 0.123 - 0.137*(log10(mstar)-9.471) ;
	double spread[2] = {slope_sparc,0.19};
	slope_sparc = error_spread(spread);

	double x1,x2;
	
	// This is where I'd begin the loop

	// Outer edge, and half of it for the slope
	x2 = RHI * 4.0/3.0;
	x1 = RHI * 2.0/3.0;

	// Calculate rotation velocities at Ropt for all vt_0, rt
	for (int i=0; i<= Mag_length;i++){
		vt[i] = vt_0[i] * (1.0 - exp(-Ropt/rt[i])) * (1.0 + a[i] * Ropt/rt[i]);
		//cout << vt_0[i] << " " << vt[i] << endl;
	}
	
	for (int i=0; i<= Mag_length;i++) {vrot_compare[i] = abs(vrot - vt[i]); }

	// Best guess for Magnitude based on vrot with other params
	// Finds index of vt that most closely matches vrot and
	// that matches the Magnitude
	int ind = argmin(vrot_compare,Mag_length);
	double Mag_guess = Mag[ind];
	double vt_guess = vt[ind];
	double rt_guess = rt[ind];
	double vt_0_guess = vt_0[ind];
	double a_guess[guess_a];
	double slope1[guess_a];
	double slope2[guess_a];
	double slope1_log[guess_a];
	double slope2_log[guess_a];
	double slope[guess_a];
	double slope_sparc_arr[guess_a];

	cout << vt[ind] << "  " << vrot << endl;

	// Consider a range of values of alpha

	for (int i=0; i<= guess_a;i++) {
		a_guess[i] =-0.04 + 0.001*i;
		slope1[i] = ((1.-exp(-x2/rt_guess))*(1.+a_guess[i]*x2/rt_guess));
		slope2[i] = ((1.-exp(-x1/rt_guess))*(1.+a_guess[i]*x1/rt_guess));
		slope_sparc_arr[i] = slope_sparc;

		if(slope1[i] > 0 and slope2[i] > 0){
			slope1_log[i] = log10(slope1[i]);
			slope2_log[i] = log10(slope2[i]);
			slope[i] = (slope1_log[i] - slope2_log[i])/(log10(x2)-log10(x1));

		} else { 
			slope1_log[i]=9999; 
			slope2_log[i]=-9999;
		}
		slope_sparc_arr[i]=abs(slope_sparc-slope[i]) ;
		cout << " " << a_guess[i]<<" "<<slope[i] << " " << slope_sparc_arr[i] << endl;
	}
	//for (int i=0; i<= guess_a;i++) {slope_sparc_arr[i]=abs(slope_sparc-slope[i]) ; cout << " " << a_guess[i]<<" "<<slope[i] << " " << slope_sparc_arr[i] << endl;}
	cout << slope[argmin(slope_sparc_arr,guess_a)]<<" " << slope_sparc_arr[argmin(slope_sparc_arr,guess_a)] << endl;
	cout << vt[ind] << "  " << vrot << endl;



	return 0;

}


int setup_relations(float mass,float beams) {
	double MHI = pow(10.0,9.5);
	double DHI = DHI_calc(MHI) ;
	double Mstar = Mstar_calc(MHI);
	double vflat = BTFR(Mstar + 1.4*MHI);
	double Rs = (DHI/2.0) * 0.18;
	double Ropt = Ropt_calc(vflat);
	/*
	cout << DHI << "\n";

	cout << Mstar << "\n";

	cout << vflat << "\n";

	cout << Ropt  << "\n";
	*/
	Mag_calc(vflat,Ropt,DHI/2.0,Mstar);

	return 0;
}

