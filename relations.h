/*
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include "relations.h"
*/
using namespace std;

int deffile_print(double *radi,double *vrot,double *sbr,double *z,int radi_len, double inclination);

int argmin(double *array,int len){
	int minloc;
	minloc = std::distance(array, std::min_element(array,array+len));
	return minloc;
}

double arcsec_to_rad(double value){
	return value*(3.14159265359/648000.);
}

double error_spread(double value[2], bool scatter_flag) {
	//construct a true random number
	random_device rd;
        mt19937 gen(rd());
	normal_distribution<double> dis(value[0],value[1]);
	
	double num ; 
	if (scatter_flag){
		return dis(gen);
	} else {
	return(value[0]);
	}
}

double DHI_calc(double MHI, bool scatter_flag) {
 	// Jing Wang, Koribalski, et al 2016 
	// HI Mass - Diameter relationship 
	// https://arxiv.org/abs/1605.01489

	double in_slope[2] = {0.506,0.003} ;
	double in_constant[2] = {-3.293,0.009} ;
	double scatr = 0.06 ;

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);
	
	return pow(10.0,slope*log10(MHI)+constant);

};

double Mstar_calc(double MHI, bool scatter_flag) {
	// Bradford et al 2015, Fig 5 
	// HI Mass - Stellar Mass Relationship 
	// https://arxiv.org/abs/1505.04819
	double Mgas = MHI*1.4;
	double in_slope[2] ;
	double in_constant[2] ; 
	double scatr[2] ;
	Mgas = log10(Mgas);

	if (Mgas  < 9.2832) {
		in_slope[0] = 1.052 ;
		in_slope[1] = 0.058 ;
		in_constant[0] = 0.236 ;
		in_constant[1] = 0.476 ;
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

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);
	double scatter = error_spread(scatr,scatter_flag);

	double Mstar = Mgas*1./slope-constant/slope;
	double spread[2] = {Mstar,scatter};

	return pow(10.0,Mstar);
	
}	

double BTFR(double Mbar, bool scatter_flag) {
	// Bradford et al 2015, Fig 6
	// Baryonic Tully-Fisher relationship
	// https://arxiv.org/abs/1505.04819
	
	double in_slope[2] = {0.277,0.004};
	double in_constant[2] = {-0.672,0.041};
	double scatr[2] = {0.075,0.002};

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);
	double scatter = error_spread(scatr,scatter_flag);

	double logv = log10(Mbar) * slope + constant;
	double spread[2] = {logv,scatter};
	logv = error_spread(spread,scatter_flag);

	return pow(10.0,logv);
}

double Ropt_calc(double vflat, bool scatter_flag) {
	// Saintonge et al 2007
	// Optical Radius (r83) - Vflat Relationship
	// https://arxiv.org/abs/0710.0760
	
	double in_slope[2] = {0.56,0.04};
	double in_constant[2] = {-0.36,0.08};
	double scatr = 0.16;
	double h = 0.70;

	double constant = error_spread(in_constant,scatter_flag);
        double slope = error_spread(in_slope,scatter_flag);

	double Ropt = constant + slope * log10(vflat);
	double spread[2] = {Ropt,scatr};
	Ropt = error_spread(spread,scatter_flag);

	return h*pow(10.0,Ropt);

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
	double a = 0.04264749;
	double b = 1.12975029;
	return a*x+b;
}

struct Mag_params {double Mag; double alpha; double slope;};


Mag_params Mag_calc(double vrot, double Ropt, double RHI, double mstar,bool scatter_flag){
	// Find Mag, and slope based on
	// Catinella et al 2007 (Mag)
    	// https://arxiv.org/abs/astro-ph/0512051
    	// Dutton et al 2018 (Slope)
    	// https://arxiv.org/abs/1807.10518
	
	// Create range of magnitudes
	// and set parameters for all mags
	int Mag_length = 24000;
	int guess_a = 4400;
	double Mag[Mag_length];
	double vt[Mag_length];
	double vt_0[Mag_length];
	double rt[Mag_length];
	double vrot_compare[Mag_length];
	double a[Mag_length];
	double a_temp;
	double Mag_guess;

	// Set slope from NIHAO 17
	double slope_sparc = 0.123 - 0.137*(log10(mstar)-9.471) ;
	double spread[2] = {slope_sparc,0.19};
	slope_sparc = error_spread(spread,scatter_flag);

	double x1,x2;
	
	// This is where I'd begin the loop
	for (int j=0; j<2; j++){
		for (int i=0;i <= Mag_length ; i++){
			Mag[i] = -24.0 + 0.001*i;
			vt_0[i] = V0_func(Mag[i]);
			rt[i] = Ropt * rt_func(Mag[i]);
			if (j == 0) {
				a[i] = a_func(Mag[i]);
			} else {
				a[i] = a_temp;
			}
		}

		// Outer edge, and half of it for the slope
		x2 = RHI * 1.0;
		x1 = RHI * 0.5;
	
		// Calculate rotation velocities at Ropt for all vt_0, rt
		for (int i=0; i<= Mag_length;i++){
			vt[i] = vt_0[i] * (1.0 - exp(-Ropt/rt[i])) * (1.0 + a[i] * Ropt/rt[i]);
			//if (abs(vt[i] - vrot) < 0.1) {cout << vt[i] << " " << vrot << " " << i << " " << Mag[i] << " " << a[i] << " " << rt[i] << " " << vt_0[i] << " " << Ropt << endl;}
		}
		
		for (int i=0; i<= Mag_length;i++) {vrot_compare[i] = abs(vrot - vt[i]); }
	
		// Best guess for Magnitude based on vrot with other params
		// Finds index of vt that most closely matches vrot and
		// that matches the Magnitude
		int ind = argmin(vrot_compare,Mag_length);
		Mag_guess = Mag[ind];
		// cout << Mag_guess << endl;
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
	
		// Consider a range of values of alpha
		for (int i=0; i<= guess_a;i++) {
			a_guess[i] =0.00 + 0.0001*i;
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
		}
		for (int i=0; i<= guess_a;i++) {
			slope_sparc_arr[i]=abs(slope_sparc-slope[i]);
			//cout << slope[i] << " slope " << slope_sparc << " target " << a_guess[i] << "a guess "<< endl;

		}
		a_temp = a_guess[argmin(slope_sparc_arr,guess_a)];
		if (argmin(slope_sparc_arr,guess_a) == 0) { a_temp = a_guess[1]; } 
		if (a_temp < 0) { a_temp = 0; }
	}
	// cout << slope_sparc << "poo poo" << slope[argmin(slope_sparc_arr,guess_a)] << endl;
	return {Mag_guess,a_temp,slope_sparc};
}

double make_vrot(double radi,double Mag,double Ropt,double alpha){
	// Returns a polyex rotation curve
	// Input a magnitude, Ropt, radii, and alpha
	// Catinella et al 2007
	// https://arxiv.org/abs/astro-ph/0512051
	double vt = V0_func(Mag);
	double rt = Ropt * rt_func(Mag);
	double vrot;
	double a = alpha;

	vrot = vt*(1.0 - exp(-radi/rt))*(1.0 + a*radi/rt);

	return vrot;}


struct Galaxy_params {
        double r_MHI, r_DHI, r_Mstar, r_Ropt, r_vflat, r_alpha;
        double r_Mag, r_slope, r_dist, r_beams;} ;

Galaxy_params setup_relations(double mass,double beams, double beam, double ring_thickness, bool scatter) {

	double MHI = pow(10.0,mass);
	double DHI = DHI_calc(MHI,scatter) ;
	double Mstar = Mstar_calc(MHI,scatter);
	//
	double vflat = BTFR(Mstar + 1.4*MHI,scatter);
	double Rs = (DHI/2.0) * 0.18;
	double Ropt = Ropt_calc(vflat,scatter);
	Mag_params Mag_stuff = Mag_calc(vflat,Ropt,DHI/2.0,Mstar,scatter);
	double Mag = Mag_stuff.Mag;
	//
	double alpha = Mag_stuff.alpha;
	double slope  = Mag_stuff.slope;
	double dist = DHI * (206265./(beam*beams));
	double delta = arcsec_to_rad(ring_thickness)*dist;
	
	return {MHI, DHI, Mstar, Ropt, vflat, alpha, Mag,slope,dist,beams};

}

