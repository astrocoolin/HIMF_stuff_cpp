#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
using namespace std;

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

int setup_relations(float mass,float beams) {
	double MHI = pow(10.0,9.5);
	double DHI = DHI_calc(MHI) ;
	double Mstar = Mstar_calc(MHI);
	double vflat = BTFR(Mstar + 1.4*MHI);
	double Rs = (DHI/2.0) * 0.18;
	double Ropt = Ropt_calc(vflat);

	cout << DHI; 
	cout << "\n";

	cout << Mstar; 
	cout << "\n";

	cout << vflat; 
	cout << "\n";

	cout << Ropt;
	cout << "\n";

	return 0;
}

