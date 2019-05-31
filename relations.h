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
	
	//double num ; 
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

double Mstar_calc_2(double MHI, bool scatter_flag) {
	// Huang et al 2012
	// HI Mass - Stellar Mass Relationship 
	// https://ui.adsabs.harvard.edu/abs/arXiv:1207.0523
	double in_slope[2] ;
	double in_constant[2] ; 
	double scatr[2] ;

	if (log10(MHI)  < 9.526) {
		in_slope[0] = 0.712 ;
		in_slope[1] = 0.0 ;
		in_constant[0] = 3.117 ;
		in_constant[1] = 0.0 ;
		scatr[0] = 0.0 ;
		scatr[1] = 0.0 ;
	} else {
		in_slope[0] = 0.276 ;
		in_slope[1] = 0.0 ;
		in_constant[0] = 7.042 ;
		in_constant[1] = 0.0 ;
		scatr[0] = 0.0 ;
		scatr[1] = 0.0 ;
	}

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);
	//double scatter = error_spread(scatr,scatter_flag);

	double Mstar = log10(MHI)/slope-constant/slope;
	//double spread[2] = {Mstar,scatter};

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

double BTFR_2_1(double Mbar, bool scatter_flag) {
	// Lelli et al 2016
	// Baryonic Tully-Fisher relationship
	// http://adsabs.harvard.edu/abs/2016ApJ...816L..14L
	
	double in_slope[2] = {3.75,0.11};
	double in_constant[2] = {2.18,0.23};

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);

	double logv = (log10(Mbar) - constant)/slope;

	return pow(10.0,logv);
}
double BTFR_2_2(double Mbar, bool scatter_flag) {
	// Lelli et al 2016
	// Baryonic Tully-Fisher relationship
	// http://adsabs.harvard.edu/abs/2016ApJ...816L..14L
	
	double in_slope[2] = {3.90,0.11};
	double in_constant[2] = {1.92,0.23};

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);

	double logv = (log10(Mbar) - constant)/slope;

	return pow(10.0,logv);
}
double BTFR_2_3(double Mbar, bool scatter_flag) {
	// Lelli et al 2016
	// Baryonic Tully-Fisher relationship
	// http://adsabs.harvard.edu/abs/2016ApJ...816L..14L
	
	double in_slope[2] = {3.71,0.11};
	double in_constant[2] = {2.27,0.23};

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);

	double logv = (log10(Mbar) - constant)/slope;

	return pow(10.0,logv);
}
double BTFR_2_4(double Mbar, bool scatter_flag) {
	// Lelli et al 2016
	// Baryonic Tully-Fisher relationship
	// http://adsabs.harvard.edu/abs/2016ApJ...816L..14L
	
	double in_slope[2] = {3.95,0.11};
	double in_constant[2] = {1.86,0.23};

	double constant = error_spread(in_constant,scatter_flag);
	double slope = error_spread(in_slope,scatter_flag);

	double logv = (log10(Mbar) - constant)/slope;

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
	//double spread[2] = {Ropt,scatr};
	//Ropt = error_spread(spread,scatter_flag);

	return h*pow(10.0,Ropt);
}

double Ropt_D_calc(double DHI, bool scatter_flag) {
	// Broeils & Rhee
	// Optical Radius (R25) - DHI relationship
	// https://ui.adsabs.harvard.edu/abs/1997A%26A...324..877B/abstract
	
	double spread[2] = {1.7,1.5};
	double scatr = error_spread(spread,scatter_flag);
	//Ropt = error_spread(spread,scatter_flag);
	double Ropt = DHI /1.7;

	return Ropt;
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
	double slope_temp;

	// Set slope from NIHAO 17
	double slope_sparc = 0.123 - 0.137*(log10(mstar)-9.471) ;
	double spread[2] = {slope_sparc,0.19};
	slope_sparc = error_spread(spread,scatter_flag);
	int jay= 0; 

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
		}
		
		for (int i=0; i<= Mag_length;i++) {vrot_compare[i] = abs(vrot - vt[i]); }
	
		// Best guess for Magnitude based on vrot with other params
		// Finds index of vt that most closely matches vrot and
		// that matches the Magnitude
		int ind = argmin(vrot_compare,Mag_length);
		Mag_guess = Mag[ind];
		double rt_guess = rt[ind];
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
			//cout << i << " "<< slope_sparc << " " <<slope[i] << endl;
		}
		a_temp = a_guess[argmin(slope_sparc_arr,guess_a)];
		jay = argmin(slope_sparc_arr,guess_a);
		slope_temp = slope[jay];

		//cout << jay << endl;

		if (argmin(slope_sparc_arr,guess_a) == 0) { a_temp = a_guess[1]; } 
		if (a_temp < 0) { a_temp = 0; }
	}
	//cout <<RHI<<" " << jay << " JaY " << slope_temp << " " << slope_sparc << " " << a_temp << endl;
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

	struct my_f_params {double RHI; double xdx; double vt; double Rs;};

double f (double x, void * p) {
	struct my_f_params * params = (struct my_f_params *)p;
	double RHI = (params->RHI);
	double xdx = (params->xdx);
	double vt  = (params->vt);
	double Rs  = (params->Rs);

	double f1 = exp(- pow((x-0.4*RHI)/(sqrt(2.0)*(xdx+0.36)*RHI),2));
	double f2 = (sqrt(vt/120.0)-1)*exp(-x/Rs);

	//double f1_RHI = exp(- pow((RHI-0.4*RHI)/(sqrt(2.0)*(xdx+0.36)*RHI),2));
	//double f2_RHI = (sqrt(vt/120.0)-1)*exp(-RHI/Rs);


	double f = (f1-f2)*x*2.0*3.1415926535;
	return f;
}

double sbr_func (double x, double RHI, double xdx, double vt, double Rs) {
	//struct my_f_params * params = (struct my_f_params *)p;
	//double RHI = (params->RHI);
	//double xdx = (params->xdx);
	//double vt  = (params->vt);
	//double Rs  = (params->Rs);

	double f1 = exp(- pow((x-  0.4*RHI)/(sqrt(2.0)*(xdx+0.36)*RHI),2));
	double f2 = (sqrt(vt/120.0)-1)*exp(-x  /Rs);

	double f = (f1-f2);

	return f;
}

double sbr_func_f (double x, void * p) {
	struct my_f_params * params = (struct my_f_params *)p;
	double RHI = (params->RHI);
	double xdx = (params->xdx);
	double vt  = (params->vt);
	double Rs  = (params->Rs);

	double f1 = exp(- pow((x-  0.4*RHI)/(sqrt(2.0)*(xdx+0.36)*RHI),2));
	double f2 = (sqrt(vt/120.0)-1)*exp(-x  /Rs);

	double f = (f1-f2)/sbr_func(RHI,RHI,xdx,vt,Rs) - 0.5;

	return f;
}

double Match_velocity (double x) {
	double a = -7.28804304e-04;
	double b = 4.11594154e-02;
	double c = -9.65306151e-01;
	double d = 1.20646885e+01;
	double e = -8.49542020e+01;
	double f = 3.20597179e+02;
	double g = -5.06959891e+02;

	return a*pow(x,6) + b*pow(x,5) + c*pow(x,4) + d*pow(x,3) + e*pow(x,2) + f*x + g;

}

double find_halfmass (double RHI,double xdx, double vt, double Rs){
	int i, status;
	gsl_root_fsolver *workspace_f;
	double x, x_l, x_r;

	workspace_f = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
	my_f_params beta = {RHI,xdx,vt,Rs};
	
	gsl_function F;
	F.function = &sbr_func_f;
	F.params = &beta;
	x_l = RHI;
	x_r = RHI*2.0;

	gsl_root_fsolver_set(workspace_f, &F, x_l, x_r);

	for(i = 0; i < 15; i++){
		status = gsl_root_fsolver_iterate(workspace_f);

		x_l = gsl_root_fsolver_x_lower(workspace_f);
        	x_r = gsl_root_fsolver_x_upper(workspace_f);

        	status = gsl_root_test_interval(x_l, x_r, 1.0e-13, 1.0e-20);
        	if(status != GSL_CONTINUE){break; }
	}
	gsl_root_fsolver_free(workspace_f);
	x = (x_r+x_l)/2.0;
	//cout << x_r << " " << x_l << endl;
	return x;
}

double integrate_inf(double RHI,double xdx, double vt, double Rs){
	gsl_integration_workspace * w
		= gsl_integration_workspace_alloc (1000);

	double result, error;
	my_f_params alpha = {RHI,xdx,vt,Rs};

	gsl_function F;
	F.function = &f;
	F.params = &alpha;
 	gsl_integration_qagiu (&F, 0, 1e-12, 1e-6, 1000,
		w, &result, &error);
	//cout << result << " Result" << endl;
	gsl_integration_workspace_free(w);
	return result;
}


struct sbr_params { double x_dx; double Mass_guess;};

sbr_params sbr_calc(double RHI,double vt,double Rs,int radi_len,double mass){
	int dx_range = (0.075+0.075)/0.0005;
	double delta[dx_range+1];
	double Mass_guess[dx_range+1];
	double Mass_compare[dx_range+1];
	for (int i=0;i<dx_range+1;i++){
		delta[i] = i*0.0005 - 0.075;
		Mass_guess[i] = log10(integrate_inf(RHI,delta[i],vt,Rs)
				*1000.0*1000.0/(sbr_func(RHI,RHI,delta[i],vt,Rs)));
		Mass_compare[i] = abs(Mass_guess[i] - mass) ;
		//cout << i <<" " << delta[i] << " " << Mass_guess[i] << " "<< Mass_compare[i] << " " << mass << " " << integrate_inf(RHI,delta[i],vt,Rs) << " " << sbr_func(RHI,RHI,delta[i],vt,Rs) << endl;
	}
	int j = argmin(Mass_compare,dx_range+1);
	return {delta[j],Mass_guess[j]};
}


struct Galaxy_params {
        double r_MHI, r_DHI, r_Mstar, r_Ropt, r_vflat, r_alpha;
        double r_dx, r_Mag, r_RHI5, r_slope, r_dist, r_beams;} ;

Galaxy_params setup_relations(double mass,double beams, double beam, double ring_thickness, bool scatter) {

	double MHI = pow(10.0,mass);
	double DHI = DHI_calc(MHI,scatter) ;
	double Mstar = Mstar_calc_2(MHI,scatter);
	//cout << log10(Mstar) << " "<< log10(MHI)<< endl;
	double vflat = pow(10,(Match_velocity(log(MHI))));
	//double vflat = BTFR_2_1(Mstar + 1.4*MHI,scatter);
	double Rs = (DHI/2.0) * 0.18;
	double Ropt = Ropt_D_calc(DHI/2.,scatter);
	Mag_params Mag_stuff = Mag_calc(vflat,Ropt,DHI/2.0,Mstar,scatter);

	double Mag = Mag_stuff.Mag;
	//
	double alpha = Mag_stuff.alpha;
	double slope  = Mag_stuff.slope;
	double dist = DHI * (206265./(beam*beams));
	double delta = arcsec_to_rad(ring_thickness)*dist;
	int radi_len = (DHI+delta) /delta;

	sbr_params sbr_stuff = sbr_calc(DHI/2.0,vflat,Rs,radi_len,mass);
	double RHI5 =  find_halfmass(DHI/2.0,sbr_stuff.x_dx,vflat,Rs);

	// my_f_params beta = {DHI/2.0,sbr_stuff.x_dx,vflat,Rs};
	// double delta = arcsec_to_rad(ring_thickness)*dist;
	
	return {MHI, DHI, Mstar, Ropt, vflat, alpha, sbr_stuff.x_dx, Mag, RHI5,slope,dist,beams};

}

