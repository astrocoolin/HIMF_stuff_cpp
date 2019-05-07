#pragma once
#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include <random>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include "NCC.h"

using namespace std;

double uniform_random() {
        //construct a trivial random generator engine from a time-based seed
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0,1.0);
        return distribution(generator);
}

double HIMF_phi(double mass, double mstar, double alpha, double phistar) {
	return log(10.0) *phistar* pow((mass/mstar),(alpha+1.)) * exp(-mass/mstar);
}

int main() {
	int i;
	double z = 0.03;
	double c = 3E5; //km/s
	double H = 70;  //km/s/mpc
	bool keep = true;

	double V ;
	double D ;
	double N ;
	double prob ; 
	double mass ;
	ofstream myfile;
	myfile.open("Glist.txt",ios::trunc | ios::out);
	myfile << "index "  << "MHI "<< "DHI " << "Mstar "<< "vflat " << "alpha "  <<  "Mag " << "dist_MPC "  << "beams" << endl;
	int max ; 
	// cout << "Number of Galaxies?: " ;
	// cin >> max ; 
	#pragma omp parallel
	//for (i=0; i < 125321; i++){
	for (i=0; i < 1006971; i++){
	//for (i=0; i < max; i++){
		keep = true;
		while (keep) {
			V = pow(((c/H) * (pow((1.0+z),2.0)-1.0)/(pow((1.0+z),2.0)+1.0)),3.0) ;
			D = pow(uniform_random() * V , 1.0/3.0)*1000.0 ;

			N = 0.11*uniform_random();
			mass = 4.0 * uniform_random()+7;

			prob = HIMF_phi(pow(10.0,mass),pow(10.0,9.96),-1.33,4.8E-3) ; 
			//cout << prob << endl;
			if (prob > N){
				keep = false;
			}
		}
		Galaxy one;
		one.reroll(mass,10.0,false);
		one.calc_dist(D);
		#pragma omp critical

		myfile << i << " "  << log10(one.MHI) << " "<< one.DHI <<" " << log10(one.Mstar) << " "<< " "<< one.vflat << " "<< one.alpha  <<  " "<< one.Mag << " "<< one.dist/1000.  << " "<< one.beams << " " << endl;
	}
	myfile.close();
	
        return 0;
}

