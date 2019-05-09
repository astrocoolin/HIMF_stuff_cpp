#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>
#include <random>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include "NCC.h"

using namespace std;

double uniform_random(double number) {
        //construct a true random number
	random_device rd;
	mt19937 gen(rd());
        uniform_real_distribution<double> dis(0.0,number);
        return dis(gen);
}

double HIMF_phi(double mass, double mstar, double alpha, double phistar) {
	mass = pow(10.0,mass);
	mstar = pow(10.0,mstar);
	double frac = mass / mstar;
	double answer = log(10.0) *phistar* pow(frac,alpha+1.) * exp(-frac);
	return answer ;
}

int main() {
	int i;
	double z = 0.03;
	double c = 3E5; //km/s
	double H = 70;  //km/s/mpc

	ofstream myfile;
	myfile.open("Glist.txt",ios::trunc | ios::out);
	myfile << "MHI "<< "DHI " << "Mstar "<< "vflat " << "alpha "  <<  "Mag " << "dist_MPC "  << "beams" << endl;
	//for (i=0; i < 125321; i++){
	#pragma omp parallel num_threads(4)
	#pragma omp for
	//for (i=0; i < 1006971; i++){
	for (i=0; i < 456971934; i++){
	//for (i=0; i < 10000; i++){
		bool keep = true;
		double V;
		double D;
		double N;
		double mass;
		double prob;
		while (keep) {
			V = pow(((c/H) * (pow((1.0+z),2.0)-1.0)/(pow((1.0+z),2.0)+1.0)),3.0) ;
			D = pow(uniform_random(V) , 1.0/3.0)*1000.0 ;

			N = uniform_random(0.11);
			mass = uniform_random(4.0)+7.0;

			prob = HIMF_phi(mass,9.96,-1.33,4.8E-3) ; 
			//cout << prob << endl;
			if (N < prob){
				keep = false;

				Galaxy one;
				one.reroll(mass,10.0,false);
				one.calc_dist(D);

				#pragma omp critical
				{
				myfile << log10(one.MHI) << " "<< one.DHI <<" " << 
					log10(one.Mstar) << " "<< " "<< one.vflat << 
					" "<< one.alpha  <<  " "<< one.Mag << " "<< 
					one.dist/1000.  << " "<< one.beams << " " << endl;
				}
			}
		
		}

	}
	myfile.close();
	
        return 0;
}

