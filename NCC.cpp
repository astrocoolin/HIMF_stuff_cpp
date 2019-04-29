#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <random>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include </home/colin/PhD/Codes/HIMF_stuff_cpp/relations.cpp>


using namespace std;
/*
struct Galaxy {
	Galaxy();
	~Galaxy();
}
*/
class Galaxy {
	public:
		double MHI;
		double DHI;
		double Mstar;
		double Ropt;
		double vflat;
		double alpha;
		double dx;
		double Mag;
		double slope;
		double dist;
		double beams;
		
		void reroll(double mass,double newbeams,bool scatter) {
			beams = newbeams;
			double beamsize = 30;
			double rwidth = 2.5;


			Galaxy_params rerolled_params = setup_relations(mass, newbeams, beamsize, rwidth,scatter);
                	MHI = rerolled_params.r_MHI;
                	DHI = rerolled_params.r_DHI;
                	Mstar = rerolled_params.r_Mstar;
                	Ropt = rerolled_params.r_Ropt;
                	vflat = rerolled_params.r_vflat;
                	alpha = rerolled_params.r_alpha;
                	dx = rerolled_params.r_dx;
                	Mag = rerolled_params.r_Mag;
                	slope = rerolled_params.r_slope;
                	dist = rerolled_params.r_dist;
			beams = rerolled_params.r_beams;
			// cout << log10(MHI) << " poop" << endl;
		}

		void calc_disc(){}
	
};

