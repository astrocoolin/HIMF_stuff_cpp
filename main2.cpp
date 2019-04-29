#pragma once
#include <iostream>
#include <fstream>
#include </home/colin/PhD/Codes/HIMF_stuff_cpp/NCC.cpp>


using namespace std;

int main() {
	int i;
	Galaxy one;
	#pragma omp parallel
	for (i=0; i < 125321; i++){
		//Galaxy * one = new Galaxy();
		one.reroll(7.5,1.0,true);
		if (i % 10 == 0){
			cout << log10(one.MHI)<< " "<< one.DHI <<" " << log10(one.Mstar) << endl;
		}	
	}
	
        return 0;
}

