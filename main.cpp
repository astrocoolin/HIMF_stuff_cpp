#include <iostream>
#include <fstream>

using namespace std;
int setup_relations(double mass,double beams, double beam, double ring_thickness, double inclination);
int emptyFits();

int main() {
	setup_relations(9.5,4,30,2.5,60.0);
	emptyFits();
	system("tirific deffile=cube_input.def");
	return 0;
}

