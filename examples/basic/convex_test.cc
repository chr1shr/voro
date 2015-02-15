// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (Harvard SEAS / LBL)
// Email    : chr@alum.mit.edu
// Date     : February 16th 2014

#include "voro++.hh"
using namespace voro;

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
//	double x,y,z,rsq,r;
	voronoicell v;

	v.init_l_shape();
	v.nplane_new(-1,3,0,10,0);

	v.draw_gnuplot(0,0,0,"single_cell.gnu");
}
