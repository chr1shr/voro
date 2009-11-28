// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,-1,1);

	v.plane(1,1,0,0.5);

	voronoicell w(v);

	v.plane(1,1,1,0.25);
	w.plane(1,1,-1,0.25);

	v.init(-1,1,-1,1,-1,1);
	w.init(v);
	v.plane(1,1,1,0.25);

	// Output the Voronoi cell to a file, in the gnuplot format
	v.draw_gnuplot("v_cell.gnu",0,0,0);
	w.draw_gnuplot("w_cell.gnu",0,0,0);
}
