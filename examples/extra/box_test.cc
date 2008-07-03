// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : June 1st 2008

#include "cell.cc"

// Parameters controlling the center of the test box
const double cx=1.5,cy=1.5,cz=1.5;

int main() {
	double x,y,z;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 20, centered
	// on the origin
	v.init(-10,10,-10,10,-10,10);

	// Cut by a grid of points in a box of width one, centered on
	// (cx,cy,cz)
	for(x=cx-0.45;x<cx+0.5;x+=0.1) {
		for(y=cy-0.45;y<cy+0.5;y+=0.1) {
			for(z=cz-0.45;z<cz+0.5;z+=0.1) {
				v.plane(x,y,z);
			}
		}
	}

	// Output the Voronoi cell in gnuplot format
	v.draw_gnuplot("box_cut.gnu",0,0,0);
}
