// Box cutting example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : June 1st 2008

#include "cell.hh"

// Parameters controlling the center of the test box
const double cx=1.5,cy=1.5,cz=1.5;

int main() {
	double x,y,z;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 16, centered
	// on the origin
	v.init(-8,8,-8,8,-8,8);

	// Cut by a grid of points in a box of width one, centered on
	// (cx,cy,cz)
	for(x=cx-0.5;x<cx+0.55;x+=0.1) for(y=cy-0.5;y<cy+0.55;y+=0.1)
		for(z=cz-0.5;z<cz+0.55;z+=0.1) v.plane(x,y,z);

	// Output the Voronoi cell in gnuplot format
	v.draw_gnuplot(0,0,0,"box_cut.gnu");

	// Now make a small file that contains the test box
	v.init(cx-0.5,cx+0.5,cy-0.5,cy+0.5,cz-0.5,cz+0.5);
	v.draw_gnuplot(0,0,0,"box_cut.points");
}
