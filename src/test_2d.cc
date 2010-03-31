// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "cell2d.cc"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	double x,y,rsq,r;
	voronoicell_2d v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1);

	v.plane(1,0.5,1);
	printf("---\n");
//	v.plane(-0.9,0.9,1.62);
//	printf("---\n");
//	v.plane(-0.9,-0.9,1.62);
//	printf("---\n");
//	v.plane(1,-1,2);

	// Cut the cell by 250 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
/*	for(int i=0;i<250;i++) {
		x=2*rnd()-1;
		y=2*rnd()-1;
		rsq=x*x+y*y;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(rsq);x*=r;y*=r;
			v.plane(x,y,1);
		}
	}*/

	// Output the Voronoi cell to a file, in the gnuplot format
	v.draw_gnuplot("single_cell.gnu",0,0);
}
