// Superellipsoid example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	double x,y,z,rsq,r;
	voronoicell v;
	
	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,-1,1);

	// Cut the cell by 1200 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
	for(int i=0;i<2500;i++) {
		x=2*rnd()-1;
		y=2*rnd()-1;
		z=2*rnd()-1;
		rsq=x*x*x*x+y*y*y*y+z*z*z*z;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(sqrt(rsq));
			x*=r;y*=r;z*=r;
			v.plane(x*x*x,y*y*y,z*z*z,x*x*x*x+y*y*y*y+z*z*z*z);
		}
	}
	
	// Output the Voronoi cell to a file, in the gnuplot format
	v.draw_gnuplot("superellipsoid.gnu",0,0,0);
}
