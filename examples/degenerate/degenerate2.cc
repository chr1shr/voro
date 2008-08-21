// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

const double pi=3.1415926535897932384626433832795;
const int n=64;
const double theta=0.04;
const double step=2*pi/n;
const int planes=100;

int main() {
	double x,y,z,rsq,r,phi;
	voronoicell v;
	int n=0;

	// Initialize the Voronoi cell to be a cube of side length 2, centered on
	// the origin
	v.init(-1,1,-1,1,-1,1);

	// Plane cutting
	while(n<planes) {
		x=double(2*rand()-1)/RAND_MAX;
		y=double(2*rand()-1)/RAND_MAX;
		z=double(2*rand()-1)/RAND_MAX;
		rsq=x*x+y*y+z*z;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
			rsq=sqrt(x*x+y*y);r=z/rsq;
			for(phi=double(rand())/RAND_MAX*step;phi<2*pi;phi+=step)
				v.plane(x*cos(theta)+sin(theta)*(-y*cos(phi)/rsq-x*r*sin(phi)),
					y*cos(theta)+sin(theta)*(x*cos(phi)/rsq-y*r*sin(phi)),
					z*cos(theta)+sin(theta)*rsq*sin(phi),1);
			n++;
		}
	}
	

	// Output the Voronoi cell to a file in Gnuplot format
	v.draw_gnuplot("degenerate2.gnu",0,0,0);

	// Optional POV-Ray output
	v.draw_pov("degenerate2_v.pov",0,0,0);
	v.draw_pov_mesh("degenerate2_m.pov",0,0,0);
}
