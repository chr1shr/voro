// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "cell.cc"

const double pi=3.1415926535897932384626433832795;
const int n=64;
const double theta=0.3;
const double step=2*pi/n;

int main() {
	double x,y,z,rsq,r,phi;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered on
	// the origin
	v.init(-1,1,-1,1,-1,1);

	// Output the initial cell	
	v.dump_gnuplot("intest",0,0,0);
	
	// Plane cutting
	for(int n=0;n<350;n++) {
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
		}
	}
	

	// Output the Voronoi cell to a file, in the gnuplot format
	v.dump_gnuplot("high_output",0,0,0);

	// Optional POV output
//	v.dump_povmesh("high_output_mesh.pov",0,0,0);
//	v.dump_pov("high_output.pov",0,0,0);
}
