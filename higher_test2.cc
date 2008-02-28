// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "container.cc"

const double pi=3.1415926535897932384626433832795;
const int n=64;
const double theta=0.2;
const double step=2*pi/n;

int main() {
	double x,y,z,rsq,r,phi;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered on
	// the origin
	v.init(-1,1,-1,1,-1,1);

	// Output the initial cell	
	ofstream file;
	file.open("intest",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0.0,0.0,0.0);
	file.close();
	
	// Plane cutting
	for(int n=0;n<500;n++) {
		x=double(2*rand()-1)/RAND_MAX;
		y=double(2*rand()-1)/RAND_MAX;
		z=double(2*rand()-1)/RAND_MAX;
		rsq=x*x+y*y+z*z;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
			rsq=sqrt(x*x+y*y);r=z/rsq;
			for(phi=double(rand())/RAND_MAX*step;phi<2*pi;phi+=step)
				v.plane(x*cos(theta)+sin(theta)*(-y*cos(phi)-x*r*sin(phi)),
					y*cos(theta)+sin(theta)*(x*cos(phi)-y*r*sin(phi)),
					z*cos(theta)+sin(theta)*rsq*sin(phi),1);
//			(x,y,z);
//			(-y,x,0);
//			(-x*z/sqrt,-y*z/sqrt,sqrt(x*x+y*y))
//			v.plane(x,y,z,1);
		}
	}

	// Output the Voronoi cell to a file, in the gnuplot format
	file.open("coolm.pov",ofstream::out|ofstream::trunc);
	v.dumppovmesh(file,0.0,0.0,0.0);
	file.close();
	
	file.open("coolp.pov",ofstream::out|ofstream::trunc);
	v.dumppov(file,0.0,0.0,0.0);
	file.close();
}
