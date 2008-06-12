// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.cc"

const double tolwidth=1e-7;

int main() {
	int gp;
	double x,y,z,rsq,r,rmin,rmax;
	voronoicell v;
	ofstream os;
	os.open("test",ofstream::out|ofstream::trunc);
	
	// Initialize the Voronoi cell to be a cube of side length 2, centered on
	// the origin
//	v.init(-1,1,-1,1,-1,1);
	v.init_octahedron(1);

	// Cut the cell by 250 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
	for(int i=0;i<10000;i++) {
		x=double(2*rand()-1)/RAND_MAX;
		y=double(2*rand()-1)/RAND_MAX;
		z=double(2*rand()-1)/RAND_MAX;
		rsq=x*x+y*y+z*z;
		if(rsq>0.01&&rsq<1) {
			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
			rmin=0;rmax=1;
			while (v.plane_intersects_guess(x,y,z,rmax,gp)) rmax*=2;
			while (rmax-rmin>tolwidth) {
				r=(rmax+rmin)*0.5;
				if (v.plane_intersects_guess(x,y,z,r,gp)) rmin=r;
				else rmax=r;
			}
			r=(rmax+rmin)*0.5;
			x*=r;y*=r;z*=r;
			os << x << " " << y << " " << z << endl;
		}
	}
	
	os.close();
}
