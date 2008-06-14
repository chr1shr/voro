// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.cc"

const double pi=3.1415926535897932384626433832795;
const double tolwidth=1e-7;
const double theta_step=pi/200;

int main() {
	int gp=0;
	double x,y,z,rsq,r,rmin,rmax;
	double theta,phi,phi_step;
	voronoicell v;
	ofstream os;
	os.open("test",ofstream::out|ofstream::trunc);
	
	// Initialize the Voronoi cell to be a cube of side length 2, centered on
	// the origin
//	v.init(-1,1,-1,1,-1,1);
	v.init_octahedron(1);
	v.plane(0.4,0.3,1,0.1);
//	v.init_tetrahedron(-3,0,-1,0,3,-1,3,0,-1,0,1,4);
//	v.check_relations();

	// Cut the cell by 250 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
	for(theta=theta_step*0.5;theta<pi;theta+=theta_step) {
		phi_step=2*pi/(2*pi*int(sin(theta)/theta_step)+1);
		for(phi=phi_step*0.5;phi<2*pi;phi+=phi_step) {
			x=sin(theta)*cos(phi);
			y=sin(theta)*sin(phi);
			z=cos(theta);
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
	
	v.dump_gnuplot("test2",0,0,0);

	os.close();
}
