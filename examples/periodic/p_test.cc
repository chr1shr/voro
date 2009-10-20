// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// Set up the number of blocks that the container is divided into
const int n_x=4,n_y=4,n_z=4;
const fpoint pi=3.1415926535897932384626433832795;

// This function returns a random double between 0 and 1
fpoint rnd() {return fpoint(rand())/RAND_MAX;}

int main() {

	int i=0;
	fpoint r,theta,x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container_periodic con(10,2.5,10,0.4,0.2,10,n_x,n_y,1,8);

	for(r=0.5;r<5;r+=1) {
		for(theta=pi/100;theta<2*pi;theta+=pi/50) {
			x=5+r*cos(theta);
			y=5+r*sin(theta);
			z=5;
			con.put(i,x,y,z);
			i++;
		}
	}

//	con.print_network("ptest.network");
	con.draw_particles("ptest.par");
//	con.draw_cells_gnuplot("ptest.out");
}
