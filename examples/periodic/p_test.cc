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
	container_periodic con(7,-4.5,10,-3.4,2.7,13,n_x,n_y,n_z,8);

	for(r=0.5;r<5;r+=1) {
		for(theta=pi/(40*r);theta<2*pi;theta+=pi/(20*r)) {
			x=5+r*cos(theta);if(x>10||x<0) continue;
			y=5+r*sin(theta);if(y>10||y<0) continue;
			z=5+2*sin(2*theta);
			con.put(i,x,y,z);
			i++;
		}
	}

//	con.print_network("ptest.network");
	con.draw_particles("ptest.par");
//	con.draw_cells_gnuplot("ptest.out");
}
