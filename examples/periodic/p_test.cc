// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// Set up the number of blocks that the container is divided into
const int n_x=4,n_y=4,n_z=4;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {

	int i;
	fpoint x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container_periodic con(1,0,1,0,0,1,n_x,n_y,n_z,8);

	for(i=0;i<20;i++) {
		x=rnd();
		y=rnd();
		z=rnd();
		con.put(i,x,y,z);
	}

	con.print_network("ptest.network");
	con.draw_particles("ptest.par");
	con.draw_cells_gnuplot("ptest.out");
}
