// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// Set up the number of blocks that the container is divided into
const int n_x=4,n_y=4,n_z=4;

int main() {

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container_periodic con(6.926,0,6.926,0,0,6.41,n_x,n_y,n_z,8);

	cout.precision(15);

	//Randomly add particles into the container
	con.import("EDI");

	// Save the Voronoi network of all the particles to text files
	// in gnuplot and POV-Ray formats
//	con.print_network("EDI.network");
	
	con.draw_particles("EDI.par");
	con.draw_cells_gnuplot("EDI.out");
}
