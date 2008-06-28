// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "cell.cc"
#include "container.cc"
#include "wall.cc"

// Set up constants for the container geometry
const fpoint x_min=-6,x_max=6;
const fpoint y_min=-6,y_max=6;
const fpoint z_min=0,z_max=12;

// Set up the number of blocks that the container is divided
// into.
const int n_x=35,n_y=39,n_z=60;

// This function returns a random fpoint between 0 and 1;
fpoint rnd() {return fpoint(rand())/RAND_MAX;}

int main() {
	int i;
	fpoint x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// 16 particles within each computational block
	container_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	con.import("six_cube_poly");
	con.put(1000,0.1,0.1,11,2);
	con.put(1000,0.1,-0.1,11,2);
	con.put(1000,-0.1,0.1,11,2);
	con.put(1000,-0.1,-0.1,11,2);
	con.put(1000,0,0.14,11,2);
	con.put(1000,0.14,0,11,2);
	con.put(1000,0,-0.14,11,2);
	con.put(1000,-0.14,0,11,2);
	wall_cylinder c1(0,0,0,0,0,1,4);con.add_wall(c1);

//	con.guess_length_scale();
	// Print out a list of the particles, and their Voronoi volumes
//	con.store_cell_volumes(bb);

	con.dump("diag");
	// Save the Voronoi network of all the particles to a text file
	// in a format ready for plotting by gnuplot
	con.draw_gnuplot("voronoi_cells5");
}
