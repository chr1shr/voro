// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "cell.cc"
#include "container.cc"
#include "wall.cc"

// Set up constants for the container geometry
const fpoint x_min=-6.5,x_max=6.5;
const fpoint y_min=-6.5,y_max=6.5;
const fpoint z_min=0,z_max=18.5;

// Set up the number of blocks that the container is divided
// into
const int n_x=7,n_y=7,n_z=14;

// Set the number of particles that are going to be randomly
// introduced
const int particles=100000;

// This function returns a random fpoint between 0 and 1
fpoint rnd() {return fpoint(rand())/RAND_MAX;}

int main() {
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// 8 particles within each computational block.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container
	wall_cylinder cyl(0,0,0,0,0,1,6);
	con.add_wall(cyl);

	// Import the particles from a file
	con.import("pack_small_cylinder");

	// Save the Voronoi network of all the particles to a text file
	// in a format ready for plotting by gnuplot
	con.draw_pov("voronoi_cells.pov");

	con.dump_pov("particles.pov");
}
