// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "voro++.cc"

// Set up constants for the container geometry
const fpoint x_min=-3,x_max=3;
const fpoint y_min=-3,y_max=3;
const fpoint z_min=0,z_max=6;

// Set up the number of blocks that the container is divided
// into.
const int n_x=3,n_y=3,n_z=3;

int main() {
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block. Import
	// the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot format.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	con.import("pack_six_cube");
	con.draw_cells_gnuplot("pack_six_cube.gnu");

	// Create a container for polydisperse particles using the same
	// geometry as above. Import the polydisperse test packing and
	// output the Voronoi radical tessellation in gnuplot format.
	container_poly conp(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	conp.import("pack_six_cube_poly");
	conp.draw_cells_gnuplot("pack_six_cube_poly.gnu");
}
