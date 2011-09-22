
#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double bx=10;
const double by=10;
const double bz=10;
const double bxy=0;
const double bxz=0;
const double byz=0;

// Set up the number of blocks that the container is divided
// into
const int n_x=3,n_y=3,n_z=3;

int main() {

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.	
        container_periodic con(bx,bxy,by,bxz,byz,bz,n_x,n_y,n_z,8);

	// Import the monodisperse test packing
        con.import("pack_ten_cube");

        // Output volume
        double vvol=con.sum_cell_volumes();
        printf("Voronoi volume   : %g\n",vvol);

        // Output particle ID
        con.print_custom("%i","pack_ten_cube_output.txt");

        // Output the Voronoi tessellation in POV-Ray format.
        con.draw_cells_pov("pack_ten_cube_v.pov");
        con.draw_particles_pov("pack_ten_cube_p.pov");

        // Output particle positions
        con.draw_particles("pack_ten_cube_traj.txt");
        con.draw_cells_gnuplot("pack_ten_cube_par.txt");

}
