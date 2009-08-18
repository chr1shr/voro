// Custom output example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 10th 2009

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
	// eight particles within each computational block.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Import the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	con.import("pack_six_cube");

	// Use the basic output routine, that saves the particle IDs,
	// positions, and Voronoi cell volumes
	con.print_all("packing.standard");

	// Use the neighbor output routine, that also includes information
	// about what particles share a face. The numbers -1 to -6 in the
	// neighbor list correspond to faces that are in contact with the
	// container walls.
	con.print_all_neighbor("packing.neighbor");

	// Do a custom output routine to store the number of vertices, edges,
	// and faces of each Voronoi cell
	con.print_all_custom(
		"ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s",
		"packing.custom1");

	// Do a custom output routine to store a variety of face-based
	// statistics. Store the particle ID and position, the number of faces
	// the total face area, the order of each face, the areas of each face,
	// the vertices making up each face, and the neighboring particle (or
	// wall) corresponding to each face.
	con.print_all_custom("%i %q %s %F %a %f %t %l %n","packing.custom2");

	// Do a custom output routine that outputs the particle IDs and
	// positions, plus the volume and the centroid position relative to the
	// particle center
	con.print_all_custom("%i %q %v %c","packing.custom3");

	// Also create POV-Ray output of the Voronoi cells for use in the
	// rendering
	con.draw_cells_pov("pack_six_cube_v.pov");
}
