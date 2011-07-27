#include "voro++_2d.hh"

int main() {

	// Set up the geometry and import the test file
	container_2d con(-2.5,2.5,-2.5,2.5,4,4,true,true,false,16);
	con.import("bd_test");

	// Output the boundaries to a file 
	con.draw_boundary("bd_test.out");

	// Output the Voronoi cells to a file, in the gnuplot format.
	//con.draw_cells_gnuplot("bd_test.gnu");
}
