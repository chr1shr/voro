#include "voro++_2d.hh"
#include <iostream>
	

int main() {
	cout << "beginning test";
	// Set up the geometry and import the test file
	container_2d con(-2,2,-2,2,4,4,false,false,false,16);
	con.import("/users/mac/voro/branches/2d_boundary/examples/bd_test");
	con.setup();

	con.debug_output();
	con.draw_cells_gnuplot("bd_test.gnu");	
	// Output the Voronoi cells to a file, in the gnuplot format.
	
	//con.draw_cells_gnuplot("bd_test.gnu");
}
