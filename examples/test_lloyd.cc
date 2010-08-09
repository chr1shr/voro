#include "cell_2d.cc"
#include "container_2d.cc"

int main() {
	int i;
	char buffer[64];

	// Initialize the container class to be the unit square, with
	// non-periodic boundary conditions. Divide it into a 6 by 6 grid, with
	// an initial memory allocation of 16 particles per grid square.
	container_2d con(0,1,0,1,6,6,false,false,16);
	
	// Import the spiral data set, and only save those particles that are
	// within the container bounds 
	con.import("particles_spiral");
	sprintf(buffer,"particles_spiral.%d",0);
	con.draw_particles(buffer);

	// Carry out sixty four iterations of Lloyd's algorithm
	for(i=0;i<64;i++) {
		con.clear();
		con.import(buffer);
		sprintf(buffer,"particles_spiral.%d",i+1);
		con.print_all_custom("%i %C",buffer);
	}
}
