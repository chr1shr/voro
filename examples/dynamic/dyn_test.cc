// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// Set up constants for the container geometry
const fpoint x_min=-6,x_max=6;
const fpoint y_min=-6,y_max=6;
const fpoint z_min=-6,z_max=6;

// Set the computational grid size
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=2000;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;
	fpoint x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container
	wall_cylinder sph(0,0,0,0,0,1,3);
	con.add_wall(sph);
	
	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (con.point_inside(x,y,z)) con.put(i,x,y,z);
	}

	// Output the particle positions in POV-Ray format
	con.wall_diagnostic();

}
