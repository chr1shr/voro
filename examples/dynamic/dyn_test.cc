// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// Set up constants for the container geometry
const fpoint x_min=-4,x_max=4;
const fpoint y_min=-4,y_max=4;
const fpoint z_min=-4,z_max=4;

// Set the computational grid size
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=325;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i=0;
	fpoint x,y,z;char q[256];

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container
	wall_sphere sph(0,0,0,4);
	con.add_wall(sph);
	
	// Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (x*x+y*y+z*z<3*3) {con.put(i,x,y,z);i++;}
	}

	for(i=0;i<100;i++) {
		con.move<velocity_brownian>();
		con.full_relax(0.8);
		con.draw_particles_pov("temp.pov");
		con.draw_cells_pov("temp2.pov");
		sprintf(q,"povray +H400 +W400 +A0.3 -J -D +Ofr_%04d.png dyn_master.pov",i);
		system(q);
	}
}
