// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// Set up constants for the container geometry
const fpoint x_min=-5,x_max=5;
const fpoint y_min=-5,y_max=5;
const fpoint z_min=-0.5,z_max=0.5;

// Set the computational grid size
const int n_x=8,n_y=8,n_z=1;

// Set the number of particles that are going to be randomly introduced
const int particles=100;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

void output_all(container_dynamic &con,int i) {
	char q[256];
	sprintf(q,"output/%04d_p.pov",i);con.draw_particles_pov(q);
	sprintf(q,"gzip -f -9 output/%04d_p.pov",i);system(q);
	sprintf(q,"output/%04d_v.pov",i);con.draw_cells_pov(q);
	sprintf(q,"gzip -f -9 output/%04d_v.pov",i);system(q);
}

int main() {
	int i=0;
	fpoint x,y;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container
//	wall_sphere sph(0,0,0,7);
//	con.add_wall(sph);
	
	// Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		con.put(i,x,y,0);
		i++;
	}

	output_all(con,0);
	for(i=0;i<=10;i++) con.full_relax(0.8);
	output_all(con,1);

	for(i=20;i<=1000;i++) {
		con.spot(0,0,0,0.005,0.005,0,2.5);
	//	con.relax(4,4,0,3.5,0.8);
		con.full_relax(0.8);
		if(i%10==0) {
			output_all(con,i/10);
		}
	}
}
