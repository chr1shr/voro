// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#define YEAST_ROUTINES
const int stickycut=1500;

#include "voro++.cc"
#include "dynamic.cc"

// Set up constants for the container geometry
const fpoint x_min=-20,x_max=20;
const fpoint y_min=-20,y_max=20;
const fpoint z_min=-20,z_max=20;

// Set the computational grid size
const int n_x=20,n_y=20,n_z=20;

// Set the number of particles that are going to be randomly introduced
const int particles=3000;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i=0,j;
	fpoint x,y,z;char q[256];

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			true,true,true,8);

	// Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if(con.count(x,y,z,1)==0) con.put(i++,x,y,z);
	}

	for(i=0;i<=800;i++) {
		cout << i << " " << con.packing_badness<cond_all>() << endl;
		sprintf(q,"output/%04d_p.pov",i);con.draw_yeast_pov(q);
		sprintf(q,"gzip -f -9 output/%04d_p.pov",i);system(q);
	//	sprintf(q,"output/%04d_v.pov",i);con.draw_cells_pov(q);
	//	sprintf(q,"gzip -f -9 output/%04d_v.pov",i);system(q);
		for(j=0;j<1000*i;j++) {
			con.move<velocity_brownian2>();con.stick(0.7);
		}
	}
}
