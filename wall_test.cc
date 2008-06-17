// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "cell.cc"
#include "container.cc"
#include "wall.cc"

// Set up constants for the container geometry
const fpoint x_min=-1,x_max=1;
const fpoint y_min=-1,y_max=1;
const fpoint z_min=-1,z_max=1;

// Set up the number of blocks that the container is divided
// into.
const int n_x=8,n_y=8,n_z=8;

// Set the number of particles that are going to be randomly
// introduced
const int particles=1000;

// This function returns a random fpoint between 0 and 1;
fpoint rnd() {return fpoint(rand())/RAND_MAX;}

int main() {
	int i=0;
	fpoint *bb;
	bb=new fpoint[particles];
	fpoint x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// 16 particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	wall_sphere sph(0,0,0,1);con.add_wall(sph);

	//Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (sph.point_inside(x,y,z)) {
			con.put(i,x,y,z);
			i++;
		}
	}

	// Save the Voronoi network of all the particles to a text file
	// in a format ready for plotting by gnuplot
	con.draw_gnuplot("voronoi_cells");
}
