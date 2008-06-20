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
const int n_x=12,n_y=12,n_z=12;

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
	wall_plane pl(1,-1,1,0.6);con.add_wall(pl);
	wall_plane pl2(1,1,-1,0.6);con.add_wall(pl2);
	wall_plane pl3(-1,-1,-1,0.6);con.add_wall(pl3);
//	wall_plane pl4(-1,1,1,0.6);con.add_wall(pl4);

	//Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (con.point_inside(x,y,z)) {
			con.put(i,x,y,z);
			i++;
		}
	}
	con.guess_length_scale();

	con.store_cell_volumes(bb);x=0;
	for(i=0;i<particles;i++) x+=bb[i];cout << x << endl;

	// Save the Voronoi network of all the particles to a text file
	// in a format ready for plotting by gnuplot
	con.draw_gnuplot("voronoi_cells");
}
