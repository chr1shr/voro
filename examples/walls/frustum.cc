// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// Set the number of particles that are going to be randomly
// introduced
const int particles=200;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i=0;
	fpoint x,y,z;
	
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container con(-1.2,1.2,-1.2,1.2,0,1,14,14,7,
			false,false,false,8);

	// Add a cylindrical wall to the container
	wall_cone cone(0,0,2,0,0,-1,atan(0.5));
	con.add_wall(cone);

	// Add random particles inside t
/*	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (con.point_inside(x,y,z)) {
			con.put(i,x,y,z);i++;
		}
	}*/

	for(z=-0.95;z<1;z+=0.2) for(y=-0.95;y<1;y+=0.2) for(x=-0.95;x<1;x+=0.2) {
		if (con.point_inside(x,y,z)) {
			con.put(i,x,y,z);i++;
		}
	}

	// Output the particle positions in POV-Ray format
	con.draw_particles("frustum_p.pov");

	// Output the Voronoi cells in POV-Ray format
	con.draw_cells_gnuplot("frustum_v.pov");
}
