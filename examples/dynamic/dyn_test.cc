// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// Set up constants for the container geometry
const fpoint x_min=-8,x_max=8;
const fpoint y_min=-8,y_max=8;
const fpoint z_min=-8,z_max=8;

// Set the computational grid size
const int n_x=8,n_y=8,n_z=8;

// Set the number of particles that are going to be randomly introduced
const int particles=1750;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

class velocity_twist {
	public:
		velocity_twist() : track_ve(false), ang(0.005) {}; 
		inline void vel(int ijk,int q,fpoint &x,fpoint &y,fpoint &z) {
			double px=x;
			x=x*cos(ang*z)+y*sin(ang*z);
			y=y*cos(ang*z)-px*sin(ang*z);
		}
		const bool track_ve;
	private:
		const fpoint ang;
};

int main() {
	int i=0,j,k;
	fpoint x,y,z;//char q[256];
	int *u[201];

	for(j=0;j<=200;j++) {
		u[j]=new int[400];
		for(k=0;k<400;k++) u[j][k]=0;
	}

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container
	wall_sphere sph(0,0,0,7);
	con.add_wall(sph);
	
	// Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (x*x+y*y+z*z<4*4) {con.put(i,x,y,z);i++;}
	}

	for(i=0;i<=200;i++) {
		con.neighbor_distribution<cond_all>(u[i],0.02,400);
		con.move<velocity_twist>();
		con.full_relax(1.6);
		con.full_relax(1.6);
		con.full_relax(1.6);
		con.full_relax(1.2);
		con.full_relax(1.2);
		con.full_relax(1.0);
		con.full_relax(0.8);
		con.full_relax(0.6);
		con.full_relax(0.4);
		con.full_relax(0.2);
	/*	sprintf(q,"output/%04d_p.pov",i);con.draw_particles_pov(q);
		sprintf(q,"gzip -f -9 output/%04d_p.pov",i);system(q);
		sprintf(q,"output/%04d_v.pov",i);con.draw_cells_pov(q);
		sprintf(q,"gzip -f -9 output/%04d_v.pov",i);system(q);*/
	}

	for(k=0;k<400;k++) {
		cout << (double(k)+0.5)*0.02;
		for(i=0;i<=200;i++) cout << " " << u[i][k];
		cout << "\n";
	}
}
