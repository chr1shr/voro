// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

#include "voro++.hh"

// Set up the number of blocks that the container is divided into
const int nx=26,ny=26,nz=26;
const double dis=1.25,mjrad=2.5,mirad=0.95;

// Set the number of particles that are going to be randomly introduced
const int particles=100000;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;
	double x,y,z,r;
	voronoicell c;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(-5,5,-5,5,-5,5,nx,ny,nz,true,true,true,8);
	
	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=10*rnd()-5;
		y=10*rnd()-5;
		z=10*rnd()-5;
		r=sqrt((x-dis)*(x-dis)+y*y);
		if((r-mjrad)*(r-mjrad)+z*z<mirad) con.put(i,x,y,z);
		else con.put(i,x,y,z);
	}

	v_loop_subset vls(con);
	vls.setup_intbox(10,15,10,15,10,15);
	
	FILE *fp(voropp_safe_fopen("subpar","w"));
	con.draw_particles(vls,fp);
	fclose(fp);
}
