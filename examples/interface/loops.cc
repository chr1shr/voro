// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

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
	container con(-5,5,-5,5,-5,5,nx,ny,nz,false,false,false,8);
	voropp_order vo;
	
	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=10*rnd()-5;
		y=10*rnd()-5;
		z=10*rnd()-5;
		r=sqrt((x-dis)*(x-dis)+y*y);
		if((r-mjrad)*(r-mjrad)+z*z<mirad) con.put(vo,i,x,y,z);
		else con.put(i,x,y,z);
	}

	con.draw_particles("draw");

	FILE *f1(voropp_safe_fopen("loop1_p.pov","w"));
	FILE *f2(voropp_safe_fopen("loop1_v.pov","w"));
	v_loop_order vlo(con,vo);
	if(vlo.start()) do if(con.compute_cell(c,vlo)) {
		vlo.pos(x,y,z);
		
		//fprintf(f1,"sphere{<%g,%g,%g>,r}\n",x,y,z);
		c.draw_pov_mesh(x,y,z,f1);
		c.draw_pov(x,y,z,f2);
	} while (vlo.inc());
	fclose(f1);
	fclose(f2);

	f1=voropp_safe_fopen("loop2_p.pov","w");
	f2=voropp_safe_fopen("loop2_v.pov","w");
	v_loop_all vla(con);
	if(vla.start()) do {
		vla.pos(x,y,z);
		r=sqrt((x+dis)*(x+dis)+z*z);
		if((r-mjrad)*(r-mjrad)+y*y<mirad&&con.compute_cell(c,vla)) {
		//	fprintf(f1,"sphere{<%g,%g,%g>,r}\n",x,y,z);
			c.draw_pov_mesh(x,y,z,f1);
			c.draw_pov(x,y,z,f2);
		}
	} while (vla.inc());
	fclose(f1);
	fclose(f2);
}
