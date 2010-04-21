// Single Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "cell_2d.cc"
#include "container_2d.cc"

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;double x,y;
	container_2d con(0,1,0,1,4,4,false,false,8);
x=rnd();
	// Cut the cell by 1000 random planes which are all a distance 1 away
	// from the origin, to make an approximation to a sphere
	for(i=0;i<1000;i++) {
		x=rnd();
		y=rnd();
//		x=0.5+(i*0.0008)*sin(i*0.01);
//		y=0.5+(i*0.0008)*cos(i*0.01);
		con.put(i,x,y);
	}

	// Output the particle positions to a file
	con.draw_particles("test_ctr.par");

	// Output the Voronoi cells to a file, in the gnuplot format
	con.draw_cells_gnuplot("test_ctr.gnu");
}
