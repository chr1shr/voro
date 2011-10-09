#include "voro++_2d.cc"
using namespace voro;

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {

	// Initialize the container class to be the unit square, with
	// non-periodic boundary conditions. Divide it into a 6 by 6 grid, with
	// an initial memory allocation of 16 particles per grid square.
	container_boundary_2d con(-1,1,-1,1,12,12,false,false,8);
	
	// Add 1000 random points to the container
	con.start_boundary();
	con.put(0,-0.4,-0.6);
	con.put(1,0.8,-0.8);
	con.put(2,0.4,0.6);
	con.put(3,-0.6,0.6);
	con.end_boundary();
	
	con.put(4,0.2,0.1);
	con.put(5,-0.2,0);

	con.draw_boundary_gnuplot("container_bd.gnu");
	con.draw_particles("container_bd.par");

	con.setup();
	con.draw_cells_gnuplot("container_bd_v.gnu");

	// Sum the Voronoi cell areas and compare to the container area
//	double carea=1,varea=con.sum_cell_areas();
//	printf("Total container area    : %g\n"
//	       "Total Voronoi cell area : %g\n"
//	       "Difference              : %g\n",carea,varea,varea-carea);
}
