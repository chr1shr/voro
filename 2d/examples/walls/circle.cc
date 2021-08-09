#define _USE_MATH_DEFINES
#include <cmath>

#include "voro++_2d.hh"
using namespace voro;

const double radius=0.7;

/** Returns a uniformly distributed random number over a given interval.
 * \param[in] (a,b) the interval bounds. */
inline double rnd(double a,double b) {return a+(b-a)/RAND_MAX*static_cast<double>(rand());}

int main() {
	int i=0;double x,y;

	// Initialize the container class to be the unit square, with
	// non-periodic boundary conditions. Divide it into a 6 by 6 grid, with
	// an initial memory allocation of 16 particles per grid square.
	container_2d con(-1,1,-1,1,10,10,false,false,16);

	// Add circular wall object
	wall_circle wc(0,0,radius);
	con.add_wall(wc);

	// Add 1000 random points to the container
	while(i<1000) {
		x=rnd(-1,1);
		y=rnd(-1,1);
		if(con.point_inside(x,y)) con.put(i++,x,y);
	}

	// Output the particle positions to a file
	con.draw_particles("circle.par");

	// Output the Voronoi cells to a file, in the gnuplot format
	con.draw_cells_gnuplot("circle.gnu");

	con.print_custom("%i %q %a %n","circle.vol");

	// Sum the Voronoi cell areas and compare to the circle area
	double carea=M_PI*radius*radius,varea=con.sum_cell_areas();
	printf("Total circle area       : %g\n"
	       "Total Voronoi cell area : %g\n"
	       "Difference              : %g\n",carea,varea,varea-carea);
}
