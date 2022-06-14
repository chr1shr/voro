// Single Voronoi cell example code
// Example code for Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#include "voro++.hh"
using namespace voro;

// This function returns a random number uniformly distributed over the range
// from a to b
inline double rnd(double a,double b) {
    return a+(b-a)/RAND_MAX*static_cast<double>(rand());
}

int main() {
    double x,y,z,rsq,r;
    voronoicell_3d v;

    // Initialize the Voronoi cell to be a cube of side length 2, centered
    // on the origin
    v.init(-1,1,-1,1,-1,1);

    // Cut the cell by 250 random planes which are all a distance 1 away
    // from the origin, to make an approximation to a sphere
    for(int i=0;i<250;i++) {
        x=rnd(-1,1);
        y=rnd(-1,1);
        z=rnd(-1,1);
        rsq=x*x+y*y+z*z;
        if(rsq>0.01&&rsq<1) {
            r=1/sqrt(rsq);x*=r;y*=r;z*=r;
            v.plane(x,y,z,1);
        }
    }

    // Output the Voronoi cell to a file, in the gnuplot format
    v.draw_gnuplot(0,0,0,"single_cell.gnu");
}
