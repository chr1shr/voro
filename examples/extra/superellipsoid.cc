// Superellipsoid example code
// By Chris H. Rycroft and the Rycroft Group

#include "voro++.hh"
using namespace voro;

// This function returns a random floating point number uniformly distributed
// over a range
inline double rnd(double a,double b) {
    return a+(b-a)/RAND_MAX*static_cast<double>(rand());
}

int main() {
    double x,y,z,rsq,r;
    voronoicell_3d v;

    // Initialize the Voronoi cell to be a cube of side length 2, centered
    // on the origin
    v.init(-1,1,-1,1,-1,1);

    // Cut the cell by 5000 random planes that are scaled to create a
    // superellipsoid
    for(int i=0;i<5000;i++) {
        x=rnd(-1,1);
        y=rnd(-1,1);
        z=rnd(-1,1);
        rsq=x*x*x*x+y*y*y*y+z*z*z*z;
        if(rsq>1e-4&&rsq<1) {
            r=1/sqrt(sqrt(rsq));
            x*=r;y*=r;z*=r;
            v.plane(x*x*x,y*y*y,z*z*z,x*x*x*x+y*y*y*y+z*z*z*z);
        }
    }

    // Output the Voronoi cell to a file, in the gnuplot format
    v.draw_gnuplot(0,0,0,"superellipsoid.gnu");

    // Output the Voronoi cell to a file in POV-Ray formats
    v.draw_pov(0,0,0,"superellipsoid_v.pov");
    v.draw_pov_mesh(0,0,0,"superellipsoid_m.pov");
}
