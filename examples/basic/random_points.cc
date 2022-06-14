// Voronoi calculation example code
// Example code for Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(z_max-z_min);

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=20;

// This function returns a random number uniformly distributed over the range
// from a to b
inline double rnd(double a,double b) {
    return a+(b-a)/RAND_MAX*static_cast<double>(rand());
}

int main() {
    int i;
    double x,y,z;

    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block
    container_3d con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                     false,false,false,8);

    // Randomly add particles into the container
    for(i=0;i<particles;i++) {
        x=rnd(x_min,x_max);
        y=rnd(y_min,y_max);
        z=rnd(z_min,z_max);
        con.put(i,x,y,z);
    }

    // Sum up the volumes, and check that this matches the container volume
    double vvol=con.sum_cell_volumes();
    printf("Container volume : %g\n"
           "Voronoi volume   : %g\n"
           "Difference       : %g\n",cvol,vvol,vvol-cvol);

    // Output the particle positions in gnuplot format
    con.draw_particles("random_points_p.gnu");

    // Output the Voronoi cells in gnuplot format
    con.draw_cells_gnuplot("random_points_v.gnu");
}
