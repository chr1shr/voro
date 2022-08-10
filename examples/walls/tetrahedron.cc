// Tetrahedron example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-2,x_max=2;
const double y_min=-2,y_max=2;
const double z_min=-2,z_max=2;

// Set up the number of blocks that the container is divided
// into
const int n_x=7,n_y=7,n_z=7;

// Set the number of particles that are going to be randomly
// introduced
const int particles=64;

// This function returns a random floating point number uniformly distributed
// over a range
inline double rnd(double a,double b) {
    return a+(b-a)/RAND_MAX*static_cast<double>(rand());
}

int main() {
    int i=0;
    double x,y,z;

    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for 8
    // particles within each computational block.
    container_3d con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                     false,false,false,8);

    // Add four plane walls to the container to make a tetrahedron
    wall_plane p1(1,1,1,1);con.add_wall(p1);
    wall_plane p2(-1,-1,1,1);con.add_wall(p2);
    wall_plane p3(1,-1,-1,1);con.add_wall(p3);
    wall_plane p4(-1,1,-1,1);con.add_wall(p4);

    // Randomly insert particles into the container, checking that they lie
    // inside the tetrahedron
    while(i<particles) {
        x=rnd(x_min,x_max);
        y=rnd(y_min,y_max);
        z=rnd(z_min,z_max);
        if(con.point_inside(x,y,z)) con.put(i++,x,y,z);
    }

    // Output the particle positions and the Voronoi cells in Gnuplot and
    // POV-Ray formats
    con.draw_particles("tetrahedron_p.gnu");
    con.draw_cells_gnuplot("tetrahedron_v.gnu");
    con.draw_particles_pov("tetrahedron_p.pov");
    con.draw_cells_pov("tetrahedron_v.pov");
}
