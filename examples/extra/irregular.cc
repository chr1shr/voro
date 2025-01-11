// Irregular packing example code
// By Chris H. Rycroft and the Rycroft Group

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-4,x_max=4;
const double y_min=-4,y_max=4;
const double z_min=-3,z_max=5;

// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

// Set up the number of blocks that the container is divided
// into.
const int n_x=5,n_y=5,n_z=5;

// Create a wall class that, whenever called, will replace the Voronoi cell
// with a prescribed shape, in this case a dodecahedron
class wall_initial_shape : public wall_3d {
    public:
        wall_initial_shape() {

            // Create a dodecahedron
            v.init(-2,2,-2,2,-2,2);
            v.plane(0,Phi,1);v.plane(0,-Phi,1);v.plane(0,Phi,-1);
            v.plane(0,-Phi,-1);v.plane(1,0,Phi);v.plane(-1,0,Phi);
            v.plane(1,0,-Phi);v.plane(-1,0,-Phi);v.plane(Phi,1,0);
            v.plane(-Phi,1,0);v.plane(Phi,-1,0);v.plane(-Phi,-1,0);
        };
        bool point_inside(double x,double y,double z) {return true;}
        bool cut_cell(voronoicell_3d &c,double x,double y,double z) {

            // Set the cell to be equal to the dodecahedron
            c=v;
            return true;
        }
        bool cut_cell(voronoicell_neighbor_3d &c,double x,double y,double z) {

            // Set the cell to be equal to the dodecahedron
            c=v;
            return true;
        }
    private:
        voronoicell_3d v;
};

int main() {

    // Create a container with the geometry given above. This is bigger
    // than the particle packing itself.
    container_3d con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
                     false,false,false,8);

    // Import the irregular particle packing
    con.import("pack_irregular");

    // Save the particles and Voronoi cells in POV-Ray format
    con.draw_particles_pov("irregular_p.pov");
    con.draw_cells_pov("irregular_v1.pov");
    con.print_custom("%i %q %n","irregular1.nei");

    // Create the "initial shape" wall class and add it to the container
    wall_initial_shape wis;
    con.add_wall(wis);
    con.draw_cells_pov("irregular_v2.pov");
    con.print_custom("%i %q %n","irregular2.nei");
}
