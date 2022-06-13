#include "voro++_2d.cc"
using namespace voro;

#include <vector>

int main() {

	// Initialize the container class to be the unit square, with
	// non-periodic boundary conditions. Divide it into a 6 by 6 grid, with
	// an initial memory allocation of 16 particles per grid square.
	container_2d con(-5,5,-5,5,6,6,false,false,16);

	// Import the spiral data set
	con.import("2dbox.pro");

    int id;
    double x,y,r;
    voronoicell_neighbor_2d c;
    std::vector<int> q;
    std::vector<double> p;
    c_loop_all_2d vl(con);
    if(vl.start()) do {
        if(con.compute_cell(c,vl)) {
            vl.pos(id,x,y,r);
            c.vertices(p);
            c.neighbors(q);
            int n=p.size()/2,j;
            for(int i=0;i<n;i++) {
                j=c.ed[2*i];
                if((q[i]<36&&id>=36)) {
                    //if(p[2*i]<p[2*j]) {
                     printf("sphere{<%g,%g,zz>,r}\ncylinder{<%g,%g,zz>,<%g,%g,zz>,r}\n",x+p[2*i],y+p[2*i+1],x+p[2*i],y+p[2*i+1],x+p[2*j],y+p[2*j+1]);
                    //}
                }
            }
        }
    } while(vl.inc());

	// Do a custom computation on the Voronoi cells, printing the IDs,
	// positions, and Voronoi cell areas to a file
	con.draw_cells_gnuplot("2dbox.vor");
}
