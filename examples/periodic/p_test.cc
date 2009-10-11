// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// Set up the number of blocks that the container is divided into
const int n_x=4,n_y=4,n_z=4;
const fpoint dx1=0.3,dy1=0.3;
const fpoint dx2=-2,dy2=0;

int main() {

	int i,j,q=0;
	fpoint x,y;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(0,10,0,10,-1,1,n_x,n_y,n_z,
			false,false,false,8);

	for(i=-20;i<20;i++) {
		for(j=-20;j<20;j++) {
			x=5.00001+i*dx1+j*dx2;
			y=5.00000459+i*dy1+j*dy2;
			if(x>0&&x<10&&y>0&&y<10) {
				cout << x << " " << y << endl;
				con.put(q,x,y,0);
				q++;
			}
		}
	}

	con.print_network("ptest.network");
	con.draw_particles("ptest.par");
	con.draw_cells_gnuplot("ptest.out");
}
