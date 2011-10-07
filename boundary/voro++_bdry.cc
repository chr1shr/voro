// Voro++, a cell-based Voronoi library
//
// Authors  : Chris H. Rycroft (LBL / UC Berkeley)
//            Cody Robert Dance (UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++_2d.hh"
using namespace voro;

int main(int argc,char **argv) {

	if(argc!=2) {
		fprintf(stderr,"Syntax: voro++_bdry <input_file>\n");
		return 1;
	}

	container_boundary_2d con(0,1,600,1100,8,8,false,false,8);
	con.import(argv[1]);
	con.setup();

	char *buf(new char[strlen(argv[1])+5]);
	sprintf(buf,"%s.gnu",argv[1]);
	con.draw_cells_gnuplot(buf);
	delete [] buf;

	return 0;
}
