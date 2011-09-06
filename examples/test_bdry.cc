#include "voro++_2d.hh"
#include <iostream>

int main(int argc,char **argv) {

	if(argc!=2) {
		fprintf(stderr,"Syntax: test_bdry <input_file>\n");
		return 1;
	}

	container_2d con(0,500,600,1100,8,8,false,false,false,16);
	con.import(argv[1]);
	con.setup();
	con.debug_output();

	char *buf(new char[strlen(argv[1])+5]);
	sprintf(buf,"%s.gnu",argv[1]);
	con.draw_cells_gnuplot(buf);

	delete [] buf;
	return 0;
}
