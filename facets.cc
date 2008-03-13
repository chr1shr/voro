// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : January 21st 2007

#include "cell.cc"
#include "container.cc"

// A guess for the memory allocation; should be bigger than 125/(diameter^3).
const int memory=256;

// This message gets displayed if the command line arguments are incorrect
void helpmessage() {
		cout << "Syntax: facets [-p] <x_min> <x_max> <y_min> <y_max> <z_min> <z_max> <filename>" << endl;
}

int main(int argc,char **argv) {
	int i;bool periodic;char buffer[256];
	
	// Check for the periodicity flag, and make sure there are the correct
	// number of arguments
	if(argc==8) {
		i=1;periodic=false;
	} else if (argc==9) {
		if (strcmp(argv[1],"-p")==0) {periodic=true;i=2;}
		else {helpmessage();return 1;}
	} else {helpmessage();return 1;}

	// Read in the dimensions of the test box, and estimate the number of
	// boxes to divide the region up into
	fpoint xmin=atof(argv[i]),xmax=atof(argv[i+1]);
	fpoint ymin=atof(argv[i+2]),ymax=atof(argv[i+3]);
	fpoint zmin=atof(argv[i+4]),zmax=atof(argv[i+5]);

	int nx=int((xmax-xmin)*0.2)+1;
	int ny=int((ymax-ymin)*0.2)+1;
	int nz=int((zmax-zmin)*0.2)+1;

	// Create a container according to the specifications above
	container con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
			periodic,periodic,periodic,memory);

	// Import the particles
	con.import(argv[i+6]);

	// Print out a list of the particles, and their Voronoi volumes
	sprintf(buffer,"%s.vol",argv[i+6]);
	con.vprintall(buffer);

	// Save the Voronoi network of all the particles to a text file in a
	// format ready for plotting by gnuplot
	sprintf(buffer,"%s.gnu",argv[i+6]);
	con.vdraw_gnuplot(buffer);
	return 0;
}
