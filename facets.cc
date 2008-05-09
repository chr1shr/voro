// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : January 21st 2008

#include "cell.cc"
#include "container.cc"
#include "container_poly.cc"

// A guess for the memory allocation; should be bigger than 27/(diameter^3).
const int memory=32;

// This message gets displayed if the command line arguments are incorrect
// or if the user requests the help flag
void help_message() {
	cout << "Syntax: facets [opts] <x_min> <x_max> <y_min> <y_max> <z_min> <z_max> <filename>" << endl;
	cout << endl;
	cout << "Available options:" << endl;
	cout << " -g         : Turn on the gnuplot output to <filename.gnu>" << endl;
	cout << " -h/--help  : Print this information" << endl;
	cout << " -l <value> : Pick a typical particle diameter length scale" << endl;
	cout << " -n         : Turn on the neighbor tracking procedure" << endl;
	cout << " -p         : Make container periodic in all three directions" << endl;
	cout << " -px        : Make container periodic in the x direction" << endl;
	cout << " -py        : Make container periodic in the y direction" << endl;
	cout << " -pz        : Make container periodic in the z direction" << endl;
	cout << " -r         : Assume the input file has an extra coordinate for radii" << endl;
}

// 
void error_message() {
	cerr << "Unrecognized command line options; type \"facets -h\" for more information." << endl;
}

int main(int argc,char **argv) {
	double ls=1,ils=1/3.0;
	int i=1;
	bool spec_length=false,gnuplot_output=false;
	bool neighbor_track=false,polydisperse=false;
	bool xperiodic=false,yperiodic=false,zperiodic=false;
	char buffer[256];
	
	// If there's one argument, check to see if it's requesting help.
	// Otherwise, bail out with an error.
	if (argc==2) {
		if (strcmp(argv[1],"-h")==0||strcmp(argv[1],"--help")==0) {
			help_message();return 0;
		} else {
			error_message();return 1;
		}
	}
	
	// If there aren't enough command line arguments, then bail out
	// with an error.
	if (argc<8) {
	       error_message();return 1;
	}

	// We have enough arguments. Now start searching for command line
	// options.
	while(i+7<argc) {
		if (strcmp(argv[i],"-g")==0) {
			gnuplot_output=true;
		} else if (strcmp(argv[i],"-h")==0||strcmp(argv[i],"--help")==0) {
			help_message();return 0;
		} else if (strcmp(argv[i],"-l")==0) {
			i++;
			spec_length=true;
			ls=atof(argv[i]);
			ils=1/(3.0*ls);
		} else if (strcmp(argv[i],"-n")==0) {
			neighbor_track=true;
		} else if (strcmp(argv[i],"-p")==0) {
			xperiodic=yperiodic=zperiodic=true;
		} else if (strcmp(argv[i],"-px")==0) {
			xperiodic=true;
		} else if (strcmp(argv[i],"-py")==0) {
			yperiodic=true;
		} else if (strcmp(argv[i],"-pz")==0) {
			zperiodic=true;
		} else if (strcmp(argv[i],"-r")==0) {
			polydisperse=true;
		} else {
			error_message();return 1;
		}
		i++;
	}

	// Read in the dimensions of the test box, and estimate the number of
	// boxes to divide the region up into
	fpoint xmin=atof(argv[i]),xmax=atof(argv[i+1]);
	fpoint ymin=atof(argv[i+2]),ymax=atof(argv[i+3]);
	fpoint zmin=atof(argv[i+4]),zmax=atof(argv[i+5]);

	int nx=int((xmax-xmin)*ils)+1;
	int ny=int((ymax-ymin)*ils)+1;
	int nz=int((zmax-zmin)*ils)+1;

	// Prepare output filename
	sprintf(buffer,"%s.vol",argv[i+6]);

	// Now switch depending on whether polydispersity was enabled
	if (polydisperse) {
		container_poly con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,memory);
		con.import(argv[i+6]);

		if (spec_length) con.length_scale=ls;
		else con.guess_length_scale();

		if (neighbor_track) con.print_all_neighbor(buffer);
		else con.print_all(buffer);

		if (gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+6]);
			con.draw_gnuplot(buffer);
		}
	} else {
		container con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,memory);
		con.import(argv[i+6]);

		if (spec_length) con.length_scale=ls;
		else con.guess_length_scale();

		if (neighbor_track) con.print_all_neighbor(buffer);
		else con.print_all(buffer);

		if (gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+6]);
			con.draw_gnuplot(buffer);
		}
	}
	return 0;
}
