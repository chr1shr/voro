// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "cell.cc"
#include "container.cc"

// A guess for the memory allocation per region
const int memory=8;

// A maximum allowed number of regions, to prevent enormous amounts of memory
// being allocated
const int max_regions=16777216;

// This message gets displayed if the command line arguments are incorrect
// or if the user requests the help flag
void help_message() {
	cout << "Syntax: voro++ [opts] <length_scale> <x_min> <x_max> <y_min>\n";
	cout << "                      <y_max> <z_min> <z_max> <filename>\n\n";
	cout << "<length_scale> should be set to a typical particle diameter,\n";
	cout << "or typical distance between particles. It is used to configure\n";
	cout << "the code for maximum efficiency.\n\n";
	cout << "Available options:\n";
	cout << " -g         : Turn on the gnuplot output to <filename.gnu>\n";
	cout << " -h/--help  : Print this information\n";
	cout << " -n         : Turn on the neighbor tracking procedure\n";
	cout << " -p         : Make container periodic in all three directions\n";
	cout << " -px        : Make container periodic in the x direction\n";
	cout << " -py        : Make container periodic in the y direction\n";
	cout << " -pz        : Make container periodic in the z direction\n";
	cout << " -r         : Assume the input file has an extra coordinate for radii" << endl;
}

// 
void error_message() {
	cerr << "Unrecognized command line options; type \"voro++ -h\" for more information." << endl;
}

int main(int argc,char **argv) {
	int i=1;
	bool gnuplot_output=false,neighbor_track=false,polydisperse=false;
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
	if (argc<9) {
	       error_message();return 1;
	}

	// We have enough arguments. Now start searching for command line
	// options.
	while(i+8<argc) {
		if (strcmp(argv[i],"-g")==0) {
			gnuplot_output=true;
		} else if (strcmp(argv[i],"-h")==0||strcmp(argv[i],"--help")==0) {
			help_message();return 0;
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
	fpoint ls=atof(argv[i]);
	fpoint xmin=atof(argv[i+1]),xmax=atof(argv[i+2]);
	fpoint ymin=atof(argv[i+3]),ymax=atof(argv[i+4]);
	fpoint zmin=atof(argv[i+5]),zmax=atof(argv[i+6]);

	// Check that the length scale is positive and reasonably large
	if (ls<tolerance) {
		if (ls<0) {
			cerr << "The length scale must be positive" << endl;
			return 0;
		} else {
			cerr << "The length scale is smaller than the safe limit of " << tolerance << ".\n";
			cerr << "Either increase the particle length scale, or recompile with a\n";
			cerr << "different limit." << endl;
		}
	}
	ls=1.8/ls;

	// Compute the number regions based on the length scale provided. If
	// the total number exceeds a cutoff then bail out, to prevent making
	// a massive memory allocation.
	int nx=int((xmax-xmin)*ls)+1;
	int ny=int((ymax-ymin)*ls)+1;
	int nz=int((zmax-zmin)*ls)+1;
	int nxyz=nx*ny*nz;
	if (nx*ny*nz>max_regions) {
		cerr << "Number of computational blocks (" << nxyz << ") exceeds the maximum\n";
		cerr << "allowed of " << max_regions << ". Either increase the particle\n";
		cerr << "length scale, or recompile with an increased maximum." << endl;
	}

	// Prepare output filename
	sprintf(buffer,"%s.vol",argv[i+7]);

	// Now switch depending on whether polydispersity was enabled
	if (polydisperse) {
		container_poly con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,memory);
		con.import(argv[i+7]);

		if (neighbor_track) con.print_all_neighbor(buffer);
		else con.print_all(buffer);

		if (gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+7]);
			con.draw_cells_gnuplot(buffer);
		}
	} else {
		container con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,memory);
		con.import(argv[i+6]);

		if (neighbor_track) con.print_all_neighbor(buffer);
		else con.print_all(buffer);

		if (gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+6]);
			con.draw_cells_gnuplot(buffer);
		}
	}
	return 0;
}
