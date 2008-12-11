// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file cmd_line.cc
 * \brief Source code for the command line utility. */

#include "voro++.cc"

// A guess for the memory allocation per region
const int memory=8;

// A maximum allowed number of regions, to prevent enormous amounts of memory
// being allocated
const int max_regions=16777216;

// This message gets displayed if the command line arguments are incorrect
// or if the user requests the help flag
void help_message() {
	cout << "Voro++ version 0.2, by Chris H. Rycroft (UC Berkeley/LBL)\n\n";
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
	cout << " -r         : Assume the input file has an extra coordinate for radii\n";
	cout << " -wc [7]    : Add a cylinder wall object, centered on (x1,x2,x3),\n";
	cout << "              pointing in (x4,x5,x6), radius x7\n";
	cout << " -wo [7]    : Add a conical wall object, apex at (x1,x2,x3), axis\n";
	cout << "              along (x4,x5,x6), angle x7 in radians\n";
	cout << " -ws [4]    : Add a sphere wall object, centered on (x1,x2,x3),\n";
	cout << "              with radius x4\n";
	cout << " -wp [4]    : Add a plane wall object, with normal (x1,x2,x3),\n";
	cout << "              and displacement x4" << endl;
}

// Prints an error message. This is called when the program is unable to make
// sense of the command line options.
void error_message() {
	cerr << "Unrecognized command line options; type \"voro++ -h\" for more information." << endl;
}

// Global variables to set the wall memory allocation, and the current number
// of allocated walls
int wall_mem=init_wall_size,wall_count=0;

// A pointer to the wall pointer array
wall **wp;

// A routine to double up the wall memory allocation if needed
void add_wall_memory() {
	wall **nwp;
	wall_mem*=2;
	if (wall_mem>max_wall_size) cerr << "Too many walls allocated. Try recompiling by boosting the value of max_wall_size in config.hh" << endl;
	nwp=new wall*[wall_mem];
	for(int i=0;i<wall_count;i++) nwp[i]=wp[i];
	delete [] wp;
	wp=nwp;
}

// A routine to deallocate the dynamically created wall objects
void wall_deallocate() {
	for(int i=0;i<wall_count;i++) delete wp[i];
	delete [] wp;
}

int main(int argc,char **argv) {
	int i=1,j=-7;
	bool gnuplot_output=false,neighbor_track=false,polydisperse=false;
	bool xperiodic=false,yperiodic=false,zperiodic=false;
	char buffer[256];
	wp=new wall*[init_wall_size];
	
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
		} else if (strcmp(argv[i],"-ws")==0) {
			if (wall_count==wall_mem) add_wall_memory();
			i++;
			fpoint w0=atof(argv[i++]),w1=atof(argv[i++]);
			fpoint w2=atof(argv[i++]),w3=atof(argv[i]);
			wp[wall_count++]=new wall_sphere(w0,w1,w2,w3,j);
			cout << "Sphere" << w0 << " " << w1 << " " << w2 << " " << w3 << endl;
			j--;
		} else if (strcmp(argv[i],"-wp")==0) {
			if (wall_count==wall_mem) add_wall_memory();
			i++;
			fpoint w0=atof(argv[i++]),w1=atof(argv[i++]);
			fpoint w2=atof(argv[i++]),w3=atof(argv[i]);
			wp[wall_count++]=new wall_plane(w0,w1,w2,w3,j);
			j--;
		} else if (strcmp(argv[i],"-wc")==0) {
			if (wall_count==wall_mem) add_wall_memory();
			i++;
			fpoint w0=atof(argv[i++]),w1=atof(argv[i++]);
			fpoint w2=atof(argv[i++]),w3=atof(argv[i++]);
			fpoint w4=atof(argv[i++]),w5=atof(argv[i++]);
			fpoint w6=atof(argv[i]);
			wp[wall_count++]=new wall_cylinder(w0,w1,w2,w3,w4,w5,w6,j);
			j--;
		} else if (strcmp(argv[i],"-wo")==0) {
			if (wall_count==wall_mem) add_wall_memory();
			i++;
			fpoint w0=atof(argv[i++]),w1=atof(argv[i++]);
			fpoint w2=atof(argv[i++]),w3=atof(argv[i++]);
			fpoint w4=atof(argv[i++]),w5=atof(argv[i++]);
			fpoint w6=atof(argv[i]);
			wp[wall_count++]=new wall_cone(w0,w1,w2,w3,w4,w5,w6,j);
			j--;
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
		} else {
			cerr << "The length scale is smaller than the safe limit of " << tolerance << ".\n";
			cerr << "Either increase the particle length scale, or recompile with a\n";
			cerr << "different limit." << endl;
		}
		wall_deallocate();
		return 0;
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
		wall_deallocate();
		return 0;
	}

	// Prepare output filename
	sprintf(buffer,"%s.vol",argv[i+7]);

	// Now switch depending on whether polydispersity was enabled
	if (polydisperse) {
		container_poly con(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,memory);
		for(j=0;j<wall_count;j++) con.add_wall(*wp[j]);
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
		for(j=0;j<wall_count;j++) con.add_wall(*wp[j]);
		con.import(argv[i+7]);

		if (neighbor_track) con.print_all_neighbor(buffer);
		else con.print_all(buffer);

		if (gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+7]);
			con.draw_cells_gnuplot(buffer);
		}
	}
	return 0;
}
