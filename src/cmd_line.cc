// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file cmd_line.cc
 * \brief Source code for the command-line utility. */

#include <cstring>
#include "voro++.hh"

enum blocks_mode {
	none,
	length_scale,
	specified
};

// A maximum allowed number of regions, to prevent enormous amounts of memory
// being allocated
const int max_regions=16777216;

// This message gets displayed if the user requests the help flag
void help_message() {
	puts("Voro++ version 0.4, by Chris H. Rycroft (UC Berkeley/LBL)\n\n"
	     "Syntax: voro++ [options] <x_min> <x_max> <y_min>\n"
	     "               <y_max> <z_min> <z_max> <filename>\n\n"
	     "Available options:\n"
	     " -c <str>   : Specify a custom output string\n"
	     " -g         : Turn on the gnuplot output to <filename.gnu>\n"
	     " -h/--help  : Print this information\n"
	     " -hc        : Print information about custom output\n"
	     " -l <len>   : Manually specify a length scale to configure the internal\n"
	     "              computational grid\n"
	     " -m <mem>   : Manually choose the memory allocation per grid block\n"
	     "              (default 8)\n"
	     " -n [3]     : Manually specify the internal grid size\n"
//	     " -o         : Ensure that the output file has the same order as the input\n"
//	     "              file\n"
	     " -p         : Make container periodic in all three directions\n"
	     " -px        : Make container periodic in the x direction\n"
	     " -py        : Make container periodic in the y direction\n"
	     " -pz        : Make container periodic in the z direction\n"
	     " -r         : Assume the input file has an extra coordinate for radii\n"
	     " -wc [7]    : Add a cylinder wall object, centered on (x1,x2,x3),\n"
	     "              pointing in (x4,x5,x6), radius x7\n"
	     " -wo [7]    : Add a conical wall object, apex at (x1,x2,x3), axis\n"
	     "              along (x4,x5,x6), angle x7 in radians\n"
	     " -ws [4]    : Add a sphere wall object, centered on (x1,x2,x3),\n"
	     "              with radius x4\n"
	     " -wp [4]    : Add a plane wall object, with normal (x1,x2,x3),\n"
	     "              and displacement x4");
}

// This message gets displayed if the user requests information about doing custom
// output
void custom_output_message() {
	puts("The \"-c\" option allows a string to be specified that will customize the output\n"
	     "file to contain a variety of statistics about each computed Voronoi cell. The\n"
	     "string is similar to the standard C printf() function, made up of text with\n"
	     "additional control sequences that begin with percentage signs that are expanded\n"
	     "to different statistics. See http://math.lbl.gov/voro++/doc/custom.html for more\n"
	     "information.\n"
	     "\nParticle-related:\n"
	     "  %i The particle ID number\n"
	     "  %x The x coordinate of the particle\n"
	     "  %y The y coordinate of the particle\n"
	     "  %z The z coordinate of the particle\n"
	     "  %q The position vector of the particle, short for \"%x %y %z\"\n"
	     "  %r The radius of the particle (only printed if -p enabled)\n"
	     "\nVertex-related:\n"
	     "  %w The number of vertices in the Voronoi cell\n"
	     "  %p A list of the vertices of the Voronoi cell in the format (x,y,z),\n"
	     "     relative to the particle center\n"
	     "  %P A list of the vertices of the Voronoi cell in the format (x,y,z),\n"
	     "     relative to the global coordinate system\n"
	     "  %o A list of the orders of each vertex\n"
	     "  %m The maximum radius squared of a vertex position, relative to the\n"
	     "     particle center\n"
	     "\nEdge-related:\n"
	     "  %g The number of edges of the Voronoi cell\n"
	     "  %E The total edge distance\n"
	     "  %e A list of perimeters of each face\n"
	     "\nFace-related:\n"
	     "  %s The number of faces of the Voronoi cell\n"
	     "  %F The total surface area of the Voronoi cell\n"
	     "  %A A frequency table of the number of edges for each face\n"
	     "  %a A list of the number of edges for each face\n"
	     "  %f A list of areas of each face\n"
	     "  %t A list of bracketed sequences of vertices that make up each face\n"
	     "  %l A list of normal vectors for each face\n"
	     "  %n A list of neighboring particle or wall IDs corresponding to each face\n"
	     "\nVolume-related:\n"
	     "  %v The volume of the Voronoi cell\n"
	     "  %c The centroid of the Voronoi cell, relative to the particle center\n"
	     "  %C The centroid of the Voronoi cell, in the global coordinate system");
}

// Prints an error message. This is called when the program is unable to make
// sense of the command-line options.
void error_message() {
	fputs("voro++: Unrecognized command-line options; type \"voro++ -h\" for more\ninformation.\n",stderr);
}

int main(int argc,char **argv) {
	int i=1,j=-7,custom_output=0,nx,ny,nz,init_mem(8);
	double ls=0;
	blocks_mode bm(none);
	bool gnuplot_output=false,polydisperse=false;
	bool xperiodic=false,yperiodic=false,zperiodic=false;
	pre_container *pcon(NULL);pre_container_poly *pconp(NULL);
	char buffer[256];
	wall_list wl;

	// If there's one argument, check to see if it's requesting help.
	// Otherwise, bail out with an error.
	if(argc==2) {
		if(strcmp(argv[1],"-h")==0||strcmp(argv[1],"--help")==0) {
			help_message();return 0;
		} else if(strcmp(argv[1],"-hc")==0) {
			custom_output_message();return 0;
		} else {
			error_message();
			return VOROPP_CMD_LINE_ERROR;
		}
	}

	// If there aren't enough command-line arguments, then bail out
	// with an error.
	if(argc<7) {
	       error_message();
	       return VOROPP_CMD_LINE_ERROR;
	}

	// We have enough arguments. Now start searching for command-line
	// options.
	while(i<argc-7) {
		if(strcmp(argv[i],"-c")==0) {
			if(i>=argc-9) {error_message();wl.deallocate();return VOROPP_CMD_LINE_ERROR;}
			if(custom_output==0) {
				custom_output=++i;
			} else {
				fputs("voro++: multiple custom output strings detected\n",stderr);
				wl.deallocate();
				return VOROPP_CMD_LINE_ERROR;
			}
		} else if(strcmp(argv[i],"-g")==0) {
			gnuplot_output=true;
		} else if(strcmp(argv[i],"-h")==0||strcmp(argv[i],"--help")==0) {
			help_message();return 0;
		} else if(strcmp(argv[i],"-hc")==0) {
			custom_output_message();return 0;
		} else if(strcmp(argv[i],"-l")==0) {
			if(bm!=none) {
				fputs("voro++: Conflicting options about grid setup (-l/-n)\n",stderr);
				wl.deallocate();
				return VOROPP_CMD_LINE_ERROR;
			}
			bm=length_scale;
			i++;ls=atof(argv[i]);
		} else if(strcmp(argv[i],"-m")==0) {
			i++;init_mem=atoi(argv[i]);
		} else if(strcmp(argv[i],"-n")==0) {
			if(bm!=none) {
				fputs("voro++: Conflicting options about grid setup (-l/-n)\n",stderr);
				wl.deallocate();
				return VOROPP_CMD_LINE_ERROR;
			}
			bm=specified;
			i++;
			nx=atoi(argv[i++]);			
			ny=atoi(argv[i++]);			
			nz=atoi(argv[i]);
			if(nx<=0||ny<=0||nz<=0) {
				fputs("voro++: Computational grid specified with -n must be greater than one\n"
				      "in each direction\n",stderr);
				wl.deallocate();
				return VOROPP_CMD_LINE_ERROR;
			}		
		} else if(strcmp(argv[i],"-p")==0) {
			xperiodic=yperiodic=zperiodic=true;
		} else if(strcmp(argv[i],"-px")==0) {
			xperiodic=true;
		} else if(strcmp(argv[i],"-py")==0) {
			yperiodic=true;
		} else if(strcmp(argv[i],"-pz")==0) {
			zperiodic=true;
		} else if(strcmp(argv[i],"-r")==0) {
			polydisperse=true;
		} else if(strcmp(argv[i],"-ws")==0) {
			if(i>=argc-12) {error_message();wl.deallocate();return VOROPP_CMD_LINE_ERROR;}
			i++;
			double w0=atof(argv[i++]),w1=atof(argv[i++]);
			double w2=atof(argv[i++]),w3=atof(argv[i]);
			wl.add_wall(new wall_sphere(w0,w1,w2,w3,j));
			j--;
		} else if(strcmp(argv[i],"-wp")==0) {
			if(i>=argc-12) {error_message();wl.deallocate();return VOROPP_CMD_LINE_ERROR;}
			i++;
			double w0=atof(argv[i++]),w1=atof(argv[i++]);
			double w2=atof(argv[i++]),w3=atof(argv[i]);
			wl.add_wall(new wall_plane(w0,w1,w2,w3,j));
			j--;
		} else if(strcmp(argv[i],"-wc")==0) {
			if(i>=argc-15) {error_message();wl.deallocate();return VOROPP_CMD_LINE_ERROR;}
			i++;
			double w0=atof(argv[i++]),w1=atof(argv[i++]);
			double w2=atof(argv[i++]),w3=atof(argv[i++]);
			double w4=atof(argv[i++]),w5=atof(argv[i++]);
			double w6=atof(argv[i]);
			wl.add_wall(new wall_cylinder(w0,w1,w2,w3,w4,w5,w6,j));
			j--;
		} else if(strcmp(argv[i],"-wo")==0) {
			if(i>=argc-15) {error_message();wl.deallocate();return VOROPP_CMD_LINE_ERROR;}
			i++;
			double w0=atof(argv[i++]),w1=atof(argv[i++]);
			double w2=atof(argv[i++]),w3=atof(argv[i++]);
			double w4=atof(argv[i++]),w5=atof(argv[i++]);
			double w6=atof(argv[i]);
			wl.add_wall(new wall_cone(w0,w1,w2,w3,w4,w5,w6,j));
			j--;
		} else {
			error_message();
			return VOROPP_CMD_LINE_ERROR;
		}
		i++;
	}

	// Check the memory guess is positive
	if(init_mem<=0) {
		fputs("voro++: The memory allocation must be positive\n",stderr);
		wl.deallocate();
		return VOROPP_CMD_LINE_ERROR;
	}

	// Read in the dimensions of the test box, and estimate the number of
	// boxes to divide the region up into
	double ax=atof(argv[i]),bx=atof(argv[i+1]);
	double ay=atof(argv[i+2]),by=atof(argv[i+3]);
	double az=atof(argv[i+4]),bz=atof(argv[i+5]);
	printf("%g %g %g %g %g %g\n",ax,bx,ay,by,az,bz);

	// Check that for each coordinate, the minimum value is smaller
	// than the maximum value
	if(bx<ax) {
		fputs("voro++: Minimum x coordinate exceeds maximum x coordinate\n",stderr);
		wl.deallocate();
		return VOROPP_CMD_LINE_ERROR;
	}
	if(by<ay) {
		fputs("voro++: Minimum y coordinate exceeds maximum y coordinate\n",stderr);
		wl.deallocate();
		return VOROPP_CMD_LINE_ERROR;
	}
	if(bz<az) {
		fputs("voro++: Minimum z coordinate exceeds maximum z coordinate\n",stderr);
		wl.deallocate();
		return VOROPP_CMD_LINE_ERROR;
	}

	if(bm==none) {
		if(polydisperse) {
			pconp=new pre_container_poly(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
			pconp->import(argv[i+6]);
			pconp->guess_optimal(nx,ny,nz);
		} else {
			pcon=new pre_container(ax,bx,ay,by,az,bz,xperiodic,yperiodic,zperiodic);
			pcon->import(argv[i+6]);
			pcon->guess_optimal(nx,ny,nz);
		}
	} else {
		double nxf,nyf,nzf;
		if(bm==length_scale) {

			// Check that the length scale is positive and
			// reasonably large
			if(ls<tolerance) {
				fputs("voro++: ",stderr);
				if(ls<0) {
					fputs("The length scale must be positive\n",stderr);
				} else {
					fprintf(stderr,"The length scale is smaller than the safe limit of %g. Either\nincrease the particle length scale, or recompile with a different limit.\n",tolerance);
				}
				wl.deallocate();
				return VOROPP_CMD_LINE_ERROR;
			}
			ls=0.6/ls;
			nxf=(bx-ax)*ls+1;
			nyf=(by-ay)*ls+1;
			nzf=(bz-az)*ls+1;

			nx=int(nxf);ny=int(nyf);nz=int(nzf);
		} else {
			nxf=nx;nyf=ny;nzf=nz;
		}

		// Compute the number regions based on the length scale
		// provided. If the total number exceeds a cutoff then bail
		// out, to prevent making a massive memory allocation. Do this
		// test using floating point numbers, since huge integers could
		// potentially wrap around to negative values.
		if(nxf*nyf*nzf>max_regions) {
			fprintf(stderr,"voro++: Number of computational blocks exceeds the maximum allowed of %d.\n"
				       "Either increase the particle length scale, or recompile with an increased\nmaximum.",max_regions);
			wl.deallocate();
			return VOROPP_MEMORY_ERROR;
		}
	}
	printf("%g %d %d %d\n",ls,nx,ny,nz);

	// Prepare output filename
	sprintf(buffer,"%s.vol",argv[i+6]);

	// Now switch depending on whether polydispersity was enabled
	if(polydisperse) {
		container_poly con(ax,bx,ay,by,az,bz,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,init_mem);
		con.add_wall(wl);
		if(bm==none) {
			pconp->setup(con);delete pconp;
		} else con.import(argv[i+6]);

		if(custom_output>0) con.print_custom(argv[custom_output],buffer);
		else con.print_custom("%i %q %v",buffer);

		if(gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+6]);
			con.draw_cells_gnuplot(buffer);
		}
	} else {
		container con(ax,bx,ay,by,az,bz,nx,ny,nz,
			      xperiodic,yperiodic,zperiodic,init_mem);
		con.add_wall(wl);
		if(bm==none) {
			pcon->setup(con);delete pcon;
		} else con.import(argv[i+6]);

		if(custom_output>0) con.print_custom(argv[custom_output],buffer);
		else con.print_custom("%i %q %v",buffer);

		if(gnuplot_output) {
			sprintf(buffer,"%s.gnu",argv[i+6]);
			con.draw_cells_gnuplot(buffer);
		}
	}
	return 0;
}
