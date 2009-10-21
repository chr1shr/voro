// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// A guess for the memory allocation per region
const int memory=8;

// A maximum allowed number of regions, to prevent enormous amounts of memory
// being allocated
const int max_regions=16777216;

// A buffer size
const int bsize=256;


int main(int argc,char **argv) {

	char buffer[bsize];
	int i,n;
	fpoint bx,bxy,by,bxz,byz,bz,x,y,z;

	if(argc!=2) {
		cerr << "Syntax: ./network <filename>" << endl;
		return VOROPP_CMD_LINE_ERROR;
	}

	ifstream is;
	is.open(argv[1],ifstream::in);
	if(is.fail()) voropp_fatal_error("Unable to open file for import",VOROPP_FILE_ERROR);

	// Read header line
	is.getline(buffer,bsize);
	cout << buffer << endl;
	if(strcmp(buffer,"Unit cell vectors:")!=0)
		voropp_fatal_error("Invalid header line",VOROPP_FILE_ERROR);
	is.width(bsize);

	is >> buffer;if(strcmp(buffer,"va=")!=0) voropp_fatal_error("Invalid first vector",VOROPP_FILE_ERROR);
	is >> bx >> x >> x;

	is >> buffer;if(strcmp(buffer,"vb=")!=0) voropp_fatal_error("Invalid second vector",VOROPP_FILE_ERROR);
	is >> bxy >> by >> x;
	
	is >> buffer;if(strcmp(buffer,"vc=")!=0) voropp_fatal_error("Invalid third vector",VOROPP_FILE_ERROR);
	is >> bxz >> byz >> bz;
	
	is >> n;

	printf("%f %f %f %f %f %f\n",bx,bxy,by,bxz,byz,bz);

	if(n<1) voropp_fatal_error("Invalid number of particles",VOROPP_FILE_ERROR);

	if(bx<tolerance||by<tolerance||bz<tolerance)
		voropp_fatal_error("Invalid box dimensions",VOROPP_FILE_ERROR);

	fpoint ls=1.8*pow(bx*by*bz,-1.0/3.0);

	fpoint nxf=bx*ls+1.5;
	fpoint nyf=by*ls+1.5;
	fpoint nzf=bz*ls+1.5;
	if (nxf*nyf*nzf>max_regions) {
		cerr << "voro++: Number of computational blocks exceeds the maximum allowed of " << max_regions << ".\n";
		cerr << "Either increase the particle length scale, or recompile with an increased\nmaximum." << endl;
		return VOROPP_MEMORY_ERROR;
	}

	// Now that we are confident that the number of regions is reasonable,
	// create integer versions of them
	int nx=int(nxf);
	int ny=int(nyf);
	int nz=int(nzf);
	printf("%f %d %d %d %d\n",ls,n,nx,ny,nz);

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_periodic con(bx,bxy,by,bxz,byz,bz,nx,ny,nz,memory);

	for(i=0;i<n;i++) {
		is >> buffer;
		is >> x >> y >> z;
		con.put(i,x,y,z);
	}

	cout << "Volume check:" << endl;
	cout << "  Total domain volume  = " << bx*by*bz << endl;
	cout << "  Total Voronoi volume = " << con.sum_cell_volumes() << endl;

	// Save the Voronoi network of all the particles to text files
	// in gnuplot and POV-Ray formats
	
	con.draw_particles("YUG.par");
	con.draw_cells_gnuplot("YUG.out");
	con.print_network("YUG.network");
}
