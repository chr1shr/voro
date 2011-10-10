// Voro++, a cell-based Voronoi library
//
// Authors  : Chris H. Rycroft (LBL / UC Berkeley)
//            Cody Robert Dance (UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++_2d.cc"
using namespace voro;

#include <vector>
using namespace std;

const double pad=0.01;

int main(int argc,char **argv) {

	if(argc!=2) {
		fprintf(stderr,"Syntax: voro++_bdry <input_file>\n");
		return 1;
	}

	// Read file to determine the size of the system and the number
	// of particles
	unsigned int i=0,j;
	int id;
	double x,y,minx=large_number,maxx=-minx,miny=minx,maxy=maxx;
	vector<int> vid;
	vector<double> vpos;
	vector<char> vbd;
	bool boundary_track=false,start=false;
	char *buf(new char[512]);

	FILE *fp=safe_fopen(argv[1],"r");

	while(fgets(buf,512,fp)!=NULL) {

		if(strcmp(buf,"#Start\n")==0||strcmp(buf,"# Start\n")==0) {

			// Check that two consecutive start tokens haven't been
			// encountered
			if(boundary_track) voro_fatal_error("File import error - two consecutive start tokens found",VOROPP_FILE_ERROR);
			start=true;boundary_track=true;

		} else if(strcmp(buf,"#End\n")==0||strcmp(buf,"# End\n")==0||
			  strcmp(buf,"#End")==0||strcmp(buf,"# End")==0) {
			
			// Check that two consecutive end tokens haven't been
			// encountered
			if(start) voro_fatal_error("File import error - end token immediately after start token",VOROPP_FILE_ERROR);
			if(!boundary_track) voro_fatal_error("File import error - found end token without start token",VOROPP_FILE_ERROR);
			vbd[i-1]|=2;boundary_track=false;
		} else {

			// Try and read three entries from the line
			if(sscanf(buf,"%d %lg %lg",&id,&x,&y)!=3) voro_fatal_error("File import error #1",VOROPP_FILE_ERROR);
			vid.push_back(id);
			vpos.push_back(x);
			vpos.push_back(y);
			vbd.push_back(start?1:0);
			i++;

			// Determine bounds
			if(x<minx) minx=x;
			if(x>maxx) maxx=x;
			if(y<miny) miny=y;
			if(y>maxy) maxy=y;

			start=false;
		}
	}
	
	if(boundary_track) voro_fatal_error("File import error - boundary not finished",VOROPP_FILE_ERROR);
	if(!feof(fp)) voro_fatal_error("File import error #2",VOROPP_FILE_ERROR);
	delete [] buf;	
	
	// Add small amount of padding to container bounds
	double dx=maxx-minx,dy=maxy-miny;
	minx-=pad*dx;maxx+=pad*dx;dx+=2*pad*dx;
	miny-=pad*dy;maxy+=pad*dy;dy+=2*pad*dy;
	
	// Guess the optimal computationl grid, aiming at eight particles per
	// grid square
	double lscale=8*sqrt(dx*dy)/i;
	int nx=int(dx/lscale)+1,ny=int(dy/lscale)+1;
	
	// Print diagnostic information
	printf("Container bounds : %g<x<%g, %g<%g\n"
	       "Total particles  : %d\n"
	       "Compute grid     : %d by %d\n",minx,maxx,miny,maxy,i,nx,ny);

	// Create container
	container_boundary_2d con(minx,maxx,miny,maxy,nx,ny,false,false,16);

	// Import data
	for(j=0;j<vid.size();j++) {
		if(vbd[j]&1) con.start_boundary();
		con.put(vid[j],vpos[2*j],vpos[2*j+1]);
		if(vbd[j]&2) con.end_boundary();
	}

	// Carry out all of the setup prior to computing any Voronoi cells
	con.setup();

	// Compute the Voronoi cells and save them to file
	char *buf2(new char[strlen(argv[1])+5]);
	sprintf(buf2,"%s.gnu",argv[1]);
	con.draw_cells_gnuplot(buf2);
	delete [] buf2;

	// Output the neighbor mesh in gnuplot format
/*	FILE *ff=safe_fopen("sphere_mesh.net","w");
	vector<int> vi;
	voronoicell_nonconvex_neighbor_2d c;
	c_loop_all_2d cl(con);
	if(cl.start()) do if(con.compute_cell(c,cl)) {
		i=cl.pid();
		c.neighbors(vi);
		for(j=0;j<vi.size();l++) if(vi[l]>i)
			fprintf(ff,"%g %g\n%g %g\n\n\n",
				p[2*i],p[2*i+1],p[2*vi[l]],p[2*vi[l]+1]
	} while (cl.inc());
	fclose(ff);*/

	return 0;
}
