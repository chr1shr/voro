// Voro++, a cell-based Voronoi library
//
// Authors  : Chris H. Rycroft (LBL / UC Berkeley)
//            Cody Robert Dance (UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++_2d.hh"
using namespace voro;

#include <vector>
#include <cstring>
using namespace std;

const double pad=0.01;

int main(int argc,char **argv) {

	if(argc<2) {
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
	int mid=0;
	bool neg_label=false,boundary_track=false,start=false;
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
			if(id<0) neg_label=true;
			if(id>mid) mid=id;
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
	double lscale=sqrt(8.0*dx*dy/i);
	int nx=int(dx/lscale)+1,ny=int(dy/lscale)+1;
	
	// Print diagnostic information
	printf("Container bounds : [%g:%g] [%g:%g]\n"
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

	// Save the boundary in a format that can be read by Gnuplot
	char *outfn(new char[strlen(argv[1])+5]);
	sprintf(outfn,"%s.bd",argv[1]);

	con.draw_boundary_gnuplot(outfn);

	//Test Global Connectivity Info
//	cout << "beginning" << endl;
//	char *connectfn(new char[strlen(argv[1])+5]);
//	sprintf(connectfn,"%s.connect",argv[1]);
//	con.print_custom(argv[2],connectfn );
//	cout << "finished" << endl;

	// Compute the Voronoi cells and save them to file
	sprintf(outfn,"%s.gnu",argv[1]);
	con.draw_cells_gnuplot(outfn);

	// Output the neighbor mesh in gnuplot format
	if(mid>16777216) puts("Network output disabled due to too high Ids");
	else if(neg_label) puts("Network output disabled due to negative IDs");
	else {

		int *mp=new int[mid+1];
		for(j=0;j<i;j++) mp[vid[j]]=j;

		sprintf(outfn,"%s.net",argv[1]);
		FILE *ff=safe_fopen(outfn,"w");
		int l1,l2;
		vector<int> vi;
		voronoicell_nonconvex_neighbor_2d c;
		c_loop_all_2d cl(con);
		if(cl.start()) do if(con.compute_cell(c,cl)) {
			id=cl.pid();l1=2*mp[id];
			c.neighbors(vi);
			for(j=0;j<vi.size();j++) if(vi[j]>id) {
				l2=2*mp[vi[j]];
				fprintf(ff,"%g %g\n%g %g\n\n\n",
					vpos[l1],vpos[l1+1],vpos[l2],vpos[l2+1]);
			}
		} while (cl.inc());
		fclose(ff);

		delete [] mp;
	}
//	delete [] connectfn;
	delete [] outfn;
	return 0;
}
