// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double boxl=1.2;

// Set up the number of blocks that the container is divided into
const int bl=24;

// Set the number of particles that are going to be randomly introduced
const int particles=4000;

const int nface=11;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

struct wall_shell : public wall {
	public:
		wall_shell(double xc_,double yc_,double zc_,double rc,double sc,int w_id_=-99)
			: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), lc(rc-sc), uc(rc+sc) {}
		bool point_inside(double x,double y,double z) {
			double rsq=(x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc);
			return rsq>lc*lc&&rsq<uc*uc;
		}
		template<class v_cell>
		bool cut_cell_base(v_cell &c,double x,double y,double z) {
			double xd=x-xc,yd=y-yc,zd=z-zc,dq=xd*xd+yd*yd+zd*zd,dq2;
			if (dq>1e-5) {
				dq2=2*(sqrt(dq)*lc-dq);
				dq=2*(sqrt(dq)*uc-dq);
				return c.nplane(xd,yd,zd,dq,w_id)&&c.nplane(-xd,-yd,-zd,-dq2,w_id);
			}
			return true;
		}
		bool cut_cell(voronoicell &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		const double xc,yc,zc,lc,uc;
};


int main() {
	int i=0,l;
	double x,y,z,r,dx,dy,dz;
	int faces[nface],*fp;
	double p[3*particles];

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(-boxl,boxl,-boxl,boxl,-boxl,boxl,bl,bl,bl,false,false,false,8);

	wall_shell ws(0,0,0,1,0.00001);
	con.add_wall(ws);

	// Randomly add particles into the container
	while(i<particles) {
		x=boxl*(2*rnd()-1);
		y=boxl*(2*rnd()-1);
		z=boxl*(2*rnd()-1);
		r=x*x+y*y+z*z;
		if(r>1e-5) {
			r=1/sqrt(r);x*=r;y*=r;z*=r;
			con.put(i,x,y,z);
			i++;
		}
	}

	for(l=4;l<10000;l++) {
		c_loop_all vl(con);
		voronoicell c;
		for(fp=faces;fp<faces+nface;fp++) *fp=0;
		if(vl.start()) do if(con.compute_cell(c,vl)) {
			vl.pos(i,x,y,z,r);
			c.centroid(dx,dy,dz);
			p[3*i]=x+dx;
			p[3*i+1]=y+dy;
			p[3*i+2]=z+dz;

			i=c.number_of_faces()-4;
			if(i<0) i=0;if(i>=nface) i=nface-1;
			faces[i]++;
		} while (vl.inc());
		con.clear();
		double fac=l<9000?0.1/sqrt(double(l)):0;
		for(i=0;i<particles;i++) con.put(i,p[3*i]+fac*(2*rnd()-1),p[3*i+1]+fac*(2*rnd()-1),p[3*i+2]+fac*(2*rnd()-1));
		printf("%d",l);
		for(fp=faces;fp<faces+nface;fp++) printf(" %d",*fp);
		puts("");
	}

	// Output the particle positions in gnuplot format
	con.draw_particles("sphere_mesh_p.gnu");

	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("sphere_mesh_v.gnu");

	FILE *ff=safe_fopen("sphere_mesh.net","w");
	vector<int> vi;
	voronoicell_neighbor c;
	c_loop_all vl(con);
	if(vl.start()) do if(con.compute_cell(c,vl)) {
		i=vl.pid();
		c.neighbors(vi);
		for(l=0;l<(signed int) vi.size();l++) if(vi[l]>i)
			fprintf(ff,"%g %g %g\n%g %g %g\n\n\n",
				p[3*i],p[3*i+1],p[3*i+2],
				p[3*vi[l]],p[3*vi[l]+1],p[3*vi[l]+2]);
	} while (vl.inc());
	fclose(ff);
}
