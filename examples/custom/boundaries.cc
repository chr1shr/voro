// Radical Voronoi tessellation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

#include "voro++.cc"

// Set up constants for the container geometry
const double x_min=-6,x_max=6;
const double y_min=-6,y_max=6;
const double z_min=-3,z_max=9;

// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

// Set up the number of blocks that the container is divided
// into.
const int n_x=5,n_y=5,n_z=5;

// The voronoicell to use as the base shape
voronoicell v;

// Modification of the container class
class container_mod : public container_base {
	public:
		container_mod(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem) : container_base(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,xperiodic_,yperiodic_,zperiodic_,init_mem,3),
		vc(*this,xperiodic_?2*nx_+1:nx_,yperiodic_?2*ny_+1:ny_,zperiodic_?2*nz_+1:nz_) {}
		template<class v_cell>
		inline bool initialize_voronoicell(v_cell &c,int ijk,int q,int ci,int cj,int ck,
			int &i,int &j,int &k,double &x,double &y,double &z,int &disp) {
			c=v;
			double *pp(p[ijk]+ps*q);
			x=*(pp++);y=*(pp++);z=*pp;
			i=xperiodic?nx:ci;
			j=xperiodic?ny:cj;
			k=xperiodic?nz:ck;
			if(!apply_walls(c,x,y,z)) return false;
			disp=ijk-i-nx*(j+ny*k);
			return true;
		}
		void put(int n,double x,double y,double z) {
			int ijk;
			if(put_locate_block(ijk,x,y,z)) {
				id[ijk][co[ijk]]=n;
				double *pp(p[ijk]+3*co[ijk]++);
				*(pp++)=x;*(pp++)=y;*pp=z;
			}
		}
		void import(FILE *fp) {
			int i,j;
			double x,y,z;
			while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(i,x,y,z);
			if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
		}
		inline void import(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(fp);
			fclose(fp);
		}		
		template<class v_cell,class v_loop>
		inline bool compute_cell(v_cell &c,v_loop &vl) {
			return vc.compute_cell(c,vl.ijk,vl.q,vl.i,vl.j,vl.k);
		}
	private:
		voropp_compute<container_mod> vc;
		inline void r_init(int ijk,int s) {};
		inline double r_cutoff(double lrs) {return lrs;}
		inline double r_max_add(double rs) {return rs;}
		inline double r_current_sub(double rs,int ijk,int q) {return rs;}
		inline double r_scale(double rs,int ijk,int q) {return rs;}
		friend class voropp_compute<container_mod>;		
};

int main() {

	// Create a container with the geometry given above. This is bigger
	// than the particle packing itself.
	container_mod con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	con.import("pack_six_cube");

	// Create an dodecahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(0,Phi,1);
	v.plane(0,-Phi,1);
	v.plane(0,Phi,-1);
	v.plane(0,-Phi,-1);
	v.plane(1,0,Phi);
	v.plane(-1,0,Phi);
	v.plane(1,0,-Phi);
	v.plane(-1,0,-Phi);
	v.plane(Phi,1,0);
	v.plane(-Phi,1,0);
	v.plane(Phi,-1,0);
	v.plane(-Phi,-1,0);

	double x,y,z;
	voronoicell c;
	FILE *f1(voropp_safe_fopen("boundaries_v.pov","w"));
	FILE *f2(voropp_safe_fopen("boundaries_p.pov","w"));
	v_loop_all vl(con);
	if(vl.start()) do if(con.compute_cell(c,vl)) {
		vl.pos(x,y,z);
		c.draw_pov(x,y,z,f1);
		fprintf(f2,"sphere{<%g,%g,%g>,s}\n",x,y,z);
	} while (vl.inc());
	fclose(f1);
	fclose(f2);
}
