// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file v_loops.hh
 * \brief Header file for the loop classes. */

#ifndef VOROPP_V_LOOPS_HH
#define VOROPP_V_LOOPS_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "config.hh"

enum v_loop_subset_mode {
	sphere,
	box,
	no_check
};

class voropp_order {
	public:
		int size;
		int *o,*op;
		voropp_order(int init_size=init_ordering_size)
			: size(init_size),o(new int[size<<1]),op(o) {}
		~voropp_order() {
			delete [] o;
		}
		inline void add(int ijk,int q) {
			if(op==o+size) add_ordering_memory();
			*(op++)=ijk;*(op++)=q;
		}
	private:
		void add_ordering_memory();
};

class v_loop_base {
	public:
		const int nx,ny,nz,nxy,nxyz,ps;
		double **p;
		int **id;
		int *co;
		int i,j,k,ijk,q;
		template<class c_class>
		v_loop_base(c_class &con) : nx(con.nx), ny(con.ny), nz(con.nz), nxy(con.nxy), nxyz(con.nxyz), ps(con.ps),
					    p(con.p), id(con.id), co(con.co) {}
		inline void pos(double &x,double &y,double &z) {
			double *pp(p[ijk]+ps*q);
			x=*(pp++);y=*(pp++);z=*pp;
		}
		inline double x() {return p[ijk][ps*q];}
		inline double y() {return p[ijk][ps*q+1];}
		inline double z() {return p[ijk][ps*q+2];}
		inline double pid() {return id[ijk][q];}
};

class v_loop_all : public v_loop_base {
	public:
		template<class c_class>
		v_loop_all(c_class &con) : v_loop_base(con) {}
		inline bool start() {
			i=j=k=ijk=q=0;
			while(co[ijk]==0) {
				ijk++;
				i++;
				if(i==nx) {
					i=0;j++;
					if(j==ny) {
						j=0;k++;
						if(ijk==nxyz) return false;
					}
				}
			}
			return true;
		}
		inline bool inc() {
			q++;
			while(q>=co[ijk]) {
				ijk++;
				i++;
				if(i==nx) {
					i=0;j++;
					if(j==ny) {
						j=0;k++;
						if(ijk==nxyz) return false;
					}
				}
				q=0;
			}
			return true;
		}
};

class v_loop_subset : public v_loop_base {
	public:
		v_loop_subset_mode mode;
		template<class c_class>
		v_loop_subset(c_class &con) : v_loop_base(con), ax(con.ax), ay(con.ay), az(con.az),
			sx(con.bx-ax), sy(con.by-ay), sz(con.bz-az), xsp(con.xsp), ysp(con.ysp), zsp(con.zsp),
			xperiodic(con.xperiodic), yperiodic(con.yperiodic), zperiodic(con.zperiodic) {}
		bool start_sphere(double vx,double vy,double vz,double r,bool bounds_test=true);
		bool start_box(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,bool bounds_test=true);
		bool start_intbox(int ai_,int bi_,int aj_,int bj_,int ak_,int bk_);
		/** Finds the next point to test.
		 * \return True if there is another point, false if no more points are
		 * available. */
		inline bool inc() {
			do {
				q++;
				while(q>=co[ijk]) {q=0;if(!next_block()) return false;}
			} while(mode!=no_check&&out_of_bounds());
			return true;
		}
	private:
		const double ax,ay,az,sx,sy,sz,xsp,ysp,zsp;
		const bool xperiodic,yperiodic,zperiodic;
		double px,py,pz,apx,apy,apz;
		double v0,v1,v2,v3,v4,v5;
		int ai,bi,aj,bj,ak,bk,s;
		int ci,cj,ck,di,dj,dk,inc1,inc2;
		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
		inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
		inline int step_int(double a) {return a<0?int(a)-1:int(a);}
		bool start_common();
		bool next_block();
		bool out_of_bounds();
};

class v_loop_order : public v_loop_base {
	public:
		voropp_order &vo;
		int *cp,*op;
		template<class c_class>
		v_loop_order(c_class &con,voropp_order &vo_)
		: v_loop_base(con), vo(vo_), nx(con.nx), nxy(con.nxy) {}
		inline bool start() {
			cp=vo.o;op=vo.op;
			if(cp!=op) {
				ijk=*(cp++);decode();
				q=*(cp++);
				return true;
			} else return false;
		}
		inline bool inc() {
			if(cp==op) return false;
			ijk=*(cp++);decode();
			q=*(cp++);
			return true;
		}
	private:
		const int nx;
		const int nxy;
		inline void decode() {
			k=ijk/nxy;
			int ijkt=ijk-nxy*k;
			j=ijkt/nx;
			i=ijkt-j*nx;
		}
};

#endif
