// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file v_loops.hh
 * \brief Header file for the container_base template and related classes. */

#ifndef VOROPP_V_LOOPS_HH
#define VOROPP_V_LOOPS_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "config.hh"

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

/*class v_loop_subset {
	public:
		template<class c_class>
		v_loop_all(c_class &con) : v_loop_base(con) {}
		inline void start();
		inline bool end();
		inline void operator++();
};*/

class v_loop_order : public v_loop_base {
	public:
		voropp_order &vo;
		int *cp,*op;
		template<class c_class>
		v_loop_order(c_class &con,voropp_order &vo_) : v_loop_base(con), vo(vo_), nx(con.nx), nxy(con.nxy) {}
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
