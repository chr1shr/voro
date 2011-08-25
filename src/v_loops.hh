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

/** A type associated with a v_loop_subset class, determining what type of
 * geometrical region to loop over. */
enum v_loop_subset_mode {
	sphere,
	box,
	no_check
};

/** \brief A class for storing ordering information when particles are added to
 * a container.
 *
 * By default, when particles are added to a container class, they are sorted
 * into an.  */
class voropp_order {
	public:
		/** A pointer to the array holding the ordering. */
		int *o;
		/** A pointer to the next position in the ordering array in
		 * which to store an entry. */
		int *op;
		/** The current memory allocation for the class, set to the
		 * number of entries which can be stored. */
		int size;
		/** The voropp_order constructor allocates memory to store the
		 * ordering information.
		 * \param[in] init_size the initial amount of memory to
		 *                      allocate. */
		voropp_order(int init_size=init_ordering_size)
			: o(new int[init_size<<1]),op(o),size(init_size) {}
		/** The voropp_order destructor frees the dynamically allocated
		 * memory used to store the ordering information. */
		~voropp_order() {
			delete [] o;
		}
		/** Adds a record to the order, corresponding to the memory
		 * address of where a particle was placed into the container.
		 * \param[in] ijk the block into which the particle was placed.
		 * \param[in] q the position within the block where the
		 * 		particle was placed. */
		inline void add(int ijk,int q) {
			if(op==o+size) add_ordering_memory();
			*(op++)=ijk;*(op++)=q;
		}
	private:
		void add_ordering_memory();
};

/** \brief Base class for looping over particles in a container.
 *
 * This class forms the base of all classes that can loop over a subset of
 * particles in a contaner in some order. When initialized, it stores constants
 * about the corresponding container geometry. It also contains a number of
 * routines for interrogating which particle currently being considered by the
 * loop, which are common between all of the derived classes. */
class v_loop_base {
	public:
		/** The number of blocks in the x direction. */
		const int nx;
		/** The number of blocks in the y direction. */
		const int ny;
		/** The number of blocks in the z direction. */
		const int nz;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines that step through blocks in
		 * sequence. */
		const int nxy;
		/** A constant, set to the value of nx*ny*nz, which is used in
		 * the routines that step through blocks in sequence. */
		const int nxyz;
		/** The number of floating point numbers per particle in the
		 * associated container data structure. */
		const int ps;
		/** A pointer to the particle position information in the
		 * associated container data structure. */
		double **p;
		/** A pointer to the particle ID information in the associated
		 * container data structure. */
		int **id;
		/** A pointer to the particle counts in the associated
		 * container data structure. */
		int *co;
		/** The current x-index of the block under consideration by the
		 * loop. */
		int i;
		/** The current y-index of the block under consideration by the
		 * loop. */
		int j;
		/** The current z-index of the block under consideration by the
		 * loop. */
		int k;
		/** The current index of the block under consideration by the
		 * loop. */
		int ijk;
		/** The index of the particle under consideration within the current
		 * block. */
		int q;
		template<class c_class>
		v_loop_base(c_class &con) : nx(con.nx), ny(con.ny), nz(con.nz), nxy(con.nxy), nxyz(con.nxyz), ps(con.ps),
					    p(con.p), id(con.id), co(con.co) {}
		/** Returns the position vector of the particle currently being
		 * considered by the loop.
		 * \param[out] (x,y,z) the position vector. */
		inline void pos(double &x,double &y,double &z) {
			double *pp(p[ijk]+ps*q);
			x=*(pp++);y=*(pp++);z=*pp;
		}
		/** Returns the ID, position vector, and radius of the particle
		 * currently being considered by the loop.
		 * \param[out] pid the particle ID.
		 * \param[out] (x,y,z) */
		inline void pos(int &pid,double &x,double &y,double &z,double &r) {
			pid=id[ijk][q];
			double *pp(p[ijk]+ps*q);
			x=*(pp++);y=*(pp++);z=*pp;
			r=ps==3?default_radius:*(++pp);
		}
		/** Returns the x position of the particle currently being
		 * considered by the loop. */
		inline double x() {return p[ijk][ps*q];}
		/** Returns the y position of the particle currently being
		 * considered by the loop. */
		inline double y() {return p[ijk][ps*q+1];}
		/** Returns the z position of the particle currently being
		 * considered by the loop. */
		inline double z() {return p[ijk][ps*q+2];}
		/** Returns the ID of the particle currently being considered
		 * by the loop. */
		inline double pid() {return id[ijk][q];}
};

/** \brief Class for looping over all of the particles in a container.
 * 
 * This is one of the simplest loop classes, that scans the computational
 * blocks in order, and scans all the particles within each block in order. */
class v_loop_all : public v_loop_base {
	public:
		template<class c_class>
		v_loop_all(c_class &con) : v_loop_base(con) {}
		/** Sets the class to consider the first particle.
		 * \return True if there is any particle to consider, false
		 * otherwise. */
		inline bool start() {
			i=j=k=ijk=q=0;
			while(co[ijk]==0) if(!next_block()) return false;
			return true;
		}
		/** Finds the next particle to test.
		 * \return True if there is another particle, false if no more
		 * particles are available. */		
		inline bool inc() {
			q++;
			if(q>=co[ijk]) {
				q=0;
				do {
					if(!next_block()) return false;
				} while(co[ijk]==0);
			}
			return true;
		}
	private:
		inline bool next_block() {
			ijk++;
			i++;
			if(i==nx) {
				i=0;j++;
				if(j==ny) {
					j=0;k++;
					if(ijk==nxyz) return false;
				}
			}
			return true;
		}
};

/** \brief Class for looping over a subset of particles in a container.
 *
 * This class can loop over a subset of particles in a certain geometrical
 * region within the container. The class can be set up to */ 
class v_loop_subset : public v_loop_base {
	public:
		/** The current mode of operation, determining whether tests
		 * should be applied to particles to ensure they are within a
		 * certain geometrical object. */
		v_loop_subset_mode mode;
		template<class c_class>
		v_loop_subset(c_class &con) : v_loop_base(con), ax(con.ax), ay(con.ay), az(con.az),
			sx(con.bx-ax), sy(con.by-ay), sz(con.bz-az), xsp(con.xsp), ysp(con.ysp), zsp(con.zsp),
			xperiodic(con.xperiodic), yperiodic(con.yperiodic), zperiodic(con.zperiodic) {}
		void setup_sphere(double vx,double vy,double vz,double r,bool bounds_test=true);
		void setup_box(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,bool bounds_test=true);
		void setup_intbox(int ai_,int bi_,int aj_,int bj_,int ak_,int bk_);
		bool start();
		/** Finds the next particle to test.
		 * \return True if there is another particle, false if no more
		 * particles are available. */
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
		void setup_common();
		bool next_block();
		bool out_of_bounds();
};

/** \brief Class for looping over all of the particles specified in a
 * pre-assembled voropp_order class.
 *
 * The voropp_order class can be used to create a specific order of particles
 * within the container. This class can then loop over these particles in this
 * order. The class is particularly useful in cases where the ordering of the
 * output must match the ordering of particles as they were inserted into the
 * container. */
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
		/** Finds the next particle to test.
		 * \return True if there is another particle, false if no more
		 * particles are available. */
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

/** \brief A class for looping over all particles in a container_periodic or
 * container_periodic_poly class.
 *
 * Since the container_periodic and container_periodic_poly classes have a
 * fundamentally different memory organization, the regular loop classes cannot
 * be used with them. */
class v_loop_all_periodic : public v_loop_base {
	public:
		template<class c_class>
		v_loop_all_periodic(c_class &con) : v_loop_base(con), ey(con.ey), ez(con.ez), wy(con.wy), wz(con.wz),
			ijk0(nx*(ey+con.oy*ez)), inc2(2*nx*con.ey+1) {}
		inline bool start() {
			i=0;
			j=ey;
			k=ez;
			ijk=ijk0;
			q=0;
			while(co[ijk]==0) if(!next_block()) return false;
			return true;
		}
		/** Finds the next particle to test.
		 * \return True if there is another particle, false if no more
		 * particles are available. */		
		inline bool inc() {
			q++;
			if(q>=co[ijk]) {
				q=0;
				do {
					if(!next_block()) return false;
				} while(co[ijk]==0);
			}
			return true;
		}
	private:
		/** The lower y index (inclusive) of the primary domain within
		 * the block structure. */
		int ey;
		/** The lower y index (inclusive) of the primary domain within
		 * the block structure. */
		int ez;
		/** The upper y index (exclusive) of the primary domain within
		 * the block structure. */
		int wy;
		/** The upper z index (exclusive) of the primary domain within
		 * the block structure. */
		int wz;
		/** The index of the (0,0,0) block within the block structure.
		 */
		int ijk0;
		/** A value to increase ijk by when the z index is increased.
		 */
		int inc2;		
		inline bool next_block() {
			i++;
			if(i==nx) {
				i=0;j++;
				if(j==wy) {
					j=ey;k++;
					if(k==wz) return false;
					ijk+=inc2;
				} else ijk++;
			} else ijk++;
			return true;
		}
};

#endif
