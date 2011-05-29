// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file pre_container.hh
 * \brief Header file for the pre_container and related classes. */

#ifndef VOROPP_PRE_CONTAINER_HH
#define VOROPP_PRE_CONTAINER_HH

#include "container.hh"

class pre_container_base {
	public:
		/** The minimum x coordinate of the container. */
		const double ax;
		/** The maximum x coordinate of the container. */
		const double bx;
		/** The minimum y coordinate of the container. */
		const double ay;
		/** The maximum y coordinate of the container. */
		const double by;
		/** The minimum z coordinate of the container. */
		const double az;
		/** The maximum z coordinate of the container. */
		const double bz;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** A boolean value that determines if the z coordinate in
		 * periodic or not. */
		const bool zperiodic;
		void guess_optimal(int &nx,int &ny,int &nz);
		pre_container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int ps_);
		~pre_container_base();
		inline int total_particles() {
			return (end_id-pre_id)*pre_container_chunk_size+(ch_id-*end_id);
		}
	protected:
		const int ps;
		void new_chunk();
		void extend_chunk_index();
		int index_sz;
		int **pre_id,**end_id,**l_id,*ch_id,*e_id;
		double **pre_p,**end_p,*ch_p;
};

class pre_container : public pre_container_base {
	public:
		pre_container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				bool xperiodic_,bool yperiodic_,bool zperiodic_)
			: pre_container_base(ax_,bx_,ay_,by_,az_,bz_,xperiodic_,yperiodic_,zperiodic_,3) {};
		void put(int n,double x,double y,double z);
		void import(FILE *fp=stdin);
		inline void import(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(fp);
			fclose(fp);
		}
		void setup(container &con);
};

class pre_container_poly : public pre_container_base {
	public:
		pre_container_poly(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				bool xperiodic_,bool yperiodic_,bool zperiodic_)
			: pre_container_base(ax_,bx_,ay_,by_,az_,bz_,xperiodic_,yperiodic_,zperiodic_,4) {};
		void put(int n,double x,double y,double z,double r);
		void import(FILE *fp=stdin);
		inline void import(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(fp);
			fclose(fp);
		}
		void setup(container_poly &con);
};

#endif
