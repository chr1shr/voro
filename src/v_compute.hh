// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file v_compute.hh
 * \brief Header file for the container_base template and related classes. */

#ifndef VOROPP_V_COMPUTE_HH
#define VOROPP_V_COMPUTE_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "config.hh"
#include "cell.hh"

template <class c_class>
class voropp_compute {
	public:
		/** A reference to the container class on which to carry out*/
		c_class &con;
		/** The size of an internal computational block in the x
		 * direction. */
		const double boxx;
		/** The size of an internal computational block in the y
		 * direction. */
		const double boxy;
		/** The size of an internal computational block in the z
		 * direction. */
		const double boxz;
		/** The inverse box length in the x direction, set to
		 * nx/(bx-ax). */
		const double xsp;
		/** The inverse box length in the y direction, set to
		 * ny/(by-ay). */
		const double ysp;
		/** The inverse box length in the z direction, set to
		 * nz/(bz-az). */
		const double zsp;
		/** The number of boxes in the x direction for the searching mask. */
		const int hx;
		/** The number of boxes in the y direction for the searching mask. */
		const int hy;
		/** The number of boxes in the z direction for the searching mask. */
		const int hz;
		/** A constant, set to the value of hx multiplied by hy, which
		 * is used in the routines which step through mask boxes in
		 * sequence. */
		const int hxy;
		/** A constant, set to the value of hx*hy*hz, which is used in
		 * the routines which step through mask boxes in sequence. */
		const int hxyz;
		/** The number of floating point entries to store for each
		 * particle. */
		const int ps;
		const int hgrid,fgrid,hgridsq,seq_length;
		/** This array holds the numerical IDs of each particle in each
		 * computational box. */
		int **id;
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		double **p;
		/** An array holding the number of particles within each
		 * computational box of the container. */
		int *co;		
		voropp_compute(c_class &con_,int hx_,int hy_,int hz_);
		/** The class destructor frees the dynamically allocated memory
		 * for the mask and queue. */
		~voropp_compute() {
			delete [] qu;
			delete [] mask;
		}
		template<class v_cell>
		bool compute_cell(v_cell &c,int ijk,int s,int i,int j,int k,double x,double y,double z);
	private:
		/** A constant set to boxx*boxx+boxy*boxy+boxz*boxz, which is
		 * frequently used in the computation. */
		const double bxsq;
		/** This sets the current value being used to mark tested blocks
		 * in the mask. */
		unsigned int mv;		
		/** The current size of the search list. */
		int qu_size;
		const unsigned int *wl;
		double *mrad;
		/** This array is used during the cell computation to determine
		 * which blocks have been considered. */
		unsigned int *mask;		
		/** An array is used to store the queue of blocks to test
		 * during the Voronoi cell computation. */
		int *qu;
		/** A pointer to the end of the queue array, used to determine
		 * when the queue is full. */
		int *qu_l;
		template<class v_cell>
		bool corner_test(v_cell &c,double xl,double yl,double zl,double xh,double yh,double zh);
		template<class v_cell>
		inline bool edge_x_test(v_cell &c,double x0,double yl,double zl,double x1,double yh,double zh);
		template<class v_cell>
		inline bool edge_y_test(v_cell &c,double xl,double y0,double zl,double xh,double y1,double zh);
		template<class v_cell>
		inline bool edge_z_test(v_cell &c,double xl,double yl,double z0,double xh,double yh,double z1);
		template<class v_cell>
		inline bool face_x_test(v_cell &c,double xl,double y0,double z0,double y1,double z1);
		template<class v_cell>
		inline bool face_y_test(v_cell &c,double x0,double yl,double z0,double x1,double z1);
		template<class v_cell>
		inline bool face_z_test(v_cell &c,double x0,double y0,double zl,double x1,double y1);
		bool compute_min_max_radius(int di,int dj,int dk,double fx,double fy,double fz,double gx,double gy,double gz,double& crs,double mrs);
		void add_list_memory(int*& qu_s,int*& qu_e);
		inline void reset_mask() {
			for(unsigned int *mp(mask);mp<mask+hxyz;mp++) *mp=0;
		}		
};

#endif
