// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file v_base.hh
 * \brief Header file for the base Voronoi container class. */

#ifndef VOROPP_V_BASE_HH
#define VOROPP_V_BASE_HH

class voropp_base {
	public:
		const int nx;
		const int ny;
		const int nz;
		const int nxy;
		const int nxyz;
		const double boxx;
		const double boxy;
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
		/** An array to hold the minimum distances associated with the
		 * worklists. This array is initialized during container
		 * construction, by the initialize_radii() routine. */
		double *mrad;
		bool contains_neighbor(const char* format);
		voropp_base(int nx_,int ny_,int nz_,double boxx_,double boxy_,double boxz_);
		~voropp_base() {delete [] mrad;}
#include "worklist.hh"
	protected:
		inline int step_int(double a) {return a<0?int(a)-1:int(a);}
		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
		inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
	private:
		void compute_minimum(double &minr,double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi,int ti,int tj,int tk);
};

#endif
