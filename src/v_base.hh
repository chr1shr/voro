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
		/** The size of a computational block in the x direction. */ 
		const double boxx;
		/** The size of a computational block in the y direction. */ 
		const double boxy;
		/** The size of a computational block in the z direction. */ 
		const double boxz;
		/** The inverse box length in the x direction. */
		const double xsp;
		/** The inverse box length in the y direction. */
		const double ysp;
		/** The inverse box length in the z direction. */
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
