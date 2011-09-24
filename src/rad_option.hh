// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file rad_option.hh
 * \brief Header file for the classes encapsulating functionality for the
 * regular and radical Voronoi tessellations. */

#ifndef VOROPP_RAD_OPTION_HH
#define VOROPP_RAD_OPTION_HH

namespace voro {

/** This class contains all of the routines that are specific to computing
 * the regular Voronoi tessellation. */
class radius_mono {
	protected:
		inline void r_init(int ijk,int s) {}
		inline void r_prime(double rv) {}
		inline bool r_ctest(double crs,double mrs) {return crs>mrs;}
		inline double r_cutoff(double lrs) {return lrs;}
		inline double r_max_add(double rs) {return rs;}
		inline double r_current_sub(double rs,int ijk,int q) {return rs;}
		inline double r_scale(double rs,int ijk,int q) {return rs;}
		inline bool r_scale_check(double &rs,double mrs,int ijk,int q) {return rs<mrs;}
};

/** This class contains all of the routines that are specific to computing
 * the radical Voronoi tessellation. It contains routines that scale the plane
 * positions based on the particle radii, as well as routines that are used
 * to cut off the computation. */
class radius_poly {
	public:
		/** A two-dimensional array holding particle positions and radii. */			
		double **ppr;
		/** The current maximum radius of any particle, used to
		 * determine when to cut off the radical Voronoi computation.
		 * */
		double max_radius;
		radius_poly() : max_radius(0) {}
	protected:
		double r_rad,r_mul,r_val;
		inline void r_init(int ijk,int s) {
			r_rad=ppr[ijk][4*s]*ppr[ijk][4*s];
			r_mul=r_rad-max_radius*max_radius;
		}
		inline void r_prime(double rv) {r_val=1+r_mul/rv;}
		inline bool r_ctest(double crs,double mrs) {return crs+r_mul>sqrt(mrs*crs);}
		inline double r_cutoff(double x) {return x*r_val;}
		inline double r_max_add(double rs) {return rs+max_radius*max_radius;}
		inline double r_current_sub(double rs,int ijk,int q) {
			return rs-ppr[ijk][4*q]*ppr[ijk][4*q];
		}
		inline double r_scale(double rs,int ijk,int q) {
			return rs+r_rad-ppr[ijk][4*q]*ppr[ijk][4*q];
		}
		inline bool r_scale_check(double &rs,double mrs,int ijk,int q) {
			double trs=rs;
			rs+=r_rad-ppr[ijk][4*q]*ppr[ijk][4*q];
			return rs<sqrt(mrs*trs);
		}
};

}
#endif
