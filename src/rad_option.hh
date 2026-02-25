// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file rad_option.hh
 * \brief Header file for the classes encapsulating functionality for the
 * regular and radical Voronoi tessellations. */

#ifndef VOROPP_RAD_OPTION_HH
#define VOROPP_RAD_OPTION_HH

#include <cstdio>

#ifdef _OPENMP
#include "omp.h"

inline int max_t(int thread_number) {
    if(thread_number!=0) {
        return thread_number;
    }
    else {
        return omp_get_max_threads();
    }
}
inline int t_num() {return omp_get_thread_num();}
inline double wtime() {return omp_get_wtime();}
#else
inline int max_t(int thread_number) {
    return 1;
}
inline int t_num() {
    return 0;
}
inline double wtime() {
    return 0.;
}
#endif

#include <cmath>

namespace voro {

/** \brief Class containing all of the routines that are specific to computing
 * the regular Voronoi tessellation.
 *
 * The container and container_periodic classes are derived from this class,
 * and during the Voronoi cell computation, these routines are used to create
 * the regular Voronoi tessellation. */
class radius_mono {
    protected:
        /** This is called prior to computing a Voronoi cell for a given
         * particle to initialize any required constants.
         * \param[in] ijk the block that the particle is within.
         * \param[in] s the index of the particle within the block. */
        inline void r_init(int ijk,int s,double &r_rad,double &r_mul) {}
        /** Sets a required constant to be used when carrying out a plane
         * bounds check. */
        inline void r_prime(double rv,double &r_mul,double &r_val) {}
        /** Carries out a radius bounds check.
         * \param[in] crs the radius squared to be tested.
         * \param[in] mrs the current maximum distance to a Voronoi
         *                vertex multiplied by two.
         * \return True if particles at this radius could not possibly
         * cut the cell, false otherwise. */
        inline bool r_ctest(double crs,double mrs,double &r_mul) {return crs>mrs;}
        /** Scales a plane displacement during a plane bounds check.
         * \param[in] lrs the plane displacement.
         * \return The scaled value. */
        inline double r_cutoff(double lrs,double &r_val) {return lrs;}
        /** Adds the maximum radius squared to a given value.
         * \param[in] rs the value to consider.
         * \return The value with the radius squared added. */
        inline double r_max_add(double rs) {return rs;}
        /** Subtracts the radius squared of a particle from a given value.
         * \param[in] rs the value to consider.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The value with the radius squared subtracted. */
        inline double r_current_sub(double rs,int ijk,int q) {return rs;}
        /** Scales a plane displacement prior to use in the plane cutting
         * algorithm.
         * \param[in] rs the initial plane displacement.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The scaled plane displacement. */
        inline double r_scale(double rs,int ijk,int q,double &r_rad) {return rs;}
        /** Scales a plane displacement prior to use in the plane cutting
         * algorithm, and also checks if it could possibly cut the cell.
         * \param[in,out] rs the plane displacement to be scaled.
         * \param[in] mrs the current maximum distance to a Voronoi vertex
         *                multiplied by two.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return True if the cell could possibly cut the cell, false
         * otherwise. */
        inline bool r_scale_check(double &rs,double mrs,int ijk,int q,double &r_rad) {return rs<mrs;}
};

/**  \brief Class containing all of the routines that are specific to computing
 * the radical Voronoi tessellation.
 *
 * The container_poly and container_periodic_poly classes are derived from this
 * class, and during the Voronoi cell computation, these routines are used to
 * create the radical Voronoi tessellation. */
class radius_poly_3d {
    public:
        /** A two-dimensional array holding particle positions and radii. */
        double **ppr;
        /** The current maximum radius of any particle, used to determine when
         * to cut off the radical Voronoi computation. */
        double max_radius;
        /** The class constructor sets the maximum particle radius to be zero.
         */
        radius_poly_3d() : max_radius(0) {}
    protected:
        /** This is called prior to computing a Voronoi cell for a given
         * particle to initialize any required constants.
         * \param[in] ijk the block that the particle is within.
         * \param[in] s the index of the particle within the block. */
        inline void r_init(int ijk,int s,double &r_rad,double &r_mul) {
            r_rad=ppr[ijk][4*s+3]*ppr[ijk][4*s+3];
            r_mul=r_rad-max_radius*max_radius;
        }
        /** Sets a required constant to be used when carrying out a plane
         * bounds check. */
        inline void r_prime(double rv,double &r_mul,double &r_val) {
            r_val=1+r_mul/rv;
        }
        /** Carries out a radius bounds check.
         * \param[in] crs the radius squared to be tested.
         * \param[in] mrs the current maximum distance to a Voronoi vertex
         *                multiplied by two.
         * \return True if particles at this radius could not possibly
         * cut the cell, false otherwise. */
        inline bool r_ctest(double crs,double mrs,double &r_mul) {
            return crs+r_mul>sqrt(mrs*crs);
        }
        /** Scales a plane displacement during a plane bounds check.
         * \param[in] lrs the plane displacement.
         * \return The scaled value. */
        inline double r_cutoff(double lrs,double &r_val) {
            return lrs*r_val;
        }
        /** Adds the maximum radius squared to a given value.
         * \param[in] rs the value to consider.
         * \return The value with the radius squared added. */
        inline double r_max_add(double rs) {return rs+max_radius*max_radius;}
        /** Subtracts the radius squared of a particle from a given value.
         * \param[in] rs the value to consider.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The value with the radius squared subtracted. */
        inline double r_current_sub(double rs,int ijk,int q) {
            return rs-ppr[ijk][4*q+3]*ppr[ijk][4*q+3];
        }
        /** Scales a plane displacement prior to use in the plane cutting
         * algorithm.
         * \param[in] rs the initial plane displacement.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The scaled plane displacement. */
        inline double r_scale(double rs,int ijk,int q,double &r_rad) {
            return rs+r_rad-ppr[ijk][4*q+3]*ppr[ijk][4*q+3];
        }

        /** Scales a plane displacement prior to use in the plane cutting
         * algorithm, and also checks if it could possibly cut the cell.
         * \param[in,out] rs the plane displacement to be scaled.
         * \param[in] mrs the current maximum distance to a Voronoi vertex
         *                multiplied by two.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return True if the cell could possibly cut the cell, false
         * otherwise. */
        inline bool r_scale_check(double &rs,double mrs,int ijk,int q,double &r_rad) {
            double trs=rs;
            rs+=r_rad-ppr[ijk][4*q+3]*ppr[ijk][4*q+3];
            return rs<sqrt(mrs*trs);
        }

};

class radius_poly_2d {
    public:
        /** A two-dimensional array holding particle positions and radii. */
        double **ppr;
        /** The current maximum radius of any particle, used to
         * determine when to cut off the radical Voronoi computation.
         * */
        double max_radius;
        /** The class constructor sets the maximum particle radius to
         * be zero. */
        radius_poly_2d() : max_radius(0) {}
    protected:
        /** This is called prior to computing a Voronoi cell for a
         * given particle to initialize any required constants.
         * \param[in] ijk the block that the particle is within.
         * \param[in] s the index of the particle within the block. */
        inline void r_init(int ijk,int s,double &r_rad,double &r_mul) {
            r_rad=ppr[ijk][3*s+2]*ppr[ijk][3*s+2];
            r_mul=r_rad-max_radius*max_radius;
        }
        /** Sets a required constant to be used when carrying out a
         * plane bounds check. */
        inline void r_prime(double rv,double &r_mul,double &r_val) {
            r_val=1+r_mul/rv;
        }
        /** Carries out a radius bounds check.
         * \param[in] crs the radius squared to be tested.
         * \param[in] mrs the current maximum distance to a Voronoi
         *                vertex multiplied by two.
         * \return True if particles at this radius could not possibly
         * cut the cell, false otherwise. */
        inline bool r_ctest(double crs,double mrs,double &r_mul) {
            return crs+r_mul>sqrt(mrs*crs);
        }
        /** Scales a plane displacement during a plane bounds check.
         * \param[in] lrs the plane displacement.
         * \return The scaled value. */
        inline double r_cutoff(double lrs,double &r_val) {
            return lrs*r_val;
        }
        /** Adds the maximum radius squared to a given value.
         * \param[in] rs the value to consider.
         * \return The value with the radius squared added. */
        inline double r_max_add(double rs) {return rs+max_radius*max_radius;}
        /** Subtracts the radius squared of a particle from a given
         * value.
         * \param[in] rs the value to consider.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The value with the radius squared subtracted. */
        inline double r_current_sub(double rs,int ijk,int q) {
            return rs-ppr[ijk][3*q+2]*ppr[ijk][3*q+2];
        }
        /** Scales a plane displacement prior to use in the plane cutting
         * algorithm.
         * \param[in] rs the initial plane displacement.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return The scaled plane displacement. */
        inline double r_scale(double rs,int ijk,int q,double &r_rad) {
            return rs+r_rad-ppr[ijk][3*q+2]*ppr[ijk][3*q+2];
        }
        /** Scales a plane displacement prior to use in the plane
         * cutting algorithm, and also checks if it could possibly cut
         * the cell.
         * \param[in,out] rs the plane displacement to be scaled.
         * \param[in] mrs the current maximum distance to a Voronoi
         *                vertex multiplied by two.
         * \param[in] ijk the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return True if the cell could possibly cut the cell, false
         * otherwise. */
        inline bool r_scale_check(double &rs,double mrs,int ijk,int q,double &r_rad) {
            double trs=rs;
            rs+=r_rad-ppr[ijk][3*q+2]*ppr[ijk][3*q+2];
            return rs<sqrt(mrs*trs);
        }
};

}
#endif
