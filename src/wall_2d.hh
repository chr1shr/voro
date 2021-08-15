// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file wall_2d.hh
 * \brief Header file for the 2D derived wall classes. */

#ifndef VOROPP_WALL_2D_HH
#define VOROPP_WALL_2D_HH

#include "cell_2d.hh"
#include "container_2d.hh"

namespace voro {

/** \brief A class representing a circular wall object.
 *
 * This class represents a circular wall object. */
struct wall_circle : public wall_2d {
    public:
        /** Constructs a spherical wall object.
         * \param[in] w_id_ an ID number to associate with the wall for
         *                  neighbor tracking.
         * \param[in] (xc_,yc_) a position vector for the circle's
         *                      center.
         * \param[in] rc_ the radius of the circle.
         * \param[in] w_id_ an ID number to associate with the wall for
         *                  neighbor tracking. */
        wall_circle(double xc_,double yc_,double rc_,int w_id_=-99)
            : w_id(w_id_), xc(xc_), yc(yc_), rc(rc_) {}
        bool point_inside(double x,double y);
        template<class v_cell_2d>
        bool cut_cell_base(v_cell_2d &c,double x,double y);
        bool cut_cell(voronoicell_2d &c,double x,double y) {return cut_cell_base(c,x,y);}
        bool cut_cell(voronoicell_neighbor_2d &c,double x,double y) {return cut_cell_base(c,x,y);}
    private:
        const int w_id;
        const double xc,yc,rc;
};

/** \brief A class representing a line wall object.
 *
 * This class represents a single line wall object. */
struct wall_line : public wall_2d {
    public:
        /** Constructs a line wall object
         * \param[in] (xc_,yc_) a normal vector to the line.
         * \param[in] ac_ a displacement along the normal vector.
         * \param[in] w_id_ an ID number to associate with the wall for
         *                  neighbor tracking. */
        wall_line(double xc_,double yc_,double ac_,int w_id_=-99)
            : w_id(w_id_), xc(xc_), yc(yc_), ac(ac_) {}
        bool point_inside(double x,double y);
        template<class v_cell_2d>
        bool cut_cell_base(v_cell_2d &c,double x,double y);
        bool cut_cell(voronoicell_2d &c,double x,double y) {return cut_cell_base(c,x,y);}
        bool cut_cell(voronoicell_neighbor_2d &c,double x,double y) {return cut_cell_base(c,x,y);}
    private:
        const int w_id;
        const double xc,yc,ac;
};

}

#endif
