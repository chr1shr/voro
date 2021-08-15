// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file wall_2d.cc
 * \brief Function implementations for the 2D derived wall classes. */

#include "wall_2d.hh"

namespace voro {

/** Tests to see whether a point is inside the sphere wall object.
 * \param[in,out] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_circle::point_inside(double x,double y) {
    return (x-xc)*(x-xc)+(y-yc)*(y-yc)<rc*rc;
}

/** Cuts a cell by the circular wall object. The circular wall is approximated
 * by a single plane applied at the point on the sphere which is closest to the
 * center of the cell. This works well for particle arrangements that are
 * packed against the wall, but loses accuracy for sparse particle
 * distributions.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
template<class v_cell_2d>
bool wall_circle::cut_cell_base(v_cell_2d &c,double x,double y) {
    double xd=x-xc,yd=y-yc,dq;
    dq=xd*xd+yd*yd;
    if(dq>1e-5) {
        dq=2*(sqrt(dq)*rc-dq);
        return c.nplane(xd,yd,dq,w_id);
    }
    return true;
}

/** Tests to see whether a point is inside the line wall object.
 * \param[in] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_line::point_inside(double x,double y) {
    return x*xc+y*yc<ac;
}

/** Cuts a cell by the line wall object.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
template<class v_cell_2d>
bool wall_line::cut_cell_base(v_cell_2d &c,double x,double y) {
    double dq=2*(ac-x*xc-y*yc);
    return c.nplane(xc,yc,dq,w_id);
}

// Explicit instantiation
template bool wall_circle::cut_cell_base(voronoicell_2d &c,double x,double y);
template bool wall_circle::cut_cell_base(voronoicell_neighbor_2d &c,double x,double y);
template bool wall_line::cut_cell_base(voronoicell_2d &c,double x,double y);
template bool wall_line::cut_cell_base(voronoicell_neighbor_2d &c,double x,double y);

}
