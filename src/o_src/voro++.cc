// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file voro++.cc
 * \brief A file that loads all of the function implementation files. */

#include "cell.cc"
#include "container.cc"
#include "wall.cc"

template voropp_loop::voropp_loop(container_base<radius_mono> *);
template voropp_loop::voropp_loop(container_base<radius_poly> *);
template class container_base<radius_mono>;
template class container_base<radius_poly>;
template class voronoicell_base<neighbor_none>;
template class voronoicell_base<neighbor_track>;

template bool wall_cylinder::cut_cell_base(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z);
template bool wall_cylinder::cut_cell_base(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z);
template bool wall_plane::cut_cell_base(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z);
template bool wall_plane::cut_cell_base(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z);
template bool wall_sphere::cut_cell_base(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z);
template bool wall_sphere::cut_cell_base(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z);
template bool wall_cone::cut_cell_base(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z);
template bool wall_cone::cut_cell_base(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z);
