// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \mainpage Voro++ class reference manual
 * \section intro Introduction
 * Voro++ is a software library for carrying out 3D cell-based calculations of
 * the Voronoi tessellation. It is primarily designed for applications in
 * physics and materials science, where the Voronoi tessellation can be a
 * useful tool in analyzing particle systems.
 *
 * Voro++ is comprised of several C++ classes, and is designed to be
 * incorporated into other programs. This manual provides a reference for every
 * function in the class structure. For a general overview of the program, see
 * the Voro++ website at http://math.lbl.gov/voro++/ and in particular the
 * example programs at http://math.lbl.gov/voro++/examples/ that demonstrate
 * many of the library's features.
 * 
 * \section class C++ class structure 
 * The code is structured around two main C++ classes. The voronoicell class
 * contains all of the routines for constructing a single Voronoi cell. It
 * represents the cells as a collection of vertices that are connected by
 * edges. The command init() can be used to initialize a cell as a large
 * rectangular box. The Voronoi cell can then be computed by repeatedly cutting
 * it with planes that correspond to the perpendicular bisectors between on
 * point and its neighbors. The command plane() will recompute the cell after
 * cutting with a single plane. Once the cell is computed, it can be drawn
 * using commands such as draw_gnuplot() and draw_pov(), or its volume can be
 * computed using the volume() command.
 * 
 * The container class represents a 3D rectangular box of particles. The
 * constructor for this class sets up the coordinate ranges, sets whether each
 * direction is periodic or not, and divides the box into a rectangular subgrid
 * of regions. Particles can be added to the container using the put() command,
 * that adds a particle's position and an integer numerical ID label, to the
 * container by adding it to the corresponding region. The command import() can
 * be used to read in large numbers of particles from a text file. The
 * compute_cell() function creates a single Voronoi cell for a particle in the
 * container, by making use of the voronoicell class, and constructing the cell
 * by making calls to the plane() routine to account for neighboring points.
 * Various commands such as store_cell_volumes() and draw_cells_gnuplot() can
 * be used to calculate and draw the cells in the entire container or in a
 * subregion.
 *
 * \section templates Extra functionality and the use of templates
 * The library also supports the Voronoi radical tessellation, by using the
 * container_poly class that is a variant of the container class where the
 * put() command accepts an additional argument for the particle radius. To
 * create this without repeating large parts of the code, the library makes use
 * of templates. The source code is structured around a general template called
 * container_base. There are then two small classes called radius_mono and
 * radius_poly that handle all of the routines for the regular and radical
 * tessellations respectively. The container class is created as the
 * instantiation of the container_base template with the radius_mono class, and
 * the container_poly class is the instantiation of the container_base template
 * with the radius_poly class. Since the different instances are created during
 * the program compilation and since all of the functions of radius_mono and
 * radius_poly are declared inline, this approach should result in minimal
 * overhead with a good compiler.
 *
 * Similarly, the voronoicell class is constructed around a template called
 * voronoicell_base. The voronoicell class is an instantiation of this template
 * using the neighbor_none class, that does not compute any neighbor
 * information. A variant called voronoicell_neighbor is also available, that
 * makes use of the neighbor_track class to additionally carry out neighbor
 * tracking during the cell construction.
 *
 * \section walls Wall computation 
 * Wall computations are handled by making use of a pure virtual wall class.
 * Specific wall types are derived from this class, and have a routine called
 * point_inside() that tests to see if a point is inside a wall or not, and a
 * routine called cut_cell() that cuts a cell according to the wall's position.
 * The walls can be added to the container using the add_wall() command, and
 * these are called each time a compute_cell() command is carried out. At
 * present, wall types for planes, spheres, cylinders, and cones are provided,
 * although custom walls can be added by creating new classes derived from the
 * pure virtual class. */

#ifndef VOROPP_HH
#define VOROPP_HH

#include "cell.hh"
#include "container.hh"
#include "wall.hh"

#endif
