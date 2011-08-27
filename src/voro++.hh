// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 28th 2011

/** \file voro++.hh
 * \brief A file that loads all of the Voro++ header files. */

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
 * represents the cell as a collection of vertices that are connected by edges,
 * and there are routines for initializing, making, and outputting the cell.
 * The container class represents a three-dimensional simulation region into
 * which particles can be added. The class can then carry out a variety of
 * Voronoi calculations by computing cells using the voronoicell class. It also
 * has a general mechanism using virtual functions to implement walls,
 * discussed in more detail below. To implement the radical Voronoi
 * tessellation and the neighbor calculations, two class variants called
 * voronoicell_neighbor and container_poly are provided by making use of
 * templates, which is discussed below.
 *
 * \section voronoicell The voronoicell class
 * The voronoicell class represents a single Voronoi cell as a convex
 * polyhedron, with a set of vertices that are connected by edges. The class
 * contains a variety of functions that can be used to compute and output the
 * Voronoi cell corresponding to a particular particle. The command init()
 * can be used to initialize a cell as a large rectangular box. The Voronoi cell
 * can then be computed by repeatedly cutting it with planes that correspond to
 * the perpendicular bisectors between that particle and its neighbors.
 *
 * This is achieved by using the plane() routine, which will recompute the
 * cell's vertices and edges after cutting it with a single plane. This is the
 * key routine in voronoicell class. It begins by exploiting the convexity
 * of the underlying cell, tracing between edges to work out if the cell
 * intersects the cutting plane. If it does not intersect, then the routine
 * immediately exits. Otherwise, it finds an edge or vertex that intersects
 * the plane, and from there, traces out a new face on the cell, recomputing
 * the edge and vertex structure accordingly.
 *
 * Once the cell is computed, it can be drawn using commands such as
 * draw_gnuplot() and draw_pov(), or its volume can be evaluated using the
 * volume() function. Many more routines are available, and are described in
 * the online reference manual.
 *
 * \subsection internal Internal data representation
 * The voronoicell class has a public member p representing the
 * number of vertices. The polyhedral structure of the cell is stored
 * in the following arrays:
 *
 * - pts[]: an array of floating point numbers, that represent the position
 *   vectors x_0, x_1, ..., x_{p-1} of the polyhedron vertices.
 * - nu[]: the order of each vertex n_0, n_1, ..., n_{p-1}, corresponding to
 *   the number of other vertices to which each is connected.
 * - ed[][]: a table of edges and relations. For the ith vertex, ed[i] has
 *   2n_i+1 elements. The first n_i elements are the edges e(j,i), where e(j,i)
 *   is the jth neighbor of vertex i. The edges are ordered according to a
 *   right-hand rule with respect to an outward-pointing normal. The next n_i
 *   elements are the relations l(j,i) which satisfy the property
 *   e(l(j,i),e(j,i)) = i. The final element of the ed[i] list is a back
 *   pointer used in memory allocation.
 *
 * In a very large number of cases, the values of n_i will be 3. This is because
 * the only way that a higher-order vertex can be created in the plane()
 * routine is if the cutting plane perfectly intersects an existing vertex. For
 * random particle arrangements with position vectors specified to double
 * precision this should happen very rarely. A preliminary version of this code
 * was quite successful with only making use of vertices of order 3. However,
 * when calculating millions of cells, it was found that this approach is not
 * robust, since a single floating point error can invalidate the computation.
 * This can also be a problem for cases featuring crystalline arrangements of
 * particles where the corresponding Voronoi cells may have high-order vertices
 * by construction.
 *
 * Because of this, Voro++ takes the approach that it if an existing vertex is
 * within a small numerical tolerance of the cutting plane, it is treated as
 * being exactly on the plane, and the polyhedral topology is recomputed
 * accordingly. However, while this improves robustness, it also adds the
 * complexity that n_i may no longer always be 3. This causes memory management
 * to be significantly more complicated, as different vertices require a
 * different number of elements in the ed[][] array. To accommodate this, the
 * voronoicell class allocated edge memory in a different array called mep[][],
 * in such a way that all vertices of order k are held in mep[k]. If vertex
 * i has order k, then ed[i] points to memory within mep[k]. The array ed[][]
 * is never directly initialized as a two-dimensional array itself, but points
 * at allocations within mep[][]. To the user, it appears as though each row of
 * ed[][] has a different number of elements. When vertices are added or
 * deleted, care must be taken to reorder and reassign elements in these
 * arrays.
 *
 * During the plane() routine, the code traces around the vertices of the cell,
 * and adds new vertices along edges which intersect the cutting plane to
 * create a new face. The values of l(j,i) are used in this computation, as
 * when the code is traversing from one vertex on the cell to another, this
 * information allows the code to immediately work out which edge of a vertex
 * points back to the one it came from. As new vertices are created, the
 * l(j,i) are also updated to ensure consistency. To ensure robustness, the
 * plane cutting algorithm should work with any possible combination of
 * vertices which are inside, outside, or exactly on the cutting plane.
 *
 * Vertices exactly on the cutting plane create some additional computational
 * difficulties. If there are two marginal vertices connected by an existing
 * edge, then it would be possible for duplicate edges to be created between
 * those two vertices, if the plane routine traces along both sides of this
 * edge while constructing the new face. The code recognizes these cases and
 * prevents the double edge from being formed. Another possibility is the
 * formation of vertices of order two or one. At the end of the plane cutting
 * routine, the code checks to see if any of these are present, removing the
 * order one vertices by just deleting them, and removing the order two
 * vertices by connecting the two neighbors of each vertex together. It is
 * possible that the removal of a single low-order vertex could result in the
 * creation of additional low-order vertices, so the process is applied
 * recursively until no more are left.
 *
 * \section container The container class
 * The container class represents a three-dimensional rectangular box of
 * particles. The constructor for this class sets up the coordinate ranges,
 * sets whether each direction is periodic or not, and divides the box into a
 * rectangular subgrid of regions. Particles can be added to the container
 * using the put() command, that adds a particle's position and an integer
 * numerical ID label to the corresponding region. Alternatively, the command
 * import() can be used to read large numbers of particles from a text file.
 *
 * The key routine in this class is compute_cell(), which makes use of the
 * voronoicell class to construct a Voronoi cell for a specific particle in the
 * container. The basic approach that this function takes is to repeatedly cut
 * the Voronoi cell by planes corresponding neighboring particles, and stop
 * when it recognizes that all the remaining particles in the container are too
 * far away to possibly influence cell's shape. The code makes use of two
 * possible methods for working out when a cell computation is complete:
 *
 * - Radius test: if the maximum distance of a Voronoi cell
 *   vertex from the cell center is R, then no particles more than a distance
 *   2R away can possibly influence the cell. This a very fast computation to
 *   do, but it has no directionality: if the cell extends a long way in one
 *   direction then particles a long distance in other directions will still
 *   need to be tested.
 * - Region test: it is possible to test whether a specific region can
 *   possibly influence the cell by applying a series of plane tests at the
 *   point on the region which is closest to the Voronoi cell center. This is a
 *   slower computation to do, but it has directionality.
 *
 * Another useful observation is that the regions that need to be tested
 * are simply connected, meaning that if a particular region does not need
 * to be tested, then neighboring regions which are further away do not
 * need to be tested.
 *
 * For maximum efficiency, it was found that a hybrid approach making use of both
 * of the above tests worked well in practice. Radius tests work well for the
 * first few blocks, but switching to region tests after then prevent the code
 * from becoming extremely slow, due to testing over very large spherical shells of
 * particles. The compute_cell() routine therefore takes the following
 * approach:
 *
 * - Initialize the voronoicell class to fill the entire computational domain.
 * - Cut the cell by any wall objects that have been added to the container.
 * - Apply plane cuts to the cell corresponding to the other particles which
 *   are within the current particle's region.
 * - Test over a pre-computed worklist of neighboring regions, that have been
 *   ordered according to the minimum distance away from the particle's
 *   position. Apply radius tests after every few regions to see if the
 *   calculation can terminate.
 * - If the code reaches the end of the worklist, add all the neighboring
 *   regions to a new list.
 * - Carry out a region test on the first item of the list. If the region needs
 *   to be tested, apply the plane() routine for all of its particles, and then
 *   add any neighboring regions to the end of the list that need to be tested.
 *   Continue until the list has no elements left.
 *
 * The compute_cell() routine forms the basis of many other routines, such as
 * store_cell_volumes() and draw_cells_gnuplot() that can be used to calculate
 * and draw the cells in the entire container or in a subdomain.
 *
 * \section walls Wall computation
 * Wall computations are handled by making use of a pure virtual wall class.
 * Specific wall types are derived from this class, and require the
 * specification of two routines: point_inside() that tests to see if a point
 * is inside a wall or not, and cut_cell() that cuts a cell according to the
 * wall's position. The walls can be added to the container using the
 * add_wall() command, and these are called each time a compute_cell() command
 * is carried out. At present, wall types for planes, spheres, cylinders, and
 * cones are provided, although custom walls can be added by creating new
 * classes derived from the pure virtual class. Currently all wall types
 * approximate the wall surface with a single plane, which produces some small
 * errors, but generally gives good results for dense particle packings in
 * direct contact with a wall surface. It would be possible to create more
 * accurate walls by making cut_cell() routines that approximate the curved
 * surface with multiple plane cuts.
 *
 * The wall objects can used for periodic calculations, although to obtain
 * valid results, the walls should also be periodic as well. For example, in a
 * domain that is periodic in the x direction, a cylinder aligned along the x
 * axis could be added. At present, the interior of all wall objects are convex
 * domains, and consequently any superposition of them will be a convex domain
 * also. Carrying out computations in non-convex domains poses some problems,
 * since this could theoretically lead to non-convex Voronoi cells, which the
 * internal data representation of the voronoicell class does not support. For
 * non-convex cases where the wall surfaces feature just a small amount of
 * negative curvature (eg. a torus) approximating the curved surface with a
 * single plane cut may give an acceptable level of accuracy. For non-convex
 * cases that feature internal angles, the best strategy may be to decompose
 * the domain into several convex subdomains, carry out a calculation in each,
 * and then add the results together. The voronoicell class cannot be easily
 * modified to handle non-convex cells as this would fundamentally alter the
 * algorithms that it uses, and cases could arise where a single plane cut
 * could create several new faces as opposed to just one.
 *
 * \section templates Extra functionality via the use of templates
 * C++ templates are often presented as a mechanism for allowing functions to
 * be coded to work with several different data types. However, they also
 * provide an extremely powerful mechanism for achieving static polymorphism,
 * allowing several variations of a program to be compiled from a single source
 * code. Voro++ makes use of templates in order to handle the radical Voronoi
 * tessellation and the neighbor calculations, both of which require only
 * relatively minimal alterations to the main body of code.
 *
 * The main body of the voronoicell class is written as template named
 * voronoicell_base. Two additional small classes are then written:
 * neighbor_track, which contains small, inlined functions that encapsulate all
 * of the neighbor calculations, and neighbor_none, which contains the same
 * function names left blank. By making use of the typedef command, two classes
 * are then created from the template:
 *
 * - voronoicell: an instance of voronoicell_base with the neighbor_none class.
 * - voronoicell_neighbor: an instance of voronoicell_base with the
 *   neighbor_track class.
 *
 * The two classes will be the same, except that the second will get all of the
 * additional neighbor-tracking functionality compiled into it through the
 * neighbor_track class. Since the two instances of the template are created
 * during the compilation, and since all of the functions in neighbor_none and
 * neighbor_track are inlined, there should be no speed overhead with this
 * construction; it should have the same efficiency as writing two completely
 * separate classes. C++ has other methods for achieving similar results, such
 * as virtual functions and class inheritance, but these are more focused on
 * dynamic polymorphism, switching between functionality at run-time, resulting
 * in a drop in performance. This would be particularly apparent in this case,
 * as the neighbor computation code, while small, is heavily integrated into
 * the low-level details of the plane() routine, and a virtual function
 * approach would require a very large number of function address look-ups.
 *
 * In a similar manner, two small classes called radius_mono and radius_poly
 * are provided. The first contains all routines suitable for calculate the
 * standard Voronoi tessellation associated with a monodisperse particle
 * packing, while the second incorporates variations to carry out the radical
 * Voronoi tessellation associated with a polydisperse particle packing. Two
 * classes are then created via typedef commands:
 *
 * - container: an instance of container_base with the radius_mono class.
 * - container_poly: an instance of container_base with the radius_poly class.
 *
 * The container_poly class accepts an additional variable in the put() command
 * for the particle's radius. These radii are then used to weight the plane
 * positions in the compute_cell() routine.
 *
 * It should be noted that the underlying template structure is largely hidden
 * from a typical user accessing the library's functionality, and as
 * demonstrated in the examples, the classes listed above behave like regular
 * C++ classes, and can be used in all the same ways. However, the template
 * structure may provide an additional method of customizing the code; for
 * example, an additional radius class could be written to implement a Voronoi
 * tessellation variant. */

#ifndef VOROPP_HH
#define VOROPP_HH

#include "config.hh"
#include "common.hh"
#include "cell.hh"
#include "v_base.hh"
#include "container.hh"
#include "unitcell.hh"
#include "container_prd.hh"
#include "pre_container.hh"
#include "v_compute.hh"
#include "v_loops.hh"
#include "wall.hh"

#endif
