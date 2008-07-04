// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#ifndef VOROPP_CONFIG_HH
#define VOROPP_CONFIG_HH

// These constants set the initial memory allocation for the Voronoi cell
const int init_vertices=256;
const int init_vertex_order=64;
const int init_3_vertices=256;
const int init_n_vertices=8;
const int init_marginal=256;
const int init_delete_size=256;
const int init_delete2_size=256;
const int init_facet_size=32;
const int init_wall_size=32;

// If the initial memory is too small, the program dynamically allocates more.
// However, if the limits below are reached, then the program bails out.
const int max_vertices=16777216;
const int max_vertex_order=2048;
const int max_n_vertices=16777216;
const int max_marginal=16777216;
const int max_delete_size=16777216;
const int max_delete2_size=16777216;
const int max_particle_memory=16777216;
const int max_wall_size=2048;

// Voro++ can print a number of different status and debugging messages to
// notify the user of special behavior, and this macro sets the amount which
// are displayed. At level 0, no messages are printed. At level 1, messages
// about unusual cases during cell construction are printed, such as when the
// plane routine bails out due to floating point problems. At level 2, general
// messages about memory expansion are printed. At level 3, technical details
// about memory management are printed.
#define VOROPP_VERBOSE 0

// This sets the numerical tolerance. Below these values, the plane cutting
// snaps to existing vertices rather than creating new ones.
#ifdef VOROPP_SINGLE_PRECISION
typedef float fpoint;
const fpoint tolerance=1e-5;
const fpoint tolerance2=2e-5;
#else
typedef double fpoint;
const fpoint tolerance=1e-10;
const fpoint tolerance2=2e-10;
#endif
const fpoint large_number=1e30;

#endif
