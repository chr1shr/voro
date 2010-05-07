/** \file config.hh
 * \brief Master configuration file for setting various compile-time options. */

#ifndef VOROPP_CONFIG_HH
#define VOROPP_CONFIG_HH

// These constants set the initial memory allocation for the Voronoi cell
/** The initial memory allocation for the number of vertices. */
const int init_vertices=256;
/** The initial size for the delete stack. */
const int init_delete_size=256;

// If the initial memory is too small, the program dynamically allocates more.
// However, if the limits below are reached, then the program bails out.
/** The maximum memory allocation for the number of vertices. */
const int max_vertices=16777216;
/** The maximum size for the delete stack. */
const int max_delete_size=16777216;

/** The maximum amount of particle memory allocated for a single region. */
const int max_particle_memory=16777216;

#ifndef VOROPP_VERBOSE
/** Voro++ can print a number of different status and debugging messages to
 * notify the user of special behavior, and this macro sets the amount which
 * are displayed. At level 0, no messages are printed. At level 1, messages
 * about unusual cases during cell construction are printed, such as when the
 * plane routine bails out due to floating point problems. At level 2, general
 * messages about memory expansion are printed. At level 3, technical details
 * about memory management are printed. */
#define VOROPP_VERBOSE 0
#endif

/** The declaration of fpoint allows that code to be compiled both using single
 * precision numbers and double precision numbers. Under normal usage fpoint is
 * set be a double precision floating point number, but defining the
 * preprocessor macro VOROPP_SINGLE_PRECISION will switch it to single
 * precision and make the code tolerances larger. */
#ifdef VOROPP_SINGLE_PRECISION
typedef float fpoint;
#else
typedef double fpoint;
#endif

/** If a point is within this distance of a cutting plane, then the code
 * assumes that point exactly lies on the plane. */
#ifdef VOROPP_SINGLE_PRECISION
const fpoint tolerance=1e-5;
#else
const fpoint tolerance=1e-10;
#endif

/** Voro++ returns this status code if there is a file-related error, such as
 * not being able to open file. */
#define VOROPP_FILE_ERROR 1

/** Voro++ returns this status code if there is a memory allocation error, if
 * one of the safe memory limits is exceeded. */
#define VOROPP_MEMORY_ERROR 2

/** Voro++ returns this status code if there is any type of internal error, if
 * it detects that representation of the Voronoi cell is inconsistent. This
 * status code will generally indicate a bug, and the developer should be
 * contacted. */
#define VOROPP_INTERNAL_ERROR 3

/** Voro++ returns this status code if it could not interpret the command line
 * arguments passed to the command line utility. */
#define VOROPP_CMD_LINE_ERROR 4

#endif
