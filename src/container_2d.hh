/** \file container_2d.hh
 * \brief Header file for the container_2d class. */

#ifndef VOROPP_CONTAINER_2D_HH
#define VOROPP_CONTAINER_2D_HH

#include "config.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class voropp_loop_2d;

/** \brief A class representing the whole 2D simulation region.
 *
 * The container class represents the whole simulation region. The
 * container constructor sets up the geometry and periodicity, and divides
 * the geometry into rectangular grid of blocks, each of which handles the
 * particles in a particular area. Routines exist for putting in particles,
 * importing particles from standard input, and carrying out Voronoi
 * calculations. */
class container_2d {
	public:
		container_2d(fpoint xa,fpoint xb,fpoint ya,fpoint yb,int xn,int yn,bool xper,bool yper,int memi);
		~container_2d();
		void draw_particles(const char *filename);
		void draw_particles();
		void draw_particles(ostream &os);
		void import(istream &is);
		inline void import();
		inline void import(const char *filename);
		void draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax);
		inline void draw_cells_gnuplot(const char *filename);
		inline bool compute_cell_sphere(voronoicell_2d &c,int i,int j,int ij,int s);
		bool compute_cell_sphere(voronoicell_2d	&c,int i,int j,int ij,int s,fpoint x,fpoint y);
		bool initialize_voronoicell(voronoicell_2d &c,fpoint x,fpoint y);
		void put(int n,fpoint x,fpoint y);
		void clear();
	private:
		void add_particle_memory(int i);
		/** The minimum x coordinate of the container. */
		const fpoint ax;
		/** The maximum x coordinate of the container. */
		const fpoint bx;
		/** The minimum y coordinate of the container. */
		const fpoint ay;
		/** The maximum y coordinate of the container. */
		const fpoint by;
		/** The inverse box length in the x direction, set to
		 * nx/(bx-ax). */
		const fpoint xsp;
		/** The inverse box length in the y direction, set to
		 * ny/(by-ay). */
		const fpoint ysp;
		/** The number of boxes in the x direction. */
		const int nx;
		/** The number of boxes in the y direction. */
		const int ny;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines which step through boxes in
		 * sequence. */
		const int nxy;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** This array holds the number of particles within each
		 * computational box of the container. */
		int *co;
		/** This array holds the maximum amount of particle memory for
		 * each computational box of the container. If the number of
		 * particles in a particular box ever approaches this limit,
		 * more is allocated using the add_particle_memory() function.
		 */
		int *mem;
		/** This array holds the numerical IDs of each particle in each
		 * computational box. */
		int **id;
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		fpoint **p;
		friend class voropp_loop_2d;	
};


/** \brief A class to handle loops on regions of the container handling
 * non-periodic and periodic boundary conditions.
 *
 * Many of the container routines require scanning over a rectangular sub-grid
 * of blocks, and the routines for handling this are stored in the
 * voropp_loop_2d class. A voropp_loop_2d class can first be initialized to
 * either calculate the subgrid which is within a distance r of a vector
 * (vx,vy,vz), or a subgrid corresponding to a rectangular box. The routine
 * inc() can then be successively called to step through all the blocks within
 * this subgrid.
 */
class voropp_loop_2d {
	public:
		voropp_loop_2d(container_2d *q);
		inline int init(fpoint vx,fpoint vy,fpoint r,fpoint &px,fpoint &py);
		inline int init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint &px,fpoint &py);
		inline int inc(fpoint &px,fpoint &py);
		/** The current block index in the x direction, referencing a
		 * real cell in the range 0 to nx-1. */
		int ip;
		/** The current block index in the y direction, referencing a
		 * real cell in the range 0 to ny-1. */
		int jp;
	private:
		int i,j,ai,bi,aj,bj,s;
		int aip,ajp,inc1;
		inline int step_mod(int a,int b);
		inline int step_div(int a,int b);
		inline int step_int(fpoint a);
		fpoint apx,apy;
		const fpoint sx,sy,xsp,ysp,ax,ay;
		const int nx,ny,nxy;
		const bool xperiodic,yperiodic;
};

#endif
