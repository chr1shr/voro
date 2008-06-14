// Voronoi calculation header file
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#ifndef FACETS_CONTAINER_HH
#define FACETS_CONTAINER_HH

#include "config.hh"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class facets_loop;
class facets_search;

/** The container class represents the whole simulation region. The
 * container constructor sets up the geometry and periodicity, and divides
 * the geometry into rectangular grid of blocks, each of which handles the
 * particles in a particular area. Routines exist for putting in particles,
 * importing particles from standard input, and carrying out Voronoi
 * calculations. */
class container {
	public:
		/** This represents a typical length scale for the diameter of
		 * the particles. It initially takes the default value of 1,
		 * but its value can be overwritten by the user. There is also
		 * a function guess_length_scale(), which picks a value based
		 * on the total number of particles in the container.
		 *
		 * The Voronoi cell calculation works by evaluating all
		 * particles in concentric spherical shells. The length scale
		 * is used to set how big these shells should be. The closer
		 * the length scale to the actual separation between particles,
		 * the faster the cell computation should be. */
		fpoint length_scale;
		container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi,int isz);
		container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi);
		~container();
		void dump(char *filename);
		void import(istream &is);
		inline void import();
		inline void import(char *filename);
		void region_count();
		void clear();
		void draw_gnuplot(char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void draw_gnuplot(char *filename);
		void draw_pov(char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void draw_pov(char *filename);
		void store_cell_volumes(fpoint *bb);
		void print_all(ostream &os);
		void print_all();
		void print_all(char *filename);
		void print_all_neighbor(ostream &os);
		void print_all_neighbor();
		void print_all_neighbor(char *filename);
		template<class n_option>
		inline void compute_cell_slow(voronoicell_base<n_option> &c,int s,int i);
		template<class n_option>
		inline void compute_cell_slow(voronoicell_base<n_option> &c,int s,int i,fpoint x,fpoint y,fpoint z);
		template<class n_option>
		inline void compute_cell(voronoicell_base<n_option> &c,int s,int i);
		template<class n_option>
		inline void compute_cell(voronoicell_base<n_option> &c,int s,int i,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z);
		void guess_length_scale();
	protected:
		/** The amount of memory in the array structure for each
		 * particle. This is set to 3 when the basic class is
		 * initialized, so that the array holds (x,y,z) positions. If
		 * the container class is initialized as part of the derived
		 * class container_poly, then this is set to 4, to also hold
		 * the particle radii. */
		const int sz;
		/** The minimum x coordinate of the container. */
		const fpoint ax;
		/** The maximum x coordinate of the container. */
		const fpoint bx;
		/** The minimum y coordinate of the container. */
		const fpoint ay;
		/** The maximum y coordinate of the container. */
		const fpoint by;
		/** The minimum z coordinate of the container. */
		const fpoint az;
		/** The maximum z coordinate of the container. */
		const fpoint bz;
		/** The inverse box length in the x direction, set to
		 * nx/(bx-ax). */
		const fpoint xsp;
		/** The inverse box length in the y direction, set to
		 * ny/(by-ay). */
		const fpoint ysp;
		/** The inverse box length in the z direction, set to
		 * nz/(bz-az). */
		const fpoint zsp;
		/** The number of boxes in the x direction. */
		const int nx;
		/** The number of boxes in the y direction. */
		const int ny;
		/** The number of boxes in the z direction. */
		const int nz;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines which step through boxes in
		 * sequence. */
		const int nxy;
		/** A constant, set to the value of nx*ny*nz, which is used in
		 * the routines which step through boxes in sequence. */
		const int nxyz;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** A boolean value that determines if the z coordinate in
		 * periodic or not. */
		const bool zperiodic;
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
		/** This contains the maximum radius of a particle in the
		 * container. For the standard container class, where all
		 * particles are assumed to have the same radius, this number
		 * does nothing. However, for the derived container_poly class,
		 * this number is computed each time a particle is added to the
		 * container. */
		fpoint max_radius;
		template<class n_option>
		inline void print_all(ostream &os,voronoicell_base<n_option> &c);
		template<class n_option>
		inline void initialize_voronoicell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
		void add_particle_memory(int i);
	private:
		friend class facets_search;
		friend class facets_loop;
		facets_search search;
};

/** This is a derived version of the container class for handling polydisperse
 * particles via a radical Voronoi tessellation. The main difference is that
 * the data structure is changed so that radii are also stored. New versions
 * of the put and import commands are provided to account for this. All other
 * routines in the container class are unaffected. */
class container_poly : public container {
	public:
		/** This constructor calls the constructor for the base class,
		 * and specifies that sz=4, so that four floating point numbers
		 * are stored for each particle in the container, for the
		 * (x,y,z) positions, and also the particle radii. */
		container_poly(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi) : container(xa,xb,ya,yb,za,zb,xn,yn,zn,xper,yper,zper,memi,4) {};
		void put(int n,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z,fpoint r);
		void import(istream &is);
		inline void import();
		inline void import(char *filename);
};

/** Many of the container routines require scanning over a rectangular sub-grid
 * of blocks, and the routines for handling this are stored in the facets_loop
 * class. A facets_loop class can first be initialized to either calculate the
 * subgrid which is within a distance r of a vector (vx,vy,vz), or a subgrid
 * corresponding to a rectangular box. The routine inc() can then be
 * successively called to step through all the blocks within this subgrid. 
 */
class facets_loop {
	public:
		facets_loop(container *q);
		inline int init(fpoint vx,fpoint vy,fpoint vz,fpoint r,fpoint &px,fpoint &py,fpoint &pz);
		inline int init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax,fpoint &px,fpoint &py,fpoint &pz);
		inline int inc(fpoint &px,fpoint &py,fpoint &pz);
	private:
		int i,j,k,ai,bi,aj,bj,ak,bk,s;
		int ip,jp,kp,aip,ajp,akp,inc1,inc2;
		inline int step_mod(int a,int b);
		inline int step_div(int a,int b);
		inline int step_int(fpoint a);
		fpoint apx,apy,apz;
		const fpoint sx,sy,sz,xsp,ysp,zsp,ax,ay,az;
		const int nx,ny,nz,nxy,nxyz;
		const bool xperiodic,yperiodic,zperiodic;
};

class facets_search {
	public:
		facets_search(container *q);
		inline int init(fpoint vx,fpoint vy,fpoint vz);
		inline int inc(fpoint &px,fpoint &py,fpoint &pz);
		inline void add_list_memory();
	private:
		int ci,cj,ck;
		unsigned int *m;
		int *s;
		unsigned int v;
		int s_start,s_end,s_size;
		const fpoint sx,sy,sz,xsp,ysp,zsp,ax,ay,az;
		const int nx,ny,nz,nxy,nxyz;
		const bool xperiodic,yperiodic,zperiodic;
		const int hx,hy,hz,hxy,hxyz;
};
#endif
