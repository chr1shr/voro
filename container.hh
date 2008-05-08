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

/** The container class represents the whole simulation region. The
 * container constructor sets up the geometry and periodicity, and divides
 * the geometry into rectangular grid of blocks, each of which handles the
 * particles in a particular area. Routines exist for putting in particles,
 * importing particles from standard input, and carrying out Voronoi
 * calculations.
 */
class container {
	public:
		container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi,int isz);
		container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi);
		virtual ~container();
		void dump(char *filename);
		virtual void import(istream &is);
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
		inline void compute_cell(voronoicell &c,int s,int i);
		virtual void compute_cell(voronoicell &c,int s,int i,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z);
	protected:
		inline void print_all(ostream &os,voronoicell &c);
		virtual void poly_clear_radius() {};
		inline void initialize_voronoicell(voronoicell &c,fpoint x,fpoint y,fpoint z);
		void add_particle_memory(int i);
		/** The amount of memory in the array structure for each particle. This
		 * is set to 3 when the basic class is initialized, so that the array holds
		 * (x,y,z) positions. If the container class is initialized as part
		 * of the derived class container_poly, then this is set to 4, to also
		 * hold the particle radii. */
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
		/** The inverse box length in the x direction, set to nx/(bx-ax). */
		const fpoint xsp;
		/** The inverse box length in the y direction, set to ny/(by-ay). */
		const fpoint ysp;
		/** The inverse box length in the z direction, set to nz/(bz-az). */
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

		int *co;
		
		int *mem;

		int **id;
		/** A two dimensional array holding particle positions. The first
		 * index labels the computational box. */
		fpoint **p;
		friend class facets_loop;
};

class container_poly : public container {
	public:
		container_poly(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi) : container(xa,xb,ya,yb,za,zb,xn,yn,zn,xper,yper,zper,memi,3), max_radius(0) {};
		void put(int n,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z,fpoint r);
		void import(istream &is);
		void clear();
		inline void compute_cell(voronoicell &c,int s,int i,fpoint x,fpoint y,fpoint z);		
	private:
		inline void poly_clear_radius();
		fpoint max_radius;
};

/** Many of the container routines require scanning over a rectangular sub-grid
 * of blocks, and the routines for handling this are stored in the facets_loop class.
 * A facets_loop class can first be initialized to either calculate the subgrid which
 * is within a distance r of a vector (vx,vy,vz), or a subgrid corresponding to
 * a rectangular box. The routine inc() can then be successively called to step
 * through all the blocks within this subgrid. 
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
#endif
