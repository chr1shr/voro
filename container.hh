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

class loop;

/** The container class represents the whole simulation region. The
 * container constructor sets up the geometry and periodicity, and divides
 * the geometry into rectangular grid of blocks, each of which handles the
 * particles in a particular area. Routines exist for putting in particles,
 * importing particles from standard input, and carrying out Voronoi
 * calculations.
 */
class container {
	public:
		container(f_point xa,f_point xb,f_point ya,f_point yb,f_point za,f_point zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi);
		~container();
		void dump(char *filename);
		void import(istream &is);
		inline void import();
		inline void import(char *filename);
		void regioncount();
		void clear();
		void vdraw_gnuplot(char *filename,f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax);
		inline void vdraw_gnuplot(char *filename);
		void vdraw_pov(char *filename,f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax);
		inline void vdraw_pov(char *filename);
		void vcomputeall(f_point *bb);
		void vprintall(ostream &of);
		inline void vprintall();
		inline void vprintall(char *filename);
		inline void compute_cell(voronoicell &c,int s,int i);
		inline void compute_cell(voronoicell &c,int s,int i,f_point x,f_point y,f_point z);
#ifdef FACETS_RADICAL
		void put(int n,f_point x,f_point y,f_point z,f_point r);
		f_point max_radius;
#else
		void put(int n,f_point x,f_point y,f_point z);
#endif
	private:
		void addparticlemem(int i);
		const f_point ax,bx,ay,by,az,bz;
		const f_point xsp,ysp,zsp;
		const int nx,ny,nz,nxy,nxyz;
		const bool xperiodic,yperiodic,zperiodic;
		int *co,*mem;
		int **id;
		f_point **p;
		friend class loop;
};

/** Many of the container routines require scanning over a rectangular sub-grid
 * of blocks, and the routines for handling this are stored in the loop class.
 * A loop class can first be initialized to either calculate the subgrid which
 * is within a distance r of a vector (vx,vy,vz), or a subgrid corresponding to
 * a rectangular box. The routine inc can then be successively called to step
 * through all the blocks which could be affected.
 */
class loop {
	public:
		loop(container *q);
		inline int init(f_point vx,f_point vy,f_point vz,f_point r,f_point &px,f_point &py,f_point &pz);
		inline int init(f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax,f_point &px,f_point &py,f_point &pz);
		inline int inc(f_point &px,f_point &py,f_point &pz);
	private:
		int i,j,k,ai,bi,aj,bj,ak,bk,s;
		int ip,jp,kp,aip,ajp,akp,inc1,inc2;
		inline int mymod(int a,int b);
		inline int mydiv(int a,int b);
		template <class T>
		inline int myint(T a);
		f_point apx,apy,apz;
		const f_point sx,sy,sz,xsp,ysp,zsp,ax,ay,az;
		const int nx,ny,nz,nxy,nxyz;
		const bool xperiodic,yperiodic,zperiodic;
};
#endif
