// Voronoi calculation header file
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#ifndef CONTAINER_HH
#define CONTAINER_HH

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#ifdef OVERFLOW_CHECKING
struct overflow {
    const char *msg;
    overflow(const char *p) : msg(p) {
	cerr << p << endl;
    }
};
#endif

enum out_type{pov,gnuplot};

class loop;

// The container class represents the whole simulation region. The container
// constructor sets up the geometry and periodicity, and divides the geometry
// into rectangular grid of blocks, each of which handles the particles in a
// particular area. Routines exist for putting in particles, importing
// particles from standard input, and carrying out Voronoi calculations.
class container {
	public:
		container(double xa,double xb,double ya,double yb,double za,double zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi);
		void dump(char *filename);
		void put(int n,double x,double y,double z);
		void import();
		void regioncount();
		void clear();
		void vdraw(char *filename,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,out_type ot);
		void vdraw(char *filename,out_type ot);
		void vcomputeall(double *bb);
		void vprintall();
	private:
		const double ax,bx,ay,by,az,bz;
		const double xsp,ysp,zsp;
		const int nx,ny,nz,nxy,nxyz,mem;
		const bool xperiodic,yperiodic,zperiodic;
		int *co;
		int **id;
		double **p;
		friend class loop;
};

// Many of the container routines require scanning over a rectangular sub-grid of
// blocks, and the routines for handling this are stored in the loop class.
// A loop class can first be initialized to either calculate the subgrid which is
// within a distance r of a vector (vx,vy,vz), or a subgrid corresponding to
// a rectangular box. The routine inc can then be successively called to
// step through all the blocks which could be affected.
class loop {
	public:
		loop(container *q);
		inline int init(double vx,double vy,double vz,double r,double &px,double &py,double &pz);
		inline int init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double &px,double &py,double &pz);
		inline int inc(double &px,double &py,double &pz);
	private:
		int i,j,k,ai,bi,aj,bj,ak,bk,s;
		int ip,jp,kp,aip,ajp,akp,inc1,inc2;
		inline int mymod(int a,int b);
		inline int mydiv(int a,int b);
		template <class T>
		inline int myint(T a);
		double apx,apy,apz;
		const double sx,sy,sz,xsp,ysp,zsp,ax,ay,az;
		const int nx,ny,nz,nxy,nxyz;
		const bool xperiodic,yperiodic,zperiodic;
};

// Floating point comparisons are notoriously unreliable. This class contains
// routines to carefully test the positions of vertices in the plane cutting
// routine. If a dubious case is encountered, the result of the comparison
// is stored in a table, so it can be accessed later, rather than risking its
// value changing
class suretest {
	public:
		suretest(double *pts) : p(pts) {};
		inline void init(double x,double y,double z,double rsq);
		inline int test(int n,double &ans);
	private:
		int sc,sn[2*maxdubiouscases];
		double px,py,pz,prsq,*p;
};

// This class encapsulates all the routines for storing and calculating a
// single Voronoi cell. The cell can first be initialized by the init function
// to be a rectangular box. The box can then be successively cut by planes
// using the plane function.  Other routines exist for outputting the cell,
// computing its volume, or finding the largest distance of a vertex from the
// cell center.  The cell is described by three arrays: pts[] which holds the
// vertex positions, ed[] which holds the table of edges, and rl[] which is a
// relational table that determines how two vertices are connected to one
// another. rl[] is redundant, but helps speed up the computation. The function
// relcheck checks that the relational table is valid.
class voronoicell {
	public:
		int mem[maxvertexorder-3];
		int *mep[maxvertexorder-3];
		int mec[maxvertexorder-3];
		int *ed[maxvertices];
		int nu[maxvertices];
		int p;
		double pts[3*maxvertices];
		suretest sure;
		voronoicell();
		inline void init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);
		inline void init_octahedron(double l);
		inline bool plane(double x,double y,double z,double rs);
		inline bool plane(double x,double y,double z);
		inline void dumppov(ofstream &of,double x,double y,double z);
		inline void dumpgnuplot(ofstream &of,double x,double y,double z);
		inline void relcheck();
		inline double volume();
		inline double maxradsq();
		inline void edgeprint(bool extend);
	private:
		void addmemory(int i);
		inline int vor_up(int a,int p);
		inline int vor_down(int a,int p);
};

#endif
