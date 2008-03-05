// Voronoi calculation header file
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#ifndef FACETS_CELL_HH
#define FACETS_CELL_HH

#include "constants.hh"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// Structure for printing fatal error messages and exiting
struct fatal_error {
    const char *msg;
    fatal_error(const char *p) : msg(p) {
	cerr << p << endl;
    }
};

// Floating point comparisons are notoriously unreliable. This class contains
// routines to carefully test the positions of vertices in the plane cutting
// routine. If a dubious case is encountered, the result of the comparison
// is stored in a table, so it can be accessed later, rather than risking its
// value changing.
class suretest {
	public:
		f_point *p;
		suretest();
		~suretest();
		inline void init(f_point x,f_point y,f_point z,f_point rsq);
		inline int test(int n,f_point &ans);
	private:
		int currentdubious;
		int sc,*sn;
		f_point px,py,pz,prsq;
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
		int *mem,**mep,*mec,**ed,*nu,*ds,*ds2;
		int currentvertices,currentvertexorder;
		int currentdeletesize,currentdeletesize2;
		f_point *pts;
		int p;
		suretest sure;
		voronoicell();
		~voronoicell();
		inline void init(f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax);
		inline void init_octahedron(f_point l);
		inline void init_test();
		inline void add_vertex(f_point x,f_point y,f_point z,int a);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b,int c);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b,int c,int d);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b,int c,int d,int e);
		inline void dumppov(ostream &of,f_point x,f_point y,f_point z);
		void dumppovmesh(ostream &of,f_point x,f_point y,f_point z);
		inline void dumpgnuplot(ostream &of,f_point x,f_point y,f_point z);
		inline void relcheck();
		inline void duplicatecheck();
		inline void relconstruct();
		inline f_point volume();
		inline f_point maxradsq();
		inline void edgeprint();
		inline bool collapseorder1();
		inline bool collapseorder2();
		inline void perturb(f_point r);
		void facets(ostream &of);
		inline void facets();
		inline void facets(char *filename);
		void facet_statistics(ostream &of);
		inline void facet_statistics();
		inline void facet_statistics(char *filename);
#ifdef FACETS_NEIGHBOR
		int **mne,**ne;
		void label_facets();
		bool nplane(f_point x,f_point y,f_point z,f_point rs,int p_id);
		inline bool nplane(f_point x,f_point y,f_point z,int p_id);
		inline bool plane(f_point x,f_point y,f_point z,f_point rs);
#else
		bool plane(f_point x,f_point y,f_point z,f_point rs);
#endif
		inline bool plane(f_point x,f_point y,f_point z);
	private:
		int stack2;
		void addmemory(int i);
		void addmemory_vertices();
		void addmemory_vorder();
		void addmemory_ds();
		void addmemory_ds2();

		inline int vor_up(int a,int p);
		inline int vor_down(int a,int p);
#ifdef FACETS_NEIGHBOR
		inline bool delete_connection(int j,int k,bool hand);
#else
		inline bool delete_connection(int j,int k);
#endif
};
#endif
