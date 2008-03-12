// Voronoi calculation header file
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#ifndef FACETS_CELL_HH
#define FACETS_CELL_HH

#include "config.hh"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

/** Structure for printing fatal error messages and exiting
 */
struct fatal_error {
	const char *msg;
	fatal_error(const char *p) : msg(p) {
		cerr << p << endl;
	}
};

/** Floating point comparisons can be unreliable on some processor architectures,
 * and can produce unpredictable results. On a number of popular Intel processors,
 * floating point numbers are held to higher precision when in registers
 * than when in memory. When a register is swapped from a register to memory,
 * a truncation error, and in some situations this can create circumstances
 * where for two numbers c and d, the program finds c>d first, but later c<d.
 * The programmer has no control over when the swaps between memory and registers
 * occur, and recompiling with slightly different code can give different results.
 * One solution to avoid this is to force the compiler to evaluate everything in
 * memory (e.g. by using the -ffloat-store option in the GNU C++ compiler) but this
 * could be viewed overkill, since it slows the code down, and the extra
 * register precision is useful.
 *
 * In the plane cutting routine of the voronoicell class, we need to reliably know whether
 * a vertex lies inside, outside, or on the cutting plane, since if it changed during
 * the tracing process there would be confusion. This class makes these tests reliable,
 * by storing the results of marginal cases, where the vertex lies within tolerance2
 * of the cutting plane. If that vertex is tested again, then code looks up the value
 * of the table in a buffer, rather than doing the floating point comparison again.
 * Only vertices which are close to the plane are stored and tested, so this routine
 * should create minimal computational overhead.
 */
class suretest {
	public:
		/** This is a pointer to the array in the voronoicell class
		 * which holds the vertex coordinates.*/
		f_point *p;
		
		suretest();
		~suretest();
		inline void init(f_point x,f_point y,f_point z,f_point rsq);
		inline int test(int n,f_point &ans);
	private:
		/** This stores the current memory allocation for the marginal
		 * cases. */
		int currentdubious;

		/** This stores the total number of marginal points which are
		 * currently in the buffer. */
		int sc;

		int *sn;
		f_point px;
		f_point py;
		f_point pz;
		f_point prsq;
};

/** This class encapsulates all the routines for storing and calculating a
 * single Voronoi cell. The cell can first be initialized by the init function
 * to be a rectangular box. The box can then be successively cut by planes
 * using the plane function.  Other routines exist for outputting the cell,
 * computing its volume, or finding the largest distance of a vertex from the
 * cell center.  The cell is described by three arrays: pts[] which holds the
 * vertex positions, ed[] which holds the table of edges, and rl[] which is a
 * relational table that determines how two vertices are connected to one
 * another. rl[] is redundant, but helps speed up the computation. The function
 * relcheck checks that the relational table is valid.
 */
class voronoicell {
	public:
		/** This is an array for holding the */
		int *mem;

		/** This is an array of pointers to different blocks of memory for storing */
		int **mep;

		/** This is an array for holding the number */
		int *mec;

		/** */
		int **ed;
		
		/** */
		int *nu;
		
		/** */
		int *ds;
		
		/** This is the auxiliary delete stack, which has size set by currentdeletesize2 .*/
		int *ds2;

		/** This holds the current size of the arrays ed and nu, which hold the vertex
		 * information. If more vertices are created than can fit in this array, then it
		 * is dynamically extended using the addmemory_vertices routine. */
		int currentvertices;

		/** This holds the current maximum allowed order of a vertex, which sets the size
		 * of the mem, mep, and mec arrays. If a vertex is created with more vertices
		 * than this, the arrays are dynamically extended using the addmemory_vorder
		 * routine. */
		int currentvertexorder;

		/** This sets the size of the main delete stack. */
		int currentdeletesize;

		/** This sets the size of the auxiliary delete stack. */
		int currentdeletesize2;

		/** This in an array with size 3*currentvertices for holding the positions of the
		 * vertices. */ 
		f_point *pts;

		/* This sets the total number of vertices in the current cell.  */
		int p;

		/* This is a class used in the plane routine for carrying out reliable comparisons
		 * of whether points in the cell are inside, outside, or on the current
		 * cutting plane. */
		suretest sure;

		voronoicell();
		~voronoicell();
		inline void init(f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax);
		inline void init_octahedron(f_point l);
		inline void init_test(int n);
		inline void add_vertex(f_point x,f_point y,f_point z,int a);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b,int c);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b,int c,int d);
		inline void add_vertex(f_point x,f_point y,f_point z,int a,int b,int c,int d,int e);
		inline void dumppov(ostream &of,f_point x,f_point y,f_point z);
		inline void dumppov(char *filename,f_point x,f_point y,f_point z);
		inline void dumppov(f_point x,f_point y,f_point z);
		inline void dumppovmesh(ostream &of,f_point x,f_point y,f_point z);		
		inline void dumppovmesh(char *filename,f_point x,f_point y,f_point z);
		inline void dumppovmesh(f_point x,f_point y,f_point z);
		inline void dumpgnuplot(ostream &of,f_point x,f_point y,f_point z);
		inline void dumpgnuplot(char *filename,f_point x,f_point y,f_point z);
		inline void dumpgnuplot(f_point x,f_point y,f_point z);
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
		/** <em>Only present for the neighbor-tracking version.</em> */
		int **mne;
		int **ne;
		void label_facets();
		void facet_check();
		void neighbors(ostream &of);
		bool nplane(f_point x,f_point y,f_point z,f_point rs,int p_id);
		inline bool nplane(f_point x,f_point y,f_point z,int p_id);
		inline bool plane(f_point x,f_point y,f_point z,f_point rs);
#else
		bool plane(f_point x,f_point y,f_point z,f_point rs);
#endif
		inline bool plane(f_point x,f_point y,f_point z);
	private:
		/** This holds the number of points currently on the auxiliary delete stack. */
		int stack2;
		void addmemory(int i);
		void addmemory_vertices();
		void addmemory_vorder();
		void addmemory_ds();
		void addmemory_ds2();

		inline int cycle_up(int a,int p);
		inline int cycle_down(int a,int p);
#ifdef FACETS_NEIGHBOR
		inline bool delete_connection(int j,int k,bool hand);
#else
		inline bool delete_connection(int j,int k);
#endif
};
#endif
