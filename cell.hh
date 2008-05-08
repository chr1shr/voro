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
	char *msg;
	fatal_error(char *p) : msg(p) {
		cerr << p << endl;
	}
};

/** Floating point comparisons can be unreliable on some processor
 * architectures, and can produce unpredictable results. On a number of popular
 * Intel processors, floating point numbers are held to higher precision when
 * in registers than when in memory. When a register is swapped from a register
 * to memory, a truncation error, and in some situations this can create
 * circumstances where for two numbers c and d, the program finds c>d first,
 * but later c<d. The programmer has no control over when the swaps between
 * memory and registers occur, and recompiling with slightly different code can
 * give different results. One solution to avoid this is to force the compiler
 * to evaluate everything in memory (e.g. by using the -ffloat-store option in
 * the GNU C++ compiler) but this could be viewed overkill, since it slows the
 * code down, and the extra register precision is useful.
 *
 * In the plane cutting routine of the voronoicell class, we need to reliably
 * know whether a vertex lies inside, outside, or on the cutting plane, since
 * if it changed during the tracing process there would be confusion. This
 * class makes these tests reliable, by storing the results of marginal cases,
 * where the vertex lies within tolerance2 of the cutting plane. If that vertex
 * is tested again, then code looks up the value of the table in a buffer,
 * rather than doing the floating point comparison again. Only vertices which
 * are close to the plane are stored and tested, so this routine should create
 * minimal computational overhead.
 */
class suretest {
	public:
		/** This is a pointer to the array in the voronoicell class
		 * which holds the vertex coordinates.*/
		fpoint *p;
		
		suretest();
		~suretest();
		inline void init(fpoint x,fpoint y,fpoint z,fpoint rsq);
		inline int test(int n,fpoint &ans);
	private:
		/** This stores the current memory allocation for the marginal
		 * cases. */
		int current_marginal;

		/** This stores the total number of marginal points which are
		 * currently in the buffer. */
		int sc;

		int *sn;
		fpoint px;
		fpoint py;
		fpoint pz;
		fpoint prsq;
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

		/** This is an array of pointers to different blocks of memory
		 * for storing */
		int **mep;

		/** This is an array for holding the number */
		int *mec;

		/** */
		int **ed;
		
		/** */
		int *nu;
		
		/** */
		int *ds;
		
		/** This is the auxiliary delete stack, which has size set by
		 * currentdeletesize2.*/
		int *ds2;

		/** This holds the current size of the arrays ed and nu, which
		 * hold the vertex information. If more vertices are created
		 * than can fit in this array, then it is dynamically extended
		 * using the add_memory_vertices routine. */
		int currentvertices;

		/** This holds the current maximum allowed order of a vertex,
		 * which sets the size of the mem, mep, and mec arrays. If a
		 * vertex is created with more vertices than this, the arrays
		 * are dynamically extended using the add_memory_vorder routine.
		 */
		int currentvertexorder;

		/** This sets the size of the main delete stack. */
		int currentdeletesize;

		/** This sets the size of the auxiliary delete stack. */
		int currentdeletesize2;

		/** This in an array with size 3*currentvertices for holding
		 * the positions of the vertices. */ 
		fpoint *pts;

		/** This sets the total number of vertices in the current cell.
		 */
		int p;

		/** This is a class used in the plane routine for carrying out
		 * reliable comparisons of whether points in the cell are
		 * inside, outside, or on the current cutting plane. */
		suretest sure;

		voronoicell();
		~voronoicell();
		void init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void init_octahedron(fpoint l);
		inline void init_test(int n);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c,int d);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c,int d,int e);
		void dump_pov(ostream &os,fpoint x,fpoint y,fpoint z);
		inline void dump_pov(char *filename,fpoint x,fpoint y,fpoint z);
		inline void dump_pov(fpoint x,fpoint y,fpoint z);
		void dump_pov_mesh(ostream &os,fpoint x,fpoint y,fpoint z);		
		inline void dump_pov_mesh(char *filename,fpoint x,fpoint y,fpoint z);
		inline void dump_pov_mesh(fpoint x,fpoint y,fpoint z);
		void dump_gnuplot(ostream &os,fpoint x,fpoint y,fpoint z);
		inline void dump_gnuplot(char *filename,fpoint x,fpoint y,fpoint z);
		inline void dump_gnuplot(fpoint x,fpoint y,fpoint z);
		inline void relcheck();
		inline void duplicate_check();
		inline void construct_relations();
		fpoint volume();
		fpoint maxradsq();
		inline void edgeprint();
		inline bool collapse_order1();
		inline bool collapse_order2();
		inline void perturb(fpoint r);
		void facets(ostream &os);
		inline void facets();
		inline void facets(char *filename);
		void facet_statistics(ostream &os);
		inline void facet_statistics();
		inline void facet_statistics(char *filename);
		virtual void label_facets();
		bool nplane(fpoint x,fpoint y,fpoint z,fpoint rs,int p_id);
		inline bool nplane(fpoint x,fpoint y,fpoint z,int p_id);
		inline bool plane(fpoint x,fpoint y,fpoint z,fpoint rs);
		inline bool plane(fpoint x,fpoint y,fpoint z);
	private:
		/** This holds the number of points currently on the auxiliary delete stack. */
		int stack2;
		void addmemory(int i);
		void add_memory_vertices();
		void add_memory_vorder();
		void add_memory_ds();
		void add_memory_ds2();

		inline int cycle_up(int a,int p);
		inline int cycle_down(int a,int p);
		inline bool delete_connection(int j,int k,bool hand);
		virtual void neighbor_print(ostream &os,int i,int j);

		virtual void neigbors(ostream &os) {};
		virtual void neighbor_main_allocate() {};
		virtual void neighbor_allocate(int i,int m) {};
		virtual void neighbor_main_deallocate() {};
		virtual void neighbor_deallocate(int i) {};
		virtual void neighbor_add_memory_vertices(int i) {};
		virtual void neighbor_add_memory_vorder(int i) {};
		virtual void neighbor_init() {};
		virtual void neighbor_init_octahedron() {};
		virtual void neighbor_set_pointer(int p,int n) {};
		virtual void neighbor_copy(int a,int b,int c,int d) {};
		virtual void neighbor_set(int a,int b,int c) {};
		virtual void neighbor_set_aux1(int k) {};
		virtual void neighbor_copy_aux1(int a,int b) {};
		virtual void neighbor_copy_aux1_shift(int a,int b) {};
		virtual void neighbor_set_aux2(int k) {};
		virtual void neighbor_copy_aux2(int a,int b) {};
		virtual void neighbor_copy_pointer(int a,int b) {};
		virtual void neighbor_set_to_aux1(int j) {};		
		virtual void neighbor_set_to_aux2(int j) {};
		virtual void neighbor_edgeprint(int i) {};
		virtual void neighbor_allocate_aux1(int i) {};
		virtual void neighbor_switch_to_aux1(int i) {};
		virtual void neighbor_copy_to_aux1(int i,int m) {};
		virtual void neighbor_set_to_aux1_offset(int k,int m) {};
};

class voronoicell_neighbor : public voronoicell {
	public :
		int **mne;
		int **ne;
		int *paux1;
		int *paux2;
		void label_facets();
		void facet_check();
		void neighbors(ostream &os);
	private:
		inline void neighbor_main_allocate();
		inline void neighbor_allocate(int i,int m);
		inline void neighbor_main_deallocate();
		inline void neighbor_deallocate(int i);
		inline void neighbor_add_memory_vertices(int i);
		inline void neighbor_add_memory_vorder(int i);
		inline void neighbor_init();
		inline void neighbor_init_octahedron();
		inline void neighbor_print(ostream &os,int i,int j);
		inline void neighbor_set_pointer(int p,int n);
		inline void neighbor_copy(int a,int b,int c,int d);
		inline void neighbor_set(int a,int b,int c);
		inline void neighbor_set_aux1(int k);
		inline void neighbor_copy_aux1(int a,int b);
		inline void neighbor_copy_aux1_shift(int a,int b);
		inline void neighbor_set_aux2(int k);
		inline void neighbor_copy_aux2(int a,int b);
		inline void neighbor_copy_pointer(int a,int b);
		inline void neighbor_set_to_aux1(int j);
		inline void neighbor_set_to_aux2(int j);
		inline void neighbor_edgeprint(int i);
		inline void neighbor_allocate_aux1(int i);
		inline void neighbor_switch_to_aux1(int i);
		inline void neighbor_copy_to_aux1(int i,int m);
		inline void neighbor_set_to_aux1_offset(int k,int m);
};
#endif
