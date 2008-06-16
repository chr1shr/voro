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

/** Structure for printing fatal error messages and exiting */
struct fatal_error {
	/** This routine prints an error message to the standard error.
	 * \param[in] p The message to print. */
	fatal_error(char *p) {cerr << p << endl;}
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
		int check_marginal(int n,fpoint &ans);
		/** This stores the current memory allocation for the marginal
		 * cases. */
		int current_marginal;
		/** This stores the total number of marginal points which are
		 * currently in the buffer. */
		int sc;
		/** This array contains a list of the marginal points, and also
		 * the outcomes of the marginal tests. */
		int *sn;
		/** The x coordinate of the normal vector to the test plane. */
		fpoint px;
		/** The y coordinate of the normal vector to the test plane. */
		fpoint py;
		/** The z coordinate of the normal vector to the test plane. */
		fpoint pz;
		/** The magnitude of the normal vector to the test plane. */
		fpoint prsq;
};

class neighbor_track;

/** This class encapsulates all the routines for storing and calculating a
 * single Voronoi cell. The cell can first be initialized by the init() function
 * to be a rectangular box. The box can then be successively cut by planes
 * using the plane function.  Other routines exist for outputting the cell,
 * computing its volume, or finding the largest distance of a vertex from the
 * cell center.  The cell is described by two arrays. pts[] is a floating point
 * array which holds the vertex positions. ed[] holds the table of edges, and
 * also a relation table that determines how two vertices are connected to one
 * another. The relation table is redundant, but helps speed up the
 * computation. The function relation_check() checks that the relational table
 * is valid. */
template <class n_option>
class voronoicell_base {
	public:
		/** */
		int *mem;
		/** */
		int **mep;
		/** */
		int *mec;
		/** */
		int **ed;
		/** */
		int *nu;
		/** */
		int *ds;
		/** This is the auxiliary delete stack, which has size set by
		 * current_delete2_size.*/
		int *ds2;
		/** This holds the current size of the arrays ed and nu, which
		 * hold the vertex information. If more vertices are created
		 * than can fit in this array, then it is dynamically extended
		 * using the add_memory_vertices routine. */
		int current_vertices;
		/** This holds the current maximum allowed order of a vertex,
		 * which sets the size of the mem, mep, and mec arrays. If a
		 * vertex is created with more vertices than this, the arrays
		 * are dynamically extended using the add_memory_vorder routine.
		 */
		int current_vertex_order;
		/** This sets the size of the main delete stack. */
		int current_delete_size;
		/** This sets the size of the auxiliary delete stack. */
		int current_delete2_size;
		/** This in an array with size 3*current_vertices for holding
		 * the positions of the vertices. */ 
		fpoint *pts;
		/** This sets the total number of vertices in the current cell.
		 */
		int p;
		/** This is a class used in the plane routine for carrying out
		 * reliable comparisons of whether points in the cell are
		 * inside, outside, or on the current cutting plane. */
		suretest sure;
		voronoicell_base();
		virtual ~voronoicell_base();
		void init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void init_octahedron(fpoint l);
		inline void init_tetrahedron(fpoint x0,fpoint y0,fpoint z0,fpoint x1,fpoint y1,fpoint z1,fpoint x2,fpoint y2,fpoint z2,fpoint x3,fpoint y3,fpoint z3);
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
		inline void check_relations();
		inline void check_duplicates();
		inline void construct_relations();
		fpoint volume();
		fpoint maxradsq();
		void print_edges();
		inline void perturb(fpoint r);
		void facets(ostream &os);
		inline void facets();
		inline void facets(char *filename);
		void facet_statistics(ostream &os);
		inline void facet_statistics();
		inline void facet_statistics(char *filename);
		bool nplane(fpoint x,fpoint y,fpoint z,fpoint rs,int p_id);
		inline bool nplane(fpoint x,fpoint y,fpoint z,int p_id);
		inline bool plane(fpoint x,fpoint y,fpoint z,fpoint rs);
		inline bool plane(fpoint x,fpoint y,fpoint z);
		bool plane_intersects(fpoint x,fpoint y,fpoint z,fpoint rs,int &gp);
		bool plane_intersects_guess(fpoint x,fpoint y,fpoint z,fpoint rs,int &gp);
		void label_facets();
		void neighbors(ostream &os);
		void check_facets();
	private:
		/** This holds the number of points currently on the auxiliary delete stack. */
		int stack2;
		/** This object contains all the functions required to carry out the neighbor
		 * computation. If the neighbor_none class is used for n_option, then all these
		 * functions are blank. If the neighbor_track class is used, then the neighbor
		 * track is enabled. All the functions for the n_option classes are declared
		 * inline, so that they should all be completely integrated into the routine
		 * during compilation. */
		n_option neighbor;
		inline int cycle_up(int a,int p);
		inline int cycle_down(int a,int p);
		void add_memory(int i);
		void add_memory_vertices();
		void add_memory_vorder();
		void add_memory_ds();
		void add_memory_ds2();
		inline bool collapse_order1();
		inline bool collapse_order2();
		inline bool delete_connection(int j,int k,bool hand);
		inline bool plane_intersects_track(fpoint x,fpoint y,fpoint z,fpoint rs,int &gp,fpoint g);
		friend class neighbor_track;
};

/** This is a class full of empty routines for neighbor computation. If the
 * voronoicell_base template is instantiated with this class, then it
 * has the effect of switching off all neighbor computation. Since all these
 * routines are declared inline, it should have the effect of a zero speed
 * overhead in the resulting code. */
class neighbor_none {
	public:
		/** This is a blank constructor. */
		neighbor_none(voronoicell_base<neighbor_none> *ivc) {};
		/** This is a blank placeholder function that does nothing. */
		inline void allocate(int i,int m) {};
		/** This is a blank placeholder function that does nothing. */
		inline void add_memory_vertices(int i) {};
		/** This is a blank placeholder function that does nothing. */
		inline void add_memory_vorder(int i) {};
		/** This is a blank placeholder function that does nothing. */
		inline void init() {};
		/** This is a blank placeholder function that does nothing. */
		inline void init_octahedron() {};
		/** This is a blank placeholder function that does nothing. */
		inline void init_tetrahedron() {};
		/** This is a blank placeholder function that does nothing. */
		inline void set_pointer(int p,int n) {};
		/** This is a blank placeholder function that does nothing. */
		inline void copy(int a,int b,int c,int d) {};
		/** This is a blank placeholder function that does nothing. */
		inline void set(int a,int b,int c) {};
		/** This is a blank placeholder function that does nothing. */
		inline void set_aux1(int k) {};
		/** This is a blank placeholder function that does nothing. */
		inline void copy_aux1(int a,int b) {};
		/** This is a blank placeholder function that does nothing. */
		inline void copy_aux1_shift(int a,int b) {};
		/** This is a blank placeholder function that does nothing. */
		inline void set_aux2_copy(int a,int b) {};
		/** This is a blank placeholder function that does nothing. */
		inline void copy_pointer(int a,int b) {};
		/** This is a blank placeholder function that does nothing. */
		inline void set_to_aux1(int j) {};		
		/** This is a blank placeholder function that does nothing. */
		inline void set_to_aux2(int j) {};
		/** This is a blank placeholder function that does nothing. */
		inline void print_edges(int i) {};
		/** This is a blank placeholder function that does nothing. */
		inline void allocate_aux1(int i) {};
		/** This is a blank placeholder function that does nothing. */
		inline void switch_to_aux1(int i) {};
		/** This is a blank placeholder function that does nothing. */
		inline void copy_to_aux1(int i,int m) {};
		/** This is a blank placeholder function that does nothing. */
		inline void set_to_aux1_offset(int k,int m) {};
		/** This is a blank placeholder function that does nothing. */
		inline void print(ostream &os,int i,int j);
		/** This is a blank placeholder function that does nothing. */
		inline void label_facets() {};
		/** This is a blank placeholder function that does nothing. */
		inline void neighbors(ostream &os) {};
		/** This is a blank placeholder function that does nothing. */
		inline void check_facets() {};
};

/** This class encapsulates all the routines which are required to carry out
 * the neighbor tracking. If the voronoicell_base template is instantiated with
 * this class, then the neighbor computation is enabled. All these routines are
 * simple and declared inline, so they should be directly integrated into the
 * functions in the voronoicell class during compilation, without zero function
 * call overhead. */
class neighbor_track {
	public:
		/** */
		int **mne;
		/** */ 
		int **ne;
		neighbor_track(voronoicell_base<neighbor_track> *ivc);
		~neighbor_track();
		/** This is a pointer back to the voronoicell class which created
		 * this class. It is used to reference the members of that
		 * class in computations. */
		voronoicell_base<neighbor_track> *vc;
		inline void allocate(int i,int m);
		inline void add_memory_vertices(int i);
		inline void add_memory_vorder(int i);
		inline void init();
		inline void init_octahedron();
		inline void init_tetrahedron();
		inline void set_pointer(int p,int n);
		inline void copy(int a,int b,int c,int d);
		inline void set(int a,int b,int c);
		inline void set_aux1(int k);
		inline void copy_aux1(int a,int b);
		inline void copy_aux1_shift(int a,int b);
		inline void set_aux2_copy(int a,int b);
		inline void copy_pointer(int a,int b);
		inline void set_to_aux1(int j);
		inline void set_to_aux2(int j);
		inline void print_edges(int i);
		inline void allocate_aux1(int i);
		inline void switch_to_aux1(int i);
		inline void copy_to_aux1(int i,int m);
		inline void set_to_aux1_offset(int k,int m);
		inline void print(ostream &os,int i,int j);
		inline void label_facets();
		inline void neighbors(ostream &os);
		inline void check_facets();
	private:
		/** This is an auxilliary pointer which is used in some of the
		 * low level neighbor operations. */
		int *paux1;
		/** This is a second auxilliary pointer which is used in some
		 * of the low level neighbor operations. */
		int *paux2;
};

/** The basic voronoicell class. */
typedef voronoicell_base<neighbor_none> voronoicell;

/** A neighbor-tracking version of the voronoicell. */
typedef voronoicell_base<neighbor_track> voronoicell_neighbor;
#endif
