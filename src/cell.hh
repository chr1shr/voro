// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file cell.hh
 * \brief Header file for the voronoicell_base template and related classes. */

#ifndef VOROPP_CELL_HH
#define VOROPP_CELL_HH

#include "config.hh"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting. */
void voropp_fatal_error(const char *p,int status) {
	cerr << "voro++: " << p << endl;
	exit(status);
}

/** \brief A class to reliably carry out floating point comparisons, storing
 * marginal cases for future reference.
 *
 * Floating point comparisons can be unreliable on some processor
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

/** \brief A class encapsulating all the routines for storing and calculating
 * a single Voronoi cell.
 *
 * This class encapsulates all the routines for storing and calculating a
 * single Voronoi cell. The cell can first be initialized by the init() function
 * to be a rectangular box. The box can then be successively cut by planes
 * using the plane function. Other routines exist for outputting the cell,
 * computing its volume, or finding the largest distance of a vertex from the
 * cell center.  The cell is described by two arrays. pts[] is a floating point
 * array which holds the vertex positions. ed[] holds the table of edges, and
 * also a relation table that determines how two vertices are connected to one
 * another. The relation table is redundant, but helps speed up the
 * computation. The function check_relations() checks that the relational table
 * is valid. */
template <class n_option>
class voronoicell_base {
	public:
		/** This is a two dimensional array that holds information about
		 * the edge connections of the vertices that make up the cell.
		 * The two dimensional array is not allocated in the usual method.
		 * To account for the fact the different vertices have different
		 * orders, and thus require different amounts of storage, the
		 * elements of ed[i] point to one-dimensional arrays in the mep[]
		 * array of different sizes.
		 *
		 * More specifically, if vertex i has order m, then ed[i]
		 * points to a one-dimensional array in mep[m] that has 2*m+1
		 * entries. The first m elements hold the neighboring edges, so
		 * that the jth edge of vertex i is held in ed[i][j]. The next
		 * m elements hold a table of relations which is redundant but
		 * helps speed up the computation. It satisfies the relation
		 * ed[ed[i][j]][ed[i][m+j]]=i. The final entry holds a back
		 * pointer, so that ed[i+2*m]=i. These are used when
		 * rearranging the memory. */
		int **ed;
		/** This array holds the order of the vertices in the Voronoi cell.
		 * This array is dynamically allocated, with its current size
		 * held by current_vertices. */
		int *nu;
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
		/** This is the index of particular point in the cell, which is used to start
		 * the tracing routines for plane intersection and cutting. These
		 * routines will work starting from any point, but it's often most
		 * efficient to start from the last point considered, since in many cases,
		 * the cell construction algorithm may consider many planes with similar
		 * vectors concurrently. */
		int up;
		/** This is a class used in the plane routine for carrying out
		 * reliable comparisons of whether points in the cell are
		 * inside, outside, or on the current cutting plane. */
		suretest sure;
		voronoicell_base();
		~voronoicell_base();
		void init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void init_octahedron(fpoint l);
		inline void init_tetrahedron(fpoint x0,fpoint y0,fpoint z0,fpoint x1,fpoint y1,fpoint z1,fpoint x2,fpoint y2,fpoint z2,fpoint x3,fpoint y3,fpoint z3);
		inline void init_test(int n);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c,int d);
		inline void add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c,int d,int e);
		void draw_pov(ostream &os,fpoint x,fpoint y,fpoint z);
		inline void draw_pov(const char *filename,fpoint x,fpoint y,fpoint z);
		inline void draw_pov(fpoint x,fpoint y,fpoint z);
		void draw_pov_mesh(ostream &os,fpoint x,fpoint y,fpoint z);
		inline void draw_pov_mesh(const char *filename,fpoint x,fpoint y,fpoint z);
		inline void draw_pov_mesh(fpoint x,fpoint y,fpoint z);
		void draw_gnuplot(ostream &os,fpoint x,fpoint y,fpoint z);
		inline void draw_gnuplot(const char *filename,fpoint x,fpoint y,fpoint z);
		inline void draw_gnuplot(fpoint x,fpoint y,fpoint z);
		inline void check_relations();
		inline void check_duplicates();
		inline void construct_relations();
		fpoint volume();
		fpoint max_radius_squared();
		fpoint total_edge_distance();
		void centroid(fpoint &cx,fpoint &cy,fpoint &cz);
		int number_of_faces();
		void output_vertices(ostream &os);
		void output_vertices(ostream &os,fpoint x,fpoint y,fpoint z);
		void print_edges();
		inline void perturb(fpoint r);
		void facets(ostream &os);
		inline void facets();
		inline void facets(const char *filename);
		void vertex_number_histogram(ostream &os);
		void neighbor_normals(ostream &os);
		bool nplane(fpoint x,fpoint y,fpoint z,fpoint rs,int p_id);
		inline bool nplane(fpoint x,fpoint y,fpoint z,int p_id);
		inline bool plane(fpoint x,fpoint y,fpoint z,fpoint rs);
		inline bool plane(fpoint x,fpoint y,fpoint z);
		bool plane_intersects(fpoint x,fpoint y,fpoint z,fpoint rs);
		bool plane_intersects_guess(fpoint x,fpoint y,fpoint z,fpoint rs);
		void label_facets();
		void neighbors(ostream &os);
		void check_facets();
	private:
		/** This a one dimensional array that holds the current sizes
		 * of the memory allocations for them mep array.*/
		int *mem;
		/** This is a two dimensional array for holding the information
		 * about the edges of the Voronoi cell. mep[p] is a
		 * one-dimensional array for holding the edge information about
		 * all vertices of order p, with each vertex holding 2*p+1
		 * integers of information. The total number of vertices held
		 * on mep[p] is stored in mem[p]. If the space runs out, the
		 * code allocates more using the add_memory() routine. */
		int **mep;
		/** This is a one dimensional array that holds the current
		 * number of vertices of order p that are stored in the mep[p]
		 * array. */
		int *mec;
		/** This is the delete stack, used to store the vertices which
		 * are going to be deleted during the plane cutting procedure.
		 */
		int *ds;
		/** This is the auxiliary delete stack, which has size set by
		 * current_delete2_size. */
		int *ds2;
		/** This holds the number of points currently on the auxiliary
		 * delete stack. */
		int stack2;
		/** This object contains all the functions required to carry
		 * out the neighbor computation. If the neighbor_none class is
		 * used for n_option, then all these functions are blank. If
		 * the neighbor_track class is used, then the neighbor tracking
		 * is enabled. All the functions for the n_option classes are
		 * declared inline, so that they should all be completely
		 * integrated into the routine during compilation. */
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
		inline bool plane_intersects_track(fpoint x,fpoint y,fpoint z,fpoint rs,fpoint g);
		inline void reset_edges();
		inline void neighbor_normals_search(ostream &os,int i,int j,int k);
		friend class neighbor_track;
};

/** \brief A class passed to the voronoicell_base template to switch off
 * neighbor computation.
 *
 * This is a class full of empty routines for neighbor computation. If the
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
		/** This is a blank placeholder function that does nothing. */
		inline void print_neighbor(ostream &os,int i,int j) {};
};

/** \brief A class passed to the voronoicell_base template to switch on the
 * neighbor computation.
 *
 * This class encapsulates all the routines which are required to carry out the
 * neighbor tracking. If the voronoicell_base template is instantiated with
 * this class, then the neighbor computation is enabled. All these routines are
 * simple and declared inline, so they should be directly integrated into the
 * functions in the voronoicell class during compilation, without zero function
 * call overhead. */
class neighbor_track {
	public:
		/** This two dimensional array holds the neighbor information
		 * associated with each vertex. mne[p] is a one dimensional
		 * array which holds all of the neighbor information for
		 * vertices of order p. */
		int **mne;
		/** This is a two dimensional array that holds the neighbor
		 * information associated with each vertex. ne[i] points to a
		 * one-dimensional array in mne[nu[i]]. ne[i][j] holds the
		 * neighbor information associated with the jth edge of vertex
		 * i. It is set to the ID number of the plane that made the
		 * face that is clockwise from the jth edge. */
		int **ne;
		neighbor_track(voronoicell_base<neighbor_track> *ivc);
		~neighbor_track();
		/** This is a pointer back to the voronoicell class which
		 * created this class. It is used to reference the members of
		 * that class in computations. */
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
		inline void print_neighbor(ostream &os,int i,int j);
		inline void check_facets();
	private:
		/** This is an auxiliary pointer which is used in some of the
		 * low level neighbor operations. */
		int *paux1;
		/** This is a second auxiliary pointer which is used in some
		 * of the low level neighbor operations. */
		int *paux2;
};

/** The basic voronoicell class. */
typedef voronoicell_base<neighbor_none> voronoicell;

/** A neighbor-tracking version of the voronoicell. */
typedef voronoicell_base<neighbor_track> voronoicell_neighbor;
#endif
