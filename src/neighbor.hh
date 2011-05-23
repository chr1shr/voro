
#ifndef VOROPP_NEIGHBOR_CC
#define VOROPP_NEIGHBOR_CC

#include "cell.hh"
#include <vector>
using namespace std;

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
		neighbor_none(voronoicell_base<neighbor_none> *vc_) {};
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
		inline void neighbors(vector<int> &v) {v.clear();};
		/** This is a blank placeholder function that does nothing. */
		inline void label_facets() {};
		/** This is a blank placeholder function that does nothing. */
		inline void check_facets() {};
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
		neighbor_track(voronoicell_base<neighbor_track> *vc_);
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
		inline void neighbors(vector<int> &v);
		inline void label_facets();
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
