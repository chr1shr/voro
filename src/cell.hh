// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file cell.hh
 * \brief Header file for the voronoicell_base template and related classes. */

#ifndef VOROPP_CELL_HH
#define VOROPP_CELL_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "config.hh"

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
		/** This sets the total number of vertices in the current cell.
		 */
		int p;
		/** This is the index of particular point in the cell, which is
		 * used to start the tracing routines for plane intersection
		 * and cutting. These routines will work starting from any
		 * point, but it's often most efficient to start from the last
		 * point considered, since in many cases, the cell construction
		 * algorithm may consider many planes with similar vectors
		 * concurrently. */
		int up;
		/** This is a two dimensional array that holds information
		 * about the edge connections of the vertices that make up the
		 * cell. The two dimensional array is not allocated in the
		 * usual method. To account for the fact the different vertices
		 * have different orders, and thus require different amounts of
		 * storage, the elements of ed[i] point to one-dimensional
		 * arrays in the mep[] array of different sizes.
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
		/** This array holds the order of the vertices in the Voronoi
		 * cell. This array is dynamically allocated, with its current
		 * size held by current_vertices. */
		int *nu;
		/** This in an array with size 3*current_vertices for holding
		 * the positions of the vertices. */
		double *pts;
		/** This object contains all the functions required to carry
		 * out the neighbor computation. If the neighbor_none class is
		 * used for n_option, then all these functions are blank. If
		 * the neighbor_track class is used, then the neighbor tracking
		 * is enabled. All the functions for the n_option classes are
		 * declared inline, so that they should all be completely
		 * integrated into the routine during compilation. */
		n_option neighbor;
		voronoicell_base();
		~voronoicell_base();
		void init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax); //XXX
		inline void init_octahedron(double l); //XXX
		inline void init_tetrahedron(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3); //XXX
		void translate(double x,double y,double z);
		void draw_pov(double x,double y,double z,FILE *fp=stdout);
		inline void draw_pov(double x,double y,double z,const char *filename) {
			FILE *fp(fopen(filename,"w"));
			if(fp==NULL) voropp_fatal_error("Unable to open file",VOROPP_FILE_ERROR);
			draw_pov(x,y,z,fp);
			fclose(fp);
		};
		void draw_pov_mesh(double x,double y,double z,FILE *fp=stdout);
		inline void draw_pov_mesh(double x,double y,double z,const char *filename) {
			FILE *fp(fopen(filename,"w"));
			if(fp==NULL) voropp_fatal_error("Unable to open file",VOROPP_FILE_ERROR);
			draw_pov_mesh(x,y,z,fp);
			fclose(fp);
		}
		void draw_gnuplot(double x,double y,double z,FILE *fp=stdout);
		inline void draw_gnuplot(double x,double y,double z,const char *filename) {
			FILE *fp(fopen(filename,"w"));
			if(fp==NULL) voropp_fatal_error("Unable to open file",VOROPP_FILE_ERROR);
			draw_gnuplot(x,y,z,fp);
			fclose(fp);
		}
		double volume();
		double max_radius_squared();
		double total_edge_distance();
		double surface_area();
		void centroid(double &cx,double &cy,double &cz);
		int number_of_faces();
		int number_of_edges();
		void vertex_orders(vector<int> &v);
		void vertices(vector<double> &v);
		void vertices(vector<double> &v,double x,double y,double z);
		void face_areas(vector<double> &v);
		void face_orders(vector<int> &v);
		void face_freq_table(vector<int> &v);
		void face_vertices(vector<int> &v);
		void face_perimeters(vector<double> &v);
		void normals(vector<double> &v);
		void neighbors(vector<int> &v);
		inline void output_custom(const char *format,FILE *fp=stdout) {output_custom(format,0,0,0,0,default_radius,fp);}
		void output_custom(const char *format,int i,double x,double y,double z,double r,FILE *fp=stdout);
		bool nplane(double x,double y,double z,double rs,int p_id); //XXX
		/** This routine calculates the modulus squared of the vector
		 * before passing it to the main nplane() routine with full
		 * arguments.
		 * \param[in] (x,y,z) the vector to cut the cell by.
		 * \param[in] p_id the plane ID (for neighbor tracking only).
		 * \return False if the plane cut deleted the cell entirely,
		 *         true otherwise. */
		inline bool nplane(double x,double y,double z,int p_id) {
			double rsq=x*x+y*y+z*z;
			return nplane(x,y,z,rsq,p_id);
		}
		/** This version of the plane routine just makes up the plane
		 * ID to be zero. It will only be referenced if neighbor
		 * tracking is enabled.
		 * \param[in] (x,y,z) the vector to cut the cell by.
		 * \param[in] rsq the modulus squared of the vector.
		 * \return False if the plane cut deleted the cell entirely,
		 *         true otherwise. */
		inline bool plane(double x,double y,double z,double rsq) {
			return nplane(x,y,z,rsq,0);
		}
		/** Cuts a Voronoi cell using the influence of a particle at
		 * (x,y,z), first calculating the modulus squared of this
		 * vector before passing it to the main nplane() routine. Zero
		 * is supplied as the plane ID, which will be ignored unless
		 * neighbor tracking is enabled.
		 * \param[in] (x,y,z) the vector to cut the cell by.
		 * \return False if the plane cut deleted the cell entirely,
		 *         true otherwise. */
		inline bool plane(double x,double y,double z) {
			double rsq=x*x+y*y+z*z;
			return nplane(x,y,z,rsq,0);
		}
		bool plane_intersects(double x,double y,double z,double rs);
		bool plane_intersects_guess(double x,double y,double z,double rs);
		void construct_relations();
		void check_relations();
		void check_duplicates();
		void print_edges();
		void check_facets();
		/** This is a simple inline function for picking out the index
		 * of the next edge counterclockwise at the current vertex.
		 * \param[in] a the index of an edge of the current vertex.
		 * \param[in] p the number of the vertex.
		 * \return 0 if a=nu[p]-1, or a+1 otherwise. */
		inline int cycle_up(int a,int p) {return a==nu[p]-1?0:a+1;}
		/** This is a simple inline function for picking out the index
		 * of the next edge clockwise from the current vertex.
		 * \param[in] a the index of an edge of the current vertex.
		 * \param[in] p the number of the vertex.
		 * \return nu[p]-1 if a=0, or a-1 otherwise. */
		inline int cycle_down(int a,int p) {return a==0?nu[p]-1:a-1;}
	private:
		/** This a one dimensional array that holds the current sizes
		 * of the memory allocations for them mep array.*/
		int *mem;
		/** This is a one dimensional array that holds the current
		 * number of vertices of order p that are stored in the mep[p]
		 * array. */
		int *mec;
		/** This is a two dimensional array for holding the information
		 * about the edges of the Voronoi cell. mep[p] is a
		 * one-dimensional array for holding the edge information about
		 * all vertices of order p, with each vertex holding 2*p+1
		 * integers of information. The total number of vertices held
		 * on mep[p] is stored in mem[p]. If the space runs out, the
		 * code allocates more using the add_memory() routine. */
		int **mep;
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
		/** This stores the current memory allocation for the marginal
		 * cases. */
		int current_marginal;
		/** This stores the total number of marginal points which are
		 * currently in the buffer. */
		int n_marg;
		/** This array contains a list of the marginal points, and also
		 * the outcomes of the marginal tests. */
		int *marg;
		/** The x coordinate of the normal vector to the test plane. */
		double px;
		/** The y coordinate of the normal vector to the test plane. */
		double py;
		/** The z coordinate of the normal vector to the test plane. */
		double pz;
		/** The magnitude of the normal vector to the test plane. */
		double prsq;
		void add_memory(int i); //XXX
		void add_memory_vertices(); //XXX
		void add_memory_vorder(); //XXX
		void add_memory_ds();
		void add_memory_ds2();
		inline bool collapse_order1(); //XXX
		inline bool collapse_order2(); //XXX
		inline bool delete_connection(int j,int k,bool hand); //XXX
		inline bool plane_intersects_track(double x,double y,double z,double rs,double g);
		inline void reset_edges();
		inline void normals_search(vector<double> &v,int i,int j,int k);
		inline bool search_edge(int l,int &m,int &k);
		inline int m_test(int n,double &ans);
		int check_marginal(int n,double &ans);
		friend class neighbor_track;
};

#endif
