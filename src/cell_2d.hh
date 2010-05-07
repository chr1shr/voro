/** \file cell_2d.hh
 * \brief Header file for the voronoicell_2d class. */

#ifndef VOROPP_CELL_2D_HH
#define VOROPP_CELL_2D_HH

#include "config.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void voropp_fatal_error(const char *p,int status) {
	cerr << "voro++: " << p << endl;
	exit(status);
}

/** \brief A class encapsulating all the routines for storing and calculating a
 * single Voronoi cell. */
class voronoicell_2d {
	public:
		/** This holds the current size of the arrays ed and nu, which
		 * hold the vertex information. If more vertices are created
		 * than can fit in this array, then it is dynamically extended
		 * using the add_memory_vertices routine. */
		int current_vertices;
		/** This sets the size of the current delete stack. */
		int current_delete_size;
		/** The total nuber of vertices in the current cell. */
		int p;
		/** This a two dimensional array that holds information about
		 * edge connections between vertices.*/
		int **ed;
		/** This in an array with size 2*current_vertices for holding
		 * the positions of the vertices. */
		fpoint *pts;
		voronoicell_2d();
		~voronoicell_2d();
		void init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax);
		void draw_gnuplot(ostream &os,fpoint x,fpoint y);
		inline void draw_gnuplot(const char *filename,fpoint x,fpoint y);
		inline void draw_gnuplot(fpoint x,fpoint y);
		void draw_pov(ostream &os,fpoint x,fpoint y,fpoint z=0);
		inline void draw_pov(const char *filename,fpoint x,fpoint y,fpoint z=0);
		inline void draw_pov(fpoint x,fpoint y,fpoint z=0);
		inline bool plane(fpoint x,fpoint y,fpoint rs);
		fpoint max_radius_squared();
		fpoint perimeter();
		fpoint area();
		void centroid(fpoint &cx,fpoint &cy);
	private:
		void add_memory_vertices();
		void add_memory_ds();
		inline fpoint pos(fpoint x,fpoint y,fpoint rsq,int qp);
		/** This is the delete stack, used to store the vertices which
		 * are going to be deleted during the plane cutting procedure.
		 */	
		int *ds;
};

#endif
