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

class voronoicell_2d {
	public:
		int current_vertices;
		int current_delete_size;
		int p;
		int up;
		int **ed;
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
		inline fpoint pos(fpoint x,fpoint y,fpoint rsq,int qp); 
		int *ds;
};

#endif
