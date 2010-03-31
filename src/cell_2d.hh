#ifndef VOROPP_CELL_2D_HH
#define VOROPP_CELL_2D_HH

#include "config.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

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
		fpoint area();
		fpoint max_radius_squared();
		fpoint perimeter();
		void centroid(fpoint &cx,fpoint &cy);
		inline bool plane(fpoint x,fpoint y,fpoint rs);
	private:
		inline fpoint pos(fpoint x,fpoint y,fpoint rsq,int qp); 
		int *ds;
};

#endif
