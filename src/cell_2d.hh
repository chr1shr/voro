/** \file cell_2d.hh
 * \brief Header file for the voronoicell_2d class. */

#ifndef VOROPP_CELL_2D_HH
#define VOROPP_CELL_2D_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "common.hh"
#include "config.hh"

/** \brief A class encapsulating all the routines for storing and calculating a
 * single Voronoi cell. */
class voronoicell_2d {
	public:
		/** This holds the current size of the ed and pts arrays. If
		 * more vertices are created than can fit in these arrays, then
		 * they are dynamically extended using the add_memory_vertices
		 * routine. */
		int current_vertices;
		/** This sets the size of the current delete stack. */
		int current_delete_size;
		/** The total nuber of vertices in the current cell. */
		int p;
		/** An array with size 2*current_vertices holding information
		 * about edge connections between vertices.*/
		int *ed;
		/** An array with size 2*current_vertices for holding
		 * the positions of the vertices. */
		double *pts;
		voronoicell_2d();
		~voronoicell_2d();
		void init(double xmin,double xmax,double ymin,double ymax);
		void draw_gnuplot(double x,double y,FILE *fp=stdout);
		/** Outputs the edges of the Voronoi cell in gnuplot format to
		 * an output stream.
		 * \param[in] (x,y) a displacement vector to be added to the
		 *            cell's position.
		 * \param[in] filename the file to write to. */
		inline void draw_gnuplot(double x,double y,const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_gnuplot(x,y,fp);
			fclose(fp);
		}
		void draw_pov(double x,double y,FILE *fp=stdout);
		/** Outputs the edges of the Voronoi cell in POV-Ray format to
		 * an open file stream, displacing the cell by given vector.
		 * \param[in] (x,y,z) a displacement vector to be added to the
		 *                    cell's position.
		 * \param[in] filename the file to write to. */
		inline void draw_pov(double x,double y,const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_pov(x,y,z,fp);
			fclose(fp);
		}
		void output_custom(const char *format,int i,double x,double y,double r,FILE *fp=stdout);
		/** Computes the Voronoi cells for all particles in the
		 * container, and for each cell, outputs a line containing
		 * custom information about the cell structure. The output
		 * format is specified using an input string with control
		 * sequences similar to the standard C printf() routine.
		 * \param[in] format the format of the output lines, using
		 *                   control sequences to denote the different
		 *                   cell statistics.
		 * \param[in] i the ID of the particle associated with this
		 *              Voronoi cell.
		 * \param[in] (x,y) the position of the particle associated
		 *                  with this Voronoi cell.
		 * \param[in] r a radius associated with the particle.
		 * \param[in] filename the file to write to. */
		inline void output_custom(const char *format,int i,double x,double y,double r,const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			output_custom(format,i,x,y,r,fp);
			fclose(fp);
		}
		bool plane(double x,double y,double rs);
		double max_radius_squared();
		double perimeter();
		double area();
		void centroid(double &cx,double &cy);
	private:
		void add_memory_vertices();
		void add_memory_ds(int *&stackp);
		/** Computes the distance of a Voronoi cell vertex to a plane.
		 * \param[in] (x,y) the normal vector to the plane.
		 * \param[in] rsq the distance along this vector of the plane.
		 * \param[in] qp the index of the vertex to consider. */
		inline double pos(double x,double y,double rsq,int qp) {
			return x*pts[2*qp]+y*pts[2*qp+1]-rsq;
		}
		/** The delete stack, used to store the vertices that are
		 * deleted during the plane cutting procedure. */
		int *ds;
		/** A pointer to the end of the delete stack, used to detect
		 * when it is full. */
		int *stacke;
};

#endif
