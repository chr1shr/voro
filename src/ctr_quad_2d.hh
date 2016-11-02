#ifndef CONTAINER_QUAD_2D_HH
#define CONTAINER_QUAD_2D_HH

#include "config.hh"
#include "common.hh"

namespace voro {

const int qt_max=8;

class quadtree {
	public:
		const double cx;
		const double cy;
		const double lx;
		const double ly;
		const int ps;
		int *id;
		double *p;
		int co;
		quadtree *qsw;
		quadtree *qse;
		quadtree *qnw;
		quadtree *qne;
		quadtree(double cx_,double cy_,double lx_,double ly_);
		~quadtree();
		void put(int i,double x,double y);
		void split();
		void draw_particles(FILE *fp=stdout);
		void draw_cross(FILE *fp=stdout);
		void setup_neighbors();
		void draw_neighbors(FILE *fp=stdout);
		inline void quick_put(int i,double x,double y) {
			id[co]=i;
			p[ps*co]=x;
			p[1+ps*co++]=y;
		}
		inline void add_neighbor(quadtree *qt) {
			if(nco==nmax) add_neighbor_memory();
			nei[nco++]=qt;
		}
	protected:
		quadtree **nei;
		int nco;
		int nmax;
		void we_neighbors(quadtree *qw,quadtree *qe);
		void ns_neighbors(quadtree *qs,quadtree *qn);
		void add_neighbor_memory();
};

class container_quad_2d : public quadtree {
	public:
		using quadtree::draw_particles;
		using quadtree::draw_neighbors;
		/** The minimum x coordinate of the container. */
		const double ax;
		/** The maximum x coordinate of the container. */
		const double bx;
		/** The minimum y coordinate of the container. */
		const double ay;
		/** The maximum y coordinate of the container. */
		const double by;
		container_quad_2d(double ax_,double bx_,double ay_,double by_);
		inline void draw_particles(const char* filename) {
			FILE *fp=safe_fopen(filename,"w");
			draw_particles(fp);
			fclose(fp);
		}
		inline void draw_neighbors(const char* filename) {
			FILE *fp=safe_fopen(filename,"w");
			draw_neighbors(fp);
			fclose(fp);
		}		
		void draw_quadtree(FILE *fp=stdout);
		inline void draw_quadtree(const char* filename) {
			FILE *fp=safe_fopen(filename,"w");
			draw_quadtree(fp);
			fclose(fp);
		}
};

}

#endif
