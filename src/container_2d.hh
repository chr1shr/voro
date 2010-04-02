#ifndef VOROPP_CONTAINER_2D_HH
#define VOROPP_CONTAINER_2D_HH

#include "config.hh"
#include "cell.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class container_2d {
	public:
		container_2d(fpoint xa,fpoint xb,fpoint ya,fpoint yb,int xn,int yn,int memi);
		~container_2d();
		void draw_particles(const char *filename);
		void draw_particles();
		void draw_particles(ostream &os);
		void import(istream &is);
		inline void import();
		inline void import(const char *filename);
		void draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax);
		inline void draw_cells_gnuplot(const char *filename);
		template<class n_option>
		inline bool compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s);
		template<class n_option>
		bool compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z);
		void put(int n,fpoint x,fpoint y,fpoint z,fpoint r);
	private:
		/** The minimum x coordinate of the container. */
		const fpoint ax;
		/** The maximum x coordinate of the container. */
		const fpoint bx;
		/** The minimum y coordinate of the container. */
		const fpoint ay;
		/** The maximum y coordinate of the container. */
		const fpoint by;
		/** The inverse box length in the x direction, set to
		 * nx/(bx-ax). */
		const fpoint xsp;
		/** The inverse box length in the y direction, set to
		 * ny/(by-ay). */
		const fpoint ysp;
		/** The number of boxes in the x direction. */
		const int nx;
		/** The number of boxes in the y direction. */
		const int ny;
		/** The number of boxes in the z direction. */
		const int nz;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines which step through boxes in
		 * sequence. */
		const int nxy;
		/** A constant, set to the value of nx*ny*nz, which is used in
		 * the routines which step through boxes in sequence. */
		const int nxyz;
		/** This array holds the number of particles within each
		 * computational box of the container. */
		int *co;
		/** This array holds the maximum amount of particle memory for
		 * each computational box of the container. If the number of
		 * particles in a particular box ever approaches this limit,
		 * more is allocated using the add_particle_memory() function.
		 */
		int *mem;
		/** This array holds the numerical IDs of each particle in each
		 * computational box. */
		int **id;
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		fpoint **p;		
};

#endif
