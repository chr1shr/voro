/** \file container_2d.hh
 * \brief Header file for the container_2d class. */

#ifndef VOROPP_CONTAINER_2D_HH
#define VOROPP_CONTAINER_2D_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

#include "config.hh"
#include "cell_2d.hh"

class voropp_loop_2d;

/** \brief A class representing the whole 2D simulation region.
 *
 * The container class represents the whole simulation region. The
 * container constructor sets up the geometry and periodicity, and divides
 * the geometry into rectangular grid of blocks, each of which handles the
 * particles in a particular area. Routines exist for putting in particles,
 * importing particles from standard input, and carrying out Voronoi
 * calculations. */
class container_2d {
	public:
		/** degbugging aid. */
		bool debug;
		/** The minimum x coordinate of the container. */
		const double ax;
		/** The maximum x coordinate of the container. */
		const double bx;
		/** The minimum y coordinate of the container. */
		const double ay;
		/** The maximum y coordinate of the container. */
		const double by;
		/** The box length in the x direction, set to (bx-ax)/nx. */
		const double boxx;
		/** The box length in the y direction, set to (by-ay)/ny. */
		const double boxy;
		/** The inverse box length in the x direction. */
		const double xsp;
		/** The inverse box length in the y direction. */
		const double ysp;
		/** The number of boxes in the x direction. */
		const int nx;
		/** The number of boxes in the y direction. */
		const int ny;
		/** A constant, set to the value of nx multiplied by ny, which
		 * is used in the routines which step through boxes in
		 * sequence. */
		const int nxy;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** A boolean value that specifies whether the domain will be convex or non-convex.
		* affected methods include cell_2d:plane, container_2d:import, container_2d:initialize_vo
		ronoicell
		* and container_2d:compute_cellsphere. */
		const bool convex;
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
		/** For each computational box, this holds the boundary walls that pass through that box.
		*/
		int **wid;  
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		double **p;
		/** A two dimensional array for holding the number of labels per
		 * particle, plus other status information. */
		unsigned int **nlab;
		/** A two dimensional array of pointers to label information.
		 */
		int ***plab;
		/** An array created from *tmp and is referenced my *soip. soi
		 * stands for "spheres of influence" */
		int *soi;
		/** The current number of boundary points. */
		int noofbnds;
		/** Specifies the memory allocation for the boundary points. */
		int bnds_size;
		/** This array holds the vertices of the non-convex domain. If
		 * the domain is convex this variable is not initialized. The
		 * format is (bnds[0],bnds[1])=(x1,y1) for some (x1,y1) 
		 * boundary point and (bnds[2],bnds[3])=(x2,y2) is the next
		 * vertex in the boundary from (x1,y1) if we are going
		 * clockwise. */
		double *bnds;
		/** Information about how the boundary vertices are connected.
		 * A single integer is stored for each vertex, giving the index
		 * of the connecting vertex in the counter-clockwise sense. */
		int *edb;
		/** Contains information about which points are boundary points,
		i.e. lie at a nonconvexity. **/
		int **bndpts;
		/** Contains information about which boundary points are problem points **/
		bool *probpts;
		container_2d(double xa,double xb,double ya,double yb,int xn,int yn,bool xper,bool yper,bool convex_,int memi);
		~container_2d();
		void setup();
		int crossproductz(double x1, double y1, double x2, double y2);
		void semi_circle_labelling(double x1, double y1, double x2, double y2, int wid);

		void import(FILE *fp=stdin);
		/** Imports a list of particles from a file.
		 * \param[in] filename the file to read from. */
		inline void import(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(fp);
			fclose(fp);
		}
		void tag_walls(double x1, double y1, double x2, double y2, int wid);
		inline void tag(int wallid, int box){
			int length=wid[box][0], *a;
			if(length%6==0){
				a=new int[length+6];
				for(int j=0; j<length; j++){
					a[j]=wid[box][j];
				}
				delete [] wid[box]; wid[box]=a;
			}
			wid[box][length]=wallid;
			wid[box][0]++;
		}
		void draw_particles(FILE *fp=stdout);
		/** Dumps all the particle positions and IDs to a file.
		 * \param[in] filename the file to write to. */
		inline void draw_particles(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_particles(fp);
			fclose(fp);
		}
		void draw_particles_pov(FILE *fp=stdout);
		/** Dumps all the particles positions in POV-Ray format.
		 * \param[in] filename the file to write to. */
		inline void draw_particles_pov(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_particles_pov(fp);
			fclose(fp);
		}
		void draw_cells_gnuplot(FILE *fp=stdout);
		/** Computes the Voronoi cells for all particles and saves the
		 * output in gnuplot format.
		 * \param[in] filename the file to write to. */
		inline void draw_cells_gnuplot(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_cells_gnuplot(fp);
			fclose(fp);
		}
		void draw_cells_pov(FILE *fp=stdout);
		/** Computes the Voronoi cells for all particles and saves the
		 * output in POV-Ray format.
		 * \param[in] filename the file to write to. */
		inline void draw_cells_pov(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_cells_pov(fp);
			fclose(fp);
		}
		void draw_boundary(FILE *fp=stdout);
		inline void draw_boundary(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_boundary(fp);
			fclose(fp);
		}
		void print_custom(const char *format,FILE *fp=stdout);
		/** Computes the Voronoi cells for all particles in the
		 * container, and for each cell, outputs a line containing
		 * custom information about the cell structure. The output
		 * format is specified using an input string with control
		 * sequences similar to the standard C printf() routine.
		 * \param[in] format the format of the output lines, using
		 *                   control sequences to denote the different
		 *                   cell statistics.
		 * \param[in] filename the file to write to. */
		inline void print_custom(const char *format,const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			print_custom(format,fp);
			fclose(fp);
		}
		double sum_cell_areas();
		void compute_all_cells();
		/** An overloaded version of the compute_cell_sphere routine,
		 * that sets up the x and y variables.
		 *\param[in,out] c a reference to a voronoicell object.
		 * \param[in] (i,j) the coordinates of the block that the test
		 *                  particle is in.
		 * \param[in] ij the index of the block that the test particle
		 *               is in, set to i+nx*j.
		 * \param[in] s the index of the particle within the test
		 *              block.
		 * \return False if the Voronoi cell was completely removed
		 * during the computation and has zero volume, true otherwise.
		 */
		inline bool compute_cell_sphere(voronoicell_2d &c,int i,int j,int ij,int s) {
			double x=p[s][2*ij],y=p[s][2*ij+1];
			return compute_cell_sphere(c,i,j,ij,s,x,y);
		}
		bool compute_cell_sphere(voronoicell_2d	&c,int i,int j,int ij,int s,double x,double y);
		bool initialize_voronoicell(voronoicell_2d &c,double x,double y);
		bool initialize_voronoicell_nonconvex(voronoicell_2d &c,double x, double y, int bid);
		bool OKCuttingParticle(double gx, double gy, int gbox, int gindex, double cx, double cy,
					int cbox, int cindex);
		void put(int n,double x,double y, int bndloc);
		void clear();
	private:
		/** A temporary array used in label creation. */
		int *tmp;
		/** The next free position in the temporary label array. */
		int *tmpp;
		/** The end of the temporary label array. */
		int *tmpe;
		inline bool put_locate_block(int &ij,double &x,double &y);
		inline bool put_remap(int &ij,double &x,double &y);
		inline double max(double a, double b) {return b<a?a:b;}
		inline double min(double a, double b) {return (b<a)?b:a;}
		inline double dist_squared(double x1,double y1,double x2,double y2) {
			return (((y2-y1)*(y2-y1))+((x2-x1)*(x2-x1)));
		}
		void create_label_table();
		void add_temporary_label_memory();
		void add_particle_memory(int i);
		void add_boundary_memory();
		/** Custom int function, that gives consistent stepping for
		 * negative numbers. With normal int, we have
		 * (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1). With this routine, we
		 * have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1). */
		inline int step_int(double a) {return a<0?int(a)-1:int(a);}
		/** Custom modulo function, that gives consistent stepping for
		 * negative numbers. */
		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
		/** Custom integer division function, that gives consistent
		 * stepping for negative numbers. */
		inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
		friend class voropp_loop_2d;
};


/** \brief A class to handle loops on regions of the container handling
 * non-periodic and periodic boundary conditions.
 *
 * Many of the container routines require scanning over a rectangular sub-grid
 * of blocks, and the routines for handling this are stored in the
 * voropp_loop_2d class. A voropp_loop_2d class can first be initialized to
 * either calculate the subgrid which is within a distance r of a vector
 * (vx,vy,vz), or a subgrid corresponding to a rectangular box. The routine
 * inc() can then be successively called to step through all the blocks within
 * this subgrid.
 */
class voropp_loop_2d {
	public:
		voropp_loop_2d(container_2d &con);
		int init(double vx,double vy,double r,double &px,double &py);
		int init(double xmin,double xmax,double ymin,double ymax,double &px,double &py);
		int inc(double &px,double &py);
		/** The current block index in the x direction, referencing a
		 * real cell in the range 0 to nx-1. */
		int ip;
		/** The current block index in the y direction, referencing a
		 * real cell in the range 0 to ny-1. */
		int jp;
	private:
		const double boxx,boxy,xsp,ysp,ax,ay;
		const int nx,ny,nxy;
		const bool xperiodic,yperiodic;
		double apx,apy;
		int i,j,ai,bi,aj,bj,s;
		int aip,ajp,inc1;
		/** Custom modulo function, that gives consistent stepping for
		 * negative numbers. */
		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
		/** Custom integer division function, that gives consistent
		 * stepping for negative numbers. */
		inline int step_div(int a,int b) {return a>=0?a/b:-1+(a+1)/b;}
		/** Custom int function, that gives consistent stepping for
		 * negative numbers. With normal int, we have
		 * (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1). With this routine, we
		 * have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1). */
		inline int step_int(double a) {return a<0?int(a)-1:int(a);}
};

#endif
