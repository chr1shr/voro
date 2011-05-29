// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file container.hh
 * \brief Header file for the container_base template and related classes. */

#ifndef VOROPP_CONTAINER_HH
#define VOROPP_CONTAINER_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "config.hh"
#include "common.hh"
#include "v_base.hh"
#include "cell.hh"
#include "v_loops.hh"
#include "v_compute.hh"

/** \brief Pure virtual class from which wall objects are derived.
 *
 * This is a pure virtual class for a generic wall object. A wall object
 * can be specified by deriving a new class from this and specifying the
 * functions.*/
class wall {
	public:
		virtual ~wall() {};
		/** A pure virtual function for testing whether a point is
		 * inside the wall object. */
		virtual bool point_inside(double x,double y,double z) = 0;
		/** A pure virtual function for cutting a cell without
		 * neighbor-tracking with a wall. */
		virtual bool cut_cell(voronoicell &c,double x,double y,double z) = 0;
		/** A pure virtual function for cutting a cell with
		 * neighbor-tracking enabled with a wall. */
		virtual bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) = 0;
};

/** \brief A class for storing a list of pointers to walls. */
class wall_list {
	public:
		/** An array holding pointers to wall objects. */
		wall **walls;
		/** A pointer to the next free position to add a wall pointer.
		 */
		wall **wep;
		wall_list();
		~wall_list();
		/** Adds a wall to the list.
		 * \param[in] w the wall to add. */
		inline void add_wall(wall *w) {
			if(wep==wel) increase_wall_memory();
			*(wep++)=w;
		}
		/** Adds a wall to the list.
		 * \param[in] w a reference to the wall to add. */
		inline void add_wall(wall &w) {add_wall(&w);}
		void add_wall(wall_list &wl);
		/** Determines whether a given position is inside all of the
		 * walls on the list.
		 * \param[in] (x,y,z) the position to test.
		 * \return True if it is inside, false if it is outside. */
		inline bool point_inside_walls(double x,double y,double z) {
			for(wall **wp=walls;wp<wep;wp++) if(!((*wp)->point_inside(x,y,z))) return false;
			return true;
		}
		/** Cuts a Voronoi cell by all of the walls currently on
		 * the list.
		 * \param[in] c a reference to the Voronoi cell class.
		 * \param[in] (x,y,z) the position of the cell.
		 * \return True if the cell still exists, false if the cell is
		 * deleted. */
		template<class c_class>
		bool apply_walls(c_class &c,double x,double y,double z) {
			for(wall **wp=walls;wp<wep;wp++) if(!((*wp)->cut_cell(c,x,y,z))) return false;
			return true;
		}
		void deallocate();
	protected:
		void increase_wall_memory();
		/** A pointer to the limit of the walls array, used to
		 * determine when array is full. */
		wall **wel;
		/** The current amount of memory allocated for walls. */
		int current_wall_size;
};

class container_base : public voropp_base, public wall_list {
	public:
		/** The minimum x coordinate of the container. */
		const double ax;
		/** The maximum x coordinate of the container. */
		const double bx;
		/** The minimum y coordinate of the container. */
		const double ay;
		/** The maximum y coordinate of the container. */
		const double by;
		/** The minimum z coordinate of the container. */
		const double az;
		/** The maximum z coordinate of the container. */
		const double bz;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** A boolean value that determines if the z coordinate in
		 * periodic or not. */
		const bool zperiodic;
		/** This array holds the numerical IDs of each particle in each
		 * computational box. */
		int **id;
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		double **p;
		/** This array holds the number of particles within each
		 * computational box of the container. */
		int *co;
		/** This array holds the maximum amount of particle memory for
		 * each computational box of the container. If the number of
		 * particles in a particular box ever approaches this limit,
		 * more is allocated using the add_particle_memory() function.
		 */
		int *mem;
		/** The amount of memory in the array structure for each
		 * particle. This is set to 3 when the basic class is
		 * initialized, so that the array holds (x,y,z) positions. If
		 * the container class is initialized as part of the derived
		 * class container_poly, then this is set to 4, to also hold
		 * the particle radii. */
		const int ps;
		container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,
				int init_mem,int ps);
		~container_base();
		bool point_inside(double x,double y,double z);
		void region_count();
		/** Initialize the Voronoi cell to be the entire container. For
		 * non-periodic coordinates, this is set by the position of the
		 * walls. For periodic coordinates, the space is equally
		 * divided in either direction from the particle's initial
		 * position. That makes sense since those boundaries would be
		 * made by the neighboring periodic images of this particle. It
		 * also applies plane cuts made by any walls that have been
		 * added to the container.
		 * \param[in,out] c a reference to a voronoicell object.
		 * \param[in] (x,y,z) the position of the particle.
		 * \return False if the plane cuts applied by walls completely
		 * removed the cell, true otherwise. */
		template<class v_cell,class v_loop>
		inline bool initialize_voronoicell(v_cell &c,v_loop &vl,int &sti,int &stj,int &stk) {
			double x1,x2,y1,y2,z1,z2;
			double *pp(p[vl.ijk]+ps*vl.q);
			cux=*(pp++);cuy=*(pp++);cuz=*pp;
			cui=vl.i;cuj=vl.j;cuk=vl.k;
			if(xperiodic) {x1=-(x2=0.5*(bx-ax));sti=nx;} else {x1=ax-cux;x2=bx-cux;sti=cui;}
			if(yperiodic) {y1=-(y2=0.5*(by-ay));stj=ny;} else {y1=ay-cuy;y2=by-cuy;stj=cuj;}
			if(zperiodic) {z1=-(z2=0.5*(bz-az));stk=nz;} else {z1=az-cuz;z2=bz-cuz;stk=cuk;}
			c.init(x1,x2,y1,y2,z1,z2);
			if(!apply_walls(c,cux,cuy,cuz)) return false;
			cuijk=vl.ijk-sti-nx*(stj+ny*stk);
			return true;
		}
		/** Returns the position of the particle currently being
		 * computed relative to the computational block that it is
		 * within. It is used to select the optimal worklist entry to
		 * use.
		 * \param[out] (fx,fy,fz) the relative position. */
		inline void frac_pos(double &fx,double &fy,double &fz) {
			fx=cux-ax-boxx*cui;
			fy=cuy-ay-boxy*cuj;
			fz=cuz-az-boxz*cuk;
		}
		inline int region_index(int ei,int ej,int ek,double &qx,double &qy,double &qz) {
			if(xperiodic) {if(cui+ei<nx) {ei+=nx;qx=-(bx-ax);} else if(cui+ei>=(nx<<1)) {ei-=nx;qx=bx-ax;} else qx=0;}
			if(yperiodic) {if(cuj+ej<ny) {ej+=ny;qy=-(by-ay);} else if(cuj+ej>=(ny<<1)) {ej-=ny;qy=by-ay;} else qy=0;}
			if(zperiodic) {if(cuk+ek<nz) {ek+=nz;qz=-(bz-az);} else if(cuk+ek>=(nz<<1)) {ek-=nz;qz=bz-az;} else qz=0;}
			return cuijk+ei+nx*(ej+ny*ek);
		}
		void draw_domain_gnuplot(FILE *fp=stdout);
		inline void draw_domain_gnuplot(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_domain_gnuplot(fp);
			fclose(fp);
		}
		void draw_domain_pov(FILE *fp=stdout);
		inline void draw_domain_pov(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_domain_pov(fp);
			fclose(fp);
		}
	protected:
		double cux;
		double cuy;
		double cuz;
		int cui,cuj,cuk,cuijk;
		void add_particle_memory(int i);
		inline bool put_locate_block(int &ijk,double &x,double &y,double &z);
		inline bool put_remap(int &ijk,double &x,double &y,double &z);
};

class container : public container_base {
	public:
		container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem);
		void clear();
		void put(int n,double x,double y,double z);
		void put(voropp_order &vo,int n,double x,double y,double z);
		void import(FILE *fp=stdin);
		void import(voropp_order &vo,FILE *fp=stdin);
		inline void import(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(fp);
			fclose(fp);
		}
		inline void import(voropp_order &vo,const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(vo,fp);
			fclose(fp);
		}
		void compute_all_cells();
		double sum_cell_volumes();
		/** Dumps particle IDs and positions to a file.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_particles(v_loop &vl,FILE *fp) {
			double *pp;
			if(vl.start()) do {
				pp=p[vl.ijk]+3*vl.q;
				fprintf(fp,"%d %g %g %g\n",id[vl.ijk][vl.q],*pp,pp[1],pp[2]);
			} while(vl.inc());
		}
		/** Dumps all of the particle IDs and positions to a file.
		 * \param[in] fp a file handle to write to. */
		inline void draw_particles(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_particles(vl,fp);
		}
		/** Dumps all of the particle IDs and positions to a file.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_particles(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_particles(fp);
			fclose(fp);
		}
		/** Dumps particle positions in POV-Ray format.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_particles_pov(v_loop &vl,FILE *fp) {
			double *pp;
			if(vl.start()) do {
				pp=p[vl.ijk]+3*vl.q;
				fprintf(fp,"// id %d\nsphere{<%g,%g,%g>,r}\n",
						id[vl.ijk][vl.q],*pp,pp[1],pp[2]);
			} while(vl.inc());
		}
		/** Dumps all particle positions in POV-Ray format.
		 * \param[in] fp a file handle to write to. */
		inline void draw_particles_pov(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_particles_pov(vl,fp);
		}
		/** Dumps all particle positions in POV-Ray format.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_particles_pov(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_particles_pov(fp);
			fclose(fp);
		}
		/** Computes Voronoi cells and saves the output in gnuplot
		 * format.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_cells_gnuplot(v_loop &vl,FILE *fp) {
			voronoicell c;double *pp;
			if(vl.start()) do if(compute_cell(c,vl)) {
				pp=p[vl.ijk]+ps*vl.q;
				c.draw_gnuplot(*pp,pp[1],pp[2],fp);
			} while(vl.inc());
		}
		/** Compute all Voronoi cells and saves the output in gnuplot
		 * format.
		 * \param[in] fp a file handle to write to. */
		inline void draw_cells_gnuplot(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_cells_gnuplot(vl,fp);
		}
		/** Compute all Voronoi cells and saves the output in gnuplot
		 * format.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_cells_gnuplot(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_cells_gnuplot(fp);
			fclose(fp);
		}
		/** Computes Voronoi cells and saves the output in POV-Ray
		 * format.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_cells_pov(v_loop &vl,FILE *fp) {
			voronoicell c;double *pp;
			if(vl.start()) do if(compute_cell(c,vl)) {
				fprintf(fp,"// cell %d\n",id[vl.ijk][vl.q]);
				pp=p[vl.ijk]+ps*vl.q;
				c.draw_pov(*pp,pp[1],pp[2],fp);
			} while(vl.inc());
		}
		/** Computes all Voronoi cells and saves the output in POV-Ray
		 * format.
		 * \param[in] fp a file handle to write to. */
		inline void draw_cells_pov(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_cells_pov(vl,fp);
		}
		/** Computes all Voronoi cells and saves the output in POV-Ray
		 * format.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_cells_pov(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_cells_pov(fp);
			fclose(fp);
		}
		/** Computes the Voronoi cells and saves customized information
		 * about them.
		 * \param[in] vl the loop class to use.
		 * \param[in] format the custom output string to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void print_custom(v_loop &vl,const char *format,FILE *fp) {
			int ijk,q;double *pp;
			if(contains_neighbor(format)) {
				voronoicell_neighbor c;
				if(vl.start()) do if(compute_cell(c,vl)) {
					ijk=vl.ijk;q=vl.q;pp=p[ijk]+ps*q;
					c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],default_radius,fp);
				} while(vl.inc());
			} else {
				voronoicell c;
				if(vl.start()) do if(compute_cell(c,vl)) {
					ijk=vl.ijk;q=vl.q;pp=p[ijk]+ps*q;
					c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],default_radius,fp);
				} while(vl.inc());
			}
		}
		void print_custom(const char *format,FILE *fp=stdout);
		void print_custom(const char *format,const char *filename);
		template<class v_cell,class v_loop>
		inline bool compute_cell(v_cell &c,v_loop &vl) {
			int sti,stj,stk;
			if(!initialize_voronoicell(c,vl,sti,stj,stk)) return false;
			return vc.compute_cell(c,vl.ijk,vl.q,sti,stj,stk,cux,cuy,cuz);
		}
	private:
		voropp_compute<container> vc;
		inline void r_init(int ijk,int s) {};
		inline double r_cutoff(double lrs) {return lrs;}
		inline double r_scale(double rs,int ijk,int q) {return rs;}
		friend class voropp_compute<container>;
};

class container_poly : public container_base {
	public:
		/** The current maximum radius of any particle, used to
		 * determine when to cut off the radical Voronoi computation.
		 * */
		double max_radius;
		container_poly(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem);
		void clear();
		void put(int n,double x,double y,double z,double r);
		void put(voropp_order &vo,int n,double x,double y,double z,double r);
		void import(FILE *fp=stdin);
		void import(voropp_order &vo,FILE *fp=stdin);
		inline void import(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(fp);
			fclose(fp);
		}
		inline void import(voropp_order &vo,const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"r"));
			import(vo,fp);
			fclose(fp);
		}
		void compute_all_cells();
		double sum_cell_volumes();
		/** Dumps particle IDs, positions and radii to a file.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_particles(v_loop &vl,FILE *fp) {
			double *pp;
			if(vl.start()) do {
				pp=p[vl.ijk]+4*vl.q;
				fprintf(fp,"%d %g %g %g %g\n",id[vl.ijk][vl.q],*pp,pp[1],pp[2],pp[3]);
			} while(vl.inc());
		}
		/** Dumps all of the particle IDs, positions and radii to a
		 * file.
		 * \param[in] fp a file handle to write to. */
		inline void draw_particles(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_particles(vl,fp);
		}
		/** Dumps all of the particle IDs, positions and radii to a
		 * file.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_particles(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_particles(fp);
			fclose(fp);
		}
		/** Dumps particle positions in POV-Ray format.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_particles_pov(v_loop &vl,FILE *fp) {
			double *pp;
			if(vl.start()) do {
				pp=p[vl.ijk]+4*vl.q;
				fprintf(fp,"// id %d\nsphere{<%g,%g,%g>,%g}\n",
						id[vl.ijk][vl.q],*pp,pp[1],pp[2],pp[3]);
			} while(vl.inc());
		}
		/** Dumps all the particle positions in POV-Ray format.
		 * \param[in] fp a file handle to write to. */
		inline void draw_particles_pov(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_particles_pov(vl,fp);
		}
		/** Dumps all the particle positions in POV-Ray format.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_particles_pov(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_particles_pov(fp);
			fclose(fp);
		}
		/** Computes Voronoi cells and saves the output in gnuplot
		 * format.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_cells_gnuplot(v_loop &vl,FILE *fp) {
			voronoicell c;double *pp;
			if(vl.start()) do if(compute_cell(c,vl)) {
				pp=p[vl.ijk]+ps*vl.q;
				c.draw_gnuplot(*pp,pp[1],pp[2],fp);
			} while(vl.inc());
		}
		/** Compute all Voronoi cells and saves the output in gnuplot
		 * format.
		 * \param[in] fp a file handle to write to. */
		inline void draw_cells_gnuplot(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_cells_gnuplot(vl,fp);
		}
		/** Compute all Voronoi cells and saves the output in gnuplot
		 * format.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_cells_gnuplot(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_cells_gnuplot(fp);
			fclose(fp);
		}
		/** Computes Voronoi cells and saves the output in POV-Ray
		 * format.
		 * \param[in] vl the loop class to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void draw_cells_pov(v_loop &vl,FILE *fp) {
			voronoicell c;double *pp;
			if(vl.start()) do if(compute_cell(c,vl)) {
				fprintf(fp,"// cell %d\n",id[vl.ijk][vl.q]);
				pp=p[vl.ijk]+ps*vl.q;
				c.draw_pov(*pp,pp[1],pp[2],fp);
			} while(vl.inc());
		}
		/** Computes all Voronoi cells and saves the output in POV-Ray
		 * format.
		 * \param[in] fp a file handle to write to. */
		inline void draw_cells_pov(FILE *fp=stdout) {
			v_loop_all vl(*this);
			draw_cells_pov(vl,fp);
		}
		/** Computes all Voronoi cells and saves the output in POV-Ray
		 * format.
		 * \param[in] filename the name of the file to write to. */
		inline void draw_cells_pov(const char *filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_cells_pov(fp);
			fclose(fp);
		}
		/** Computes the Voronoi cells and saves customized information
		 * about them.
		 * \param[in] vl the loop class to use.
		 * \param[in] format the custom output string to use.
		 * \param[in] fp a file handle to write to. */
		template<class v_loop>
		void print_custom(v_loop &vl,const char *format,FILE *fp) {
			int ijk,q;double *pp;
			if(contains_neighbor(format)) {
				voronoicell_neighbor c;
				if(vl.start()) do if(compute_cell(c,vl)) {
					ijk=vl.ijk;q=vl.q;pp=p[ijk]+ps*q;
					c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],pp[3],fp);
				} while(vl.inc());
			} else {
				voronoicell c;
				if(vl.start()) do if(compute_cell(c,vl)) {
					ijk=vl.ijk;q=vl.q;pp=p[ijk]+ps*q;
					c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],pp[3],fp);
				} while(vl.inc());
			}
		}

		template<class v_cell,class v_loop>
		inline bool compute_cell(v_cell &c,v_loop &vl) {
			int sti,stj,stk;
			if(!initialize_voronoicell(c,vl,sti,stj,stk)) return false;
			return vc.compute_cell(c,vl.ijk,vl.q,sti,stj,stk,cux,cuy,cuz);
		}
		void print_custom(const char *format,FILE *fp=stdout);
		void print_custom(const char *format,const char *filename);
	private:
		voropp_compute<container_poly> vc;
		double r_rad,r_mul;
		inline void r_init(int ijk,int s) {
			r_rad=p[ijk][4*s+3];
			r_mul=1+(r_rad*r_rad-max_radius*max_radius)/((max_radius+r_rad)*(max_radius+r_rad));
			r_rad*=r_rad;
		}
		inline double r_cutoff(double lrs) {
			return r_mul*lrs;
		}
		inline double r_scale(double rs,int ijk,int q) {
			return rs+r_rad-p[ijk][4*q+3]*p[ijk][4*q+3];
		}
		friend class voropp_compute<container_poly>;
};

#endif
