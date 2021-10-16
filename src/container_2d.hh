// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file container_2d.hh
 * \brief Header file for the container_base_2d and related classes. */

#ifndef VOROPP_CONTAINER_2D_HH
#define VOROPP_CONTAINER_2D_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "config.hh"
#include "common.hh"
#include "rad_option.hh"
#include "particle_order.hh"
#include "cell_2d.hh"
#include "v_base_2d.hh"
#include "v_compute_2d.hh"
#include "wall.hh"

namespace voro {

class subset_info_2d;

/** \brief Class for representing a particle system in a three-dimensional
 * rectangular box.
 *
 * This class represents a system of particles in a three-dimensional
 * rectangular box. Any combination of non-periodic and periodic coordinates
 * can be used in the three coordinate directions. The class is not intended
 * for direct use, but instead forms the base of the container and
 * container_poly classes that add specialized routines for computing the
 * regular and radical Voronoi tessellations respectively. It contains routines
 * that are commonly between these two classes, such as those for drawing the
 * domain, and placing particles within the internal data structure.
 *
 * The class is derived from the wall_list class, which encapsulates routines
 * for associating walls with the container, and the voro_base class, which
 * encapsulates routines about the underlying computational grid. */
class container_base_2d : public voro_base_2d, public wall_list_2d {
    public:
        /** The minimum x coordinate of the container. */
        const double ax;
        /** The maximum x coordinate of the container. */
        const double bx;
        /** The minimum y coordinate of the container. */
        const double ay;
        /** The maximum y coordinate of the container. */
        const double by;
        /** A boolean value that determines if the x coordinate in periodic or
         * not. */
        const bool x_prd;
        /** A boolean value that determines if the y coordinate in periodic or
         * not. */
        const bool y_prd;
        /** This array holds the numerical IDs of each particle in each
         * computational box. */
        int **id;
        /** A two dimensional array holding particle positions. For the derived
         * container_poly class, this also holds particle radii. */
        double **p;
        /** This array holds the number of particles within each computational
         * box of the container. */
        int *co;
        /** This array holds the maximum amount of particle memory for each
         * computational box of the container. If the number of particles in a
         * particular box ever approaches this limit, more is allocated using
         * the add_particle_memory() function. */
        int *mem;
        /** The amount of memory in the array structure for each particle. This
         * is set to 2 when the basic class is initialized, so that the array
         * holds (x,y) positions. If the 2D container class is initialized as
         * part of the derived class container_poly_2d, then this is set to 3,
         * to also hold the particle radii. */
        const int ps;
        container_base_2d(double ax_,double bx_,double ay_,double by_,
                int nx_,int ny_,bool x_prd_,bool y_prd_,
                int init_mem,int ps_);
        ~container_base_2d();
        bool point_inside(double x,double y);
        void region_count();
        inline bool skip(int ij,int l,double x,double y) {return false;}
        template<class v_cell_2d>
        inline bool boundary_cuts(v_cell_2d &c,int ij,double x,double y) {return true;}
        /** Initializes the Voronoi cell prior to a compute_cell operation for
         * a specific particle being carried out by a voro_compute class. The
         * cell is initialized to fill the entire container. For non-periodic
         * coordinates, this is set by the position of the walls. For periodic
         * coordinates, the space is equally divided in either direction from
         * the particle's initial position. Plane cuts made by any walls that
         * have been added are then applied to the cell.
         * \param[in,out] c a reference to a voronoicell_2d object.
         * \param[in] ij the block that the particle is within.
         * \param[in] q the index of the particle within its block.
         * \param[in] (ci,cj) the coordinates of the block in the container
         *                    coordinate system.
         * \param[out] (i,j) the coordinates of the test block relative to the
         *                   voro_compute coordinate system.
         * \param[out] (x,y) the position of the particle.
         * \param[out] disp a block displacement used internally by the
         *                  compute_cell routine.
         * \return False if the plane cuts applied by walls completely removed
         * the cell, true otherwise. */
        template<class v_cell_2d>
        inline bool initialize_voronoicell(v_cell_2d &c,int ij,int q,int ci,int cj,
                int &i,int &j,double &x,double &y,int &disp) {
            double x1,x2,y1,y2,*pp=p[ij]+ps*q;
            x=*(pp++);y=*(pp++);
            if(x_prd) {x1=-(x2=0.5*(bx-ax));i=nx;} else {x1=ax-x;x2=bx-x;i=ci;}
            if(y_prd) {y1=-(y2=0.5*(by-ay));j=ny;} else {y1=ay-y;y2=by-y;j=cj;}
            c.init(x1,x2,y1,y2);
            if(!apply_walls(c,x,y)) return false;
            disp=ij-i-nx*j;
            return true;
        }
        /** Initializes parameters for a find_voronoi_cell call within the
         * voro_compute template.
         * \param[in] (ci,cj) the coordinates of the test block in the
         *                    container coordinate system.
         * \param[in] ij the index of the test block
         * \param[out] (i,j) the coordinates of the test block relative to the
         *                   voro_compute coordinate system.
         * \param[out] disp a block displacement used internally by the
         *                  find_voronoi_cell routine. */
        inline void initialize_search(int ci,int cj,int ij,int &i,int &j,int &disp) {
            i=x_prd?nx:ci;
            j=y_prd?ny:cj;
            disp=ij-i-nx*j;
        }
        /** Returns the position of a particle currently being computed
         * relative to the computational block that it is within. It is used to
         * select the optimal worklist entry to use.
         * \param[in] (x,y) the position of the particle.
         * \param[in] (ci,cj) the block that the particle is within.
         * \param[out] (fx,fy) the position relative to the block. */
        inline void frac_pos(double x,double y,double ci,double cj,
                double &fx,double &fy) {
            fx=x-ax-boxx*ci;
            fy=y-ay-boxy*cj;
        }
        /** Calculates the index of block in the container structure
         * corresponding to given coordinates.
         * \param[in] (ci,cj) the coordinates of the original block in the
         *                    current computation, relative to the container
         *                    coordinate system.
         * \param[in] (ei,ej) the displacement of the current block from the
         *                    original block.
         * \param[in,out] (qx,qy) the periodic displacement that must be added
         *                        to the particles within the computed block.
         * \param[in] disp a block displacement used internally by the
         *                 find_voronoi_cell and compute_cell routines.
         * \return The block index. */
        inline int region_index(int ci,int cj,int ei,int ej,double &qx,double &qy,int &disp) {
            if(x_prd) {if(ci+ei<nx) {ei+=nx;qx=-(bx-ax);} else if(ci+ei>=(nx<<1)) {ei-=nx;qx=bx-ax;} else qx=0;}
            if(y_prd) {if(cj+ej<ny) {ej+=ny;qy=-(by-ay);} else if(cj+ej>=(ny<<1)) {ej-=ny;qy=by-ay;} else qy=0;}
            return disp+ei+nx*ej;
        }
        void draw_domain_gnuplot(FILE *fp=stdout);
        /** Draws an outline of the domain in Gnuplot format.
         * \param[in] filename the filename to write to. */
        inline void draw_domain_gnuplot(const char* filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_domain_gnuplot(fp);
            fclose(fp);
        }
        void draw_domain_pov(FILE *fp=stdout);
        /** Draws an outline of the domain in Gnuplot format.
         * \param[in] filename the filename to write to. */
        inline void draw_domain_pov(const char* filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_domain_pov(fp);
            fclose(fp);
        }
        /** Sums up the total number of stored particles.
         * \return The number of particles. */
        inline int total_particles() {
            int tp=*co;
            for(int *cop=co+1;cop<co+nxy;cop++) tp+=*cop;
            return tp;
        }
        friend class iterator;
        class iterator;
        iterator begin();
        iterator end();
        friend class iterator_subset;
        class iterator_subset;
        iterator_subset begin(subset_info_2d& si);
        iterator_subset end(subset_info_2d& si);
        friend class iterator_order;
        class iterator_order;
        iterator_order begin(particle_order &vo);
        iterator_order end(particle_order &vo);
    protected:
        void add_particle_memory(int i);
        inline bool put_locate_block(int &ij,double &x,double &y);
        inline bool put_remap(int &ij,double &x,double &y);
        inline bool remap(int &ai,int &aj,int &ci,int &cj,double &x,double &y,int &ij);
};

/** \brief Extension of the container_base class for computing regular Voronoi
 * tessellations.
 *
 * This class is an extension of the container_base class that has routines
 * specifically for computing the regular Voronoi tessellation with no
 * dependence on particle radii. */
class container_2d : public container_base_2d, public radius_mono {
    public:
        container_2d(double ax_,double bx_,double ay_,double by_,
                 int nx_,int ny_,bool x_prd_,bool y_prd_,int init_mem, int number_thread);
        void clear();
        ~container_2d(); 
        void change_number_thread(int number_thread);
        void put(double *pt_list, int num_pt, int num_thread);
        void put(int n,double x,double y);
        void put(particle_order &vo,int n,double x,double y);
        void put_reconcile_overflow();
        void import(FILE *fp=stdin);
        void import(particle_order &vo,FILE *fp=stdin);
        /** Imports a list of particles from an open file stream into the
         * container. Entries of three numbers (Particle ID, x position, y
         * position) are searched for. If the file cannot be successfully read,
         * then the routine causes a fatal error.
         * \param[in] filename the name of the file to open and read from. */
        inline void import(const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(fp);
            fclose(fp);
        }
        /** Imports a list of particles from an open file stream into the
         * container. Entries of three numbers (Particle ID, x position, y
         * position) are searched for. In addition, the order in which
         * particles are read is saved into an ordering class. If the file
         * cannot be successfully read, then the routine causes a fatal error.
         * \param[in,out] vo the ordering class to use.
         * \param[in] filename the name of the file to open and read from. */
        inline void import(particle_order &vo,const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(vo,fp);
            fclose(fp);
        }
        void compute_all_cells();
        double sum_cell_areas();
        void draw_particles(FILE *fp=stdout);
        /** Dumps all of the particle IDs and positions to a file.
         * \param[in] filename the name of the file to write to. */
        inline void draw_particles(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_particles(fp);
            fclose(fp);
        }
        void draw_particles_pov(FILE *fp=stdout);
        /** Dumps all particle positions in POV-Ray format.
         * \param[in] filename the name of the file to write to. */
        inline void draw_particles_pov(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_particles_pov(fp);
            fclose(fp);
        }
        void draw_cells_gnuplot(FILE *fp=stdout);
        /** Computes all Voronoi cells and saves the output in Gnuplot format.
         * \param[in] filename the name of the file to write to. */
        inline void draw_cells_gnuplot(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_cells_gnuplot(fp);
            fclose(fp);
        }
        void draw_cells_pov(FILE *fp=stdout);
        /** Computes all Voronoi cells and saves the output in POV-Ray format.
         * \param[in] filename the name of the file to write to. */
        inline void draw_cells_pov(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_cells_pov(fp);
            fclose(fp);
        }
        void print_custom(const char *format,FILE *fp=stdout);
        /** Computes all the Voronoi cells and saves customized information about them.
         * \param[in] format the custom output string to use.
         * \param[in] filename the name of the file to write to. */
        inline void print_custom(const char *format,const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            print_custom(format,fp);
            fclose(fp);
        }
        bool find_voronoi_cell(double x,double y,double &rx,double &ry,int &pid);
        /** Computes the Voronoi cell for given particle.
         * \param[out] c a Voronoi cell class in which to store the computed
         *               cell.
         * \param[in] ij the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return True if the cell was computed. If the cell cannot be
         * computed, if it is removed entirely by a wall or boundary condition,
         * then the routine returns false. */
        template<class v_cell_2d>
        inline bool compute_cell(v_cell_2d &c,int ij,int q) {
            const int tn=t_num();
            return vc[tn]->compute_cell(c,ij,q,ij%nx,ij/nx);
        }
        /** Computes the Voronoi cell for a particle currently being referenced
         * by a loop class.
         * \param[out] c a Voronoi cell class in which to store the computed
         *               cell.
         * \param[in] cli the iterator class to use.
         * \return True if the cell was computed. If the cell cannot be
         * computed, if it is removed entirely by a wall or boundary condition,
         * then the routine returns false. */
        template<class v_cell_2d,class c_iter_2d>
        inline bool compute_cell(v_cell_2d &c,c_iter_2d &cli) {
            return compute_cell(c,cli->ijk,cli->q);
        }
    private:
        int nt;
        voro_compute_2d<container_2d> **vc;
        friend class voro_compute_2d<container_2d>;
        int overflow_mem_numPt;
        int overflowPtCt;
        int *ijk_m_id_overflow;
        double *p_overflow; 
};

/** \brief Extension of the container_base class for computing radical Voronoi
 * tessellations.
 *
 * This class is an extension of container_base class that has routines
 * specifically for computing the radical Voronoi tessellation that depends on
 * the particle radii. */
class container_poly_2d : public container_base_2d, public radius_poly_2d {
    public:
        container_poly_2d(double ax_,double bx_,double ay_,double by_,
                   int nx_,int ny_,bool x_prd_,bool y_prd_,int init_mem, int number_thread);
        ~container_poly_2d();
        void change_number_thread(int number_thread);
        void clear();
        void put(int n,double x,double y,double r);
        void put(double *pt_r_list, int num_pt, int num_thread);
        void put(particle_order &vo,int n,double x,double y,double r);
        void put_reconcile_overflow();
        void import(FILE *fp=stdin);
        void import(particle_order &vo,FILE *fp=stdin);
        /** Imports a list of particles from an open file stream into the
         * container_poly class. Entries of four numbers (Particle ID, x
         * position, y position, radius) are searched for. If the file cannot
         * be successfully read, then the routine causes a fatal error.
         * \param[in] filename the name of the file to open and read
         *                     from. */
        inline void import(const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(fp);
            fclose(fp);
        }
        /** Imports a list of particles from an open file stream into the
         * container_poly class. Entries of four numbers (Particle ID, x
         * position, y position, radius) are searched for. In addition, the
         * order in which particles are read is saved into an ordering class.
         * If the file cannot be successfully read, then the routine causes a
         * fatal error.
         * \param[in,out] vo the ordering class to use.
         * \param[in] filename the name of the file to open and read from. */
        inline void import(particle_order &vo,const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(vo,fp);
            fclose(fp);
        }
        void compute_all_cells();
        double sum_cell_areas();
        void draw_particles(FILE *fp=stdout);
        /** Dumps all of the particle IDs, positions and radii to a file.
         * \param[in] filename the name of the file to write to. */
        inline void draw_particles(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_particles(fp);
            fclose(fp);
        }
        void draw_particles_pov(FILE *fp=stdout);
        /** Dumps all the particle positions in POV-Ray format.
         * \param[in] filename the name of the file to write to. */
        inline void draw_particles_pov(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_particles_pov(fp);
            fclose(fp);
        }
        void draw_cells_gnuplot(FILE *fp=stdout);
        /** Compute all Voronoi cells and saves the output in Gnuplot format.
         * \param[in] filename the name of the file to write to. */
        inline void draw_cells_gnuplot(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_cells_gnuplot(fp);
            fclose(fp);
        }
        void draw_cells_pov(FILE *fp=stdout);
        /** Computes all Voronoi cells and saves the output in POV-Ray format.
         * \param[in] filename the name of the file to write to. */
        inline void draw_cells_pov(const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_cells_pov(fp);
            fclose(fp);
        }
        void print_custom(const char *format,FILE *fp=stdout);
        /** Computes all the Voronoi cells and saves customized information
         * about them.
         * \param[in] format the custom output string to use.
         * \param[in] filename the name of the file to write to. */
        inline void print_custom(const char *format,const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            print_custom(format,fp);
            fclose(fp);
        }
        /** Computes the Voronoi cell for given particle.
         * \param[out] c a Voronoi cell class in which to store the computed
         *               cell.
         * \param[in] ij the block that the particle is within.
         * \param[in] q the index of the particle within the block.
         * \return True if the cell was computed. If the cell cannot be
         * computed, if it is removed entirely by a wall or boundary condition,
         * then the routine returns false. */
        template<class v_cell_2d>
        inline bool compute_cell(v_cell_2d &c,int ij,int q) {
            const int tn=t_num();
            return vc[tn]->compute_cell(c,ij,q,ij%nx,ij/nx);
        }
        /** Computes the Voronoi cell for a particle currently being referenced
         * by a loop class.
         * \param[out] c a Voronoi cell class in which to store the computed
         *               cell.
         * \param[in] vl the loop class to use.
         * \return True if the cell was computed. If the cell cannot be
         * computed, if it is removed entirely by a wall or boundary condition,
         * then the routine returns false. */
        template<class v_cell_2d,class c_iter_2d>
        inline bool compute_cell(v_cell_2d &c,c_iter_2d &cli) {
            return compute_cell(c,cli->ijk,cli->q);
        }
        bool find_voronoi_cell(double x,double y,double &rx,double &ry,int &pid);
    private:
        double *max_r;
        int nt;
        voro_compute_2d<container_poly_2d> **vc;
        friend class voro_compute_2d<container_poly_2d>;
        int overflow_mem_numPt;
        int overflowPtCt;
        int *ijk_m_id_overflow;
        double *p_overflow; 
};

}

#endif
