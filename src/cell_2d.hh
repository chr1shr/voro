// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file cell_2d.hh
 * \brief Header file for the voronoicell_2d and related classes. */

#ifndef VOROPP_CELL_2D_HH
#define VOROPP_CELL_2D_HH

#include <cstdio>
#include <vector>

#include "config.hh"
#include "common.hh"

namespace voro {

/** \brief A class encapsulating all the routines for storing and calculating a
 * single two-dimensional Voronoi cell. */
class voronoicell_base_2d {
    public:
        /** This holds the current size of the ed and pts arrays. If more
         * vertices are created than can fit in these arrays, then they
         * are dynamically extended using the add_memory_vertices routine. */
        int current_vertices;
        /** This sets the size of the current delete stack. */
        int current_delete_size;
        /** The total nuber of vertices in the current cell. */
        int p;
        /** A tolerance (specified as a squared length) used to identify when a
         * vertex should be treated as being exactly on a cutting plane. */
		const double tol;
        /** An array with size 2*current_vertices holding information about
         * edge connections between vertices.*/
        int *ed;
        /** An array with size 2*current_vertices for holding the positions of
         * the vertices. */
        double *pts;
        voronoicell_base_2d(double max_len_sq);
        ~voronoicell_base_2d();
        void init_base(double xmin,double xmax,double ymin,double ymax);
        void draw_gnuplot(double x,double y,FILE *fp=stdout);
        /** Outputs the edges of the Voronoi cell in Gnuplot format to an
         * output stream.
         * \param[in] (x,y) a displacement vector to be added to the cell's
         *                  position.
         * \param[in] filename the file to write to. */
        inline void draw_gnuplot(double x,double y,const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_gnuplot(x,y,fp);
            fclose(fp);
        }
        void draw_pov(double x,double y,FILE *fp=stdout);
        /** Outputs the edges of the Voronoi cell in POV-Ray format to an open
         * file stream, displacing the cell by given vector.
         * \param[in] (x,y,z) a displacement vector to be added to the cell's
         *                    position.
         * \param[in] filename the file to write to. */
        inline void draw_pov(double x,double y,const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            draw_pov(x,y,fp);
            fclose(fp);
        }
        void output_custom(const char *format,int i,double x,double y,double r,FILE *fp=stdout);
        /** Computes the Voronoi cells for all particles in the container, and
         * for each cell, outputs a line containing custom information about
         * the cell structure. The output format is specified using an input
         * string with control sequences similar to the standard C printf()
         * routine.
         * \param[in] format the format of the output lines, using control
         *                   sequences to denote the different cell statistics.
         * \param[in] i the ID of the particle associated with this Voronoi
         *              cell.
         * \param[in] (x,y) the position of the particle associated with this
         *                  Voronoi cell.
         * \param[in] r a radius associated with the particle.
         * \param[in] filename the file to write to. */
        inline void output_custom(const char *format,int i,double x,double y,double r,const char *filename) {
            FILE *fp=safe_fopen(filename,"w");
            output_custom(format,i,x,y,r,fp);
            fclose(fp);
        }
        template<class vc_class>
        bool nplane(vc_class &vc,double x,double y,double rs,int p_id);
        template<class vc_class>
        bool nplane_cut(vc_class &vc,double x,double y,double rsq,int p_id,double u,int up);
        bool plane_intersects(double x,double y,double rs);
        inline bool plane_intersects_guess(double x,double y,double rs) {
            return plane_intersects(x,y,rs);
        }
        double max_radius_squared();
        double perimeter();
        double area();
        void vertices(std::vector<double> &v);
        void output_vertices(FILE *fp=stdout);
        void output_vertices(int pr,FILE *fp=stdout);
        void vertices(double x,double y,std::vector<double> &v);
        void output_vertices(double x,double y,FILE *fp=stdout);
        void output_vertices(int pr,double x,double y,FILE *fp=stdout);
        void edge_lengths(std::vector<double> &vd);
        void normals(std::vector<double> &vd);
        void centroid(double &cx,double &cy);
        virtual void neighbors(std::vector<int> &v) {v.clear();}
    protected:
        /** Computes the distance of a Voronoi cell vertex to a plane.
         * \param[in] (x,y) the normal vector to the plane.
         * \param[in] rsq the distance along this vector of the plane.
         * \param[in] qp the index of the vertex to consider. */
        inline double pos(double x,double y,double rsq,int qp) {
            return x*pts[2*qp]+y*pts[2*qp+1]-rsq;
        }
    private:
        template<class vc_class>
        void add_memory_vertices(vc_class &vc);
        void add_memory_ds(int *&stackp);
        /** The delete stack, used to store the vertices that are deleted
         * during the plane cutting procedure. */
        int *ds;
        /** A pointer to the end of the delete stack, used to detect when it is
         * full. */
        int *stacke;
};

class voronoicell_2d : public voronoicell_base_2d {
    public:
        using voronoicell_base_2d::nplane;
        voronoicell_2d() : voronoicell_base_2d(default_length*default_length) {}
        voronoicell_2d(double max_len_sq_) : voronoicell_base_2d(max_len_sq_) {}
        template<class c_class>
        voronoicell_2d(c_class &con) : voronoicell_base_2d(con.max_len_sq) {}
        inline bool nplane(double x,double y,double rs,int p_id) {
            return nplane(*this,x,y,rs,0);
        }
        inline bool nplane(double x,double y,int p_id) {
            double rs=x*x+y*y;
            return nplane(*this,x,y,rs,0);
        }
        inline bool plane(double x,double y,double rs) {
            return nplane(*this,x,y,rs,0);
        }
        inline bool plane(double x,double y) {
            double rs=x*x+y*y;
            return nplane(*this,x,y,rs,0);
        }
        inline void init(double xmin,double xmax,double ymin,double ymax) {
            init_base(xmin,xmax,ymin,ymax);
        }
    private:
        inline void n_add_memory_vertices(int ocv) {}
        inline void n_copy(int a,int b) {}
        inline void n_set(int a,int id) {}
        friend class voronoicell_base_2d;
};

class voronoicell_neighbor_2d : public voronoicell_base_2d {
    public:
        using voronoicell_base_2d::nplane;
        voronoicell_neighbor_2d() :
            voronoicell_base_2d(default_length*default_length),
            ne(new int[init_vertices]) {}
        voronoicell_neighbor_2d(double max_len_sq_) :
            voronoicell_base_2d(max_len_sq_),
            ne(new int[init_vertices]) {}
        template<class c_class>
        voronoicell_neighbor_2d(c_class &con) :
            voronoicell_base_2d(con.max_len_sq),
            ne(new int[init_vertices]) {}
        int *ne;
        ~voronoicell_neighbor_2d() {delete [] ne;}
        inline bool nplane(double x,double y,double rs,int p_id) {
            return nplane(*this,x,y,rs,p_id);
        }
        inline bool nplane(double x,double y,int p_id) {
            double rs=x*x+y*y;
            return nplane(*this,x,y,rs,p_id);
        }
        inline bool plane(double x,double y,double rs) {
            return nplane(*this,x,y,rs,0);
        }
        inline bool plane(double x,double y) {
            double rs=x*x+y*y;
            return nplane(*this,x,y,rs,0);
        }
        void init(double xmin,double xmax,double ymin,double ymax);
        virtual void neighbors(std::vector<int> &v);
    private:
        inline void n_add_memory_vertices(int ocv);
        inline void n_copy(int a,int b) {ne[a]=ne[b];}
        inline void n_set(int a,int id) {ne[a]=id;}
        friend class voronoicell_base_2d;
};

}
#endif
