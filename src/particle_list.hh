// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file particle_list.hh
 * \brief Header file for the particle_list_base and related classes. */

#ifndef VOROPP_PARTICLE_LIST_HH
#define VOROPP_PARTICLE_LIST_HH

#include <cstdio>

#include "c_loops.hh"
#include "container.hh"

namespace voro {

/** \brief A class for storing an arbitrary number of particles, prior to
 * setting up a container geometry.
 *
 * The particle_list class can dynamically import and store an arbitrary number
 * of particles. Once the particles have been read in, an appropriate container
 * class can be set up with the optimal grid size, and the particles can be
 * transferred. */
class particle_list_base {
    public:
        void guess_optimal(double lx,double ly,int &nx,int &ny);
        void guess_optimal(double lx,double ly,double lz,int &nx,int &ny,int &nz);
        particle_list()
        ~particle_list();
        /** Calculates and returns the total number of particles stored within
         * the class.
         * \return The number of particles. */
        inline int total_particles() {
            return (end_id-pre_id)*pre_container_chunk_size+(ch_id-*end_id);
        }
    protected:
        /** The number of doubles associated with a single particle (either
         * two, three, or four). */
        const int ps;
        void new_chunk();
        void extend_chunk_index();
        /** The size of the chunk index. */
        int index_sz;
        /** A pointer to the chunk index to store the integer particle
         * IDs. */
        int **pre_id;
        /** A pointer to the last allocated integer ID chunk. */
        int **end_id;
        /** A pointer to the end of the integer ID chunk index, used to
         * determine when the chunk index is full. */
        int **l_id;
        /** A pointer to the next available slot on the current
         * particle ID chunk. */
        int *ch_id;
        /** A pointer to the end of the current integer chunk. */
        int *e_id;
        /** A pointer to the chunk index to store the floating point
         * information associated with particles. */
        double **pre_p;
        /** A pointer to the last allocated chunk of floating point
         * information. */
        double **end_p;
        /** A pointer to the next available slot on the current
         * floating point chunk. */
        double *ch_p;
};

/** \brief A class for storing an arbitrary number of particles with three
 * floating point number records.
 *
 * The particle_list3 class is an extension of the particle_list_base class for
 * storing particles with three floating number records (either 3D, or 2D with
 * radius information). */
class particle_list2 : public particle_list_base {
    public:
        particle_list2() : particle_list_base(3) {}
        void put(int n,double x,double y,double c);
        void import(FILE *fp=stdin);
        /** Imports particles from a file.
         * \param[in] filename the name of the file to read from. */
        inline void import(const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(fp);
            fclose(fp);
        }
        void setup(container_2d &con);
        void setup(particle_order &po,con_class &con);
};

/** \brief A class for storing an arbitrary number of particles with three
 * floating point number records.
 *
 * The particle_list3 class is an extension of the particle_list_base class for
 * storing particles with three floating number records (either 3D, or 2D with
 * radius information). */
class particle_list3 : public particle_list_base {
    public:
        particle_list3() : particle_list_base(3) {}
        void put(int n,double x,double y,double c);
        void import(FILE *fp=stdin);
        /** Imports particles from a file.
         * \param[in] filename the name of the file to read from. */
        inline void import(const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(fp);
            fclose(fp);
        }
        template<class con_class>
        void setup(con_class &con);
        template<class con_class>
        void setup(particle_order &po,con_class &con);
};

/** \brief A class for storing an arbitrary number of particles with four
 * floating point number records.
 *
 * The particle_list4 class is an extension of the particle_list_base class for
 * storing particles with four floating number records, corresponding to 3D
 * particles with additional radius information. */
class particle_list3 : public particle_list_base {
    public:
        particle_list4() : particle_list_base(4) {}
        void put(int n,double x,double y,double z,double r);
        void import(FILE *fp=stdin);
        /** Imports particles from a file.
         * \param[in] filename the name of the file to read from. */
        inline void import(const char* filename) {
            FILE *fp=safe_fopen(filename,"r");
            import(fp);
            fclose(fp);
        }
        template<class con_class>
        void setup(con_class &con);
        template<class con_class>
        void setup(particle_order &po,con_class &con);
};

}

#endif
