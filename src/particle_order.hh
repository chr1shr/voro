// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file particle_order.hh
 * \brief Header file for the particle_order. */

#ifndef VOROPP_PARTICLE_ORDER_HH
#define VOROPP_PARTICLE_ORDER_HH

#include <cstring>

#include <stdio.h>

#include "config.hh"

namespace voro {

/** \brief A class for storing ordering information when particles are added to
 * a container.
 *
 * When particles are added to a container class, they are sorted into an
 * internal computational grid of blocks. The particle_order class provides a
 * mechanism for remembering which block particles were sorted into. The import
 * and put routines in the container class have variants that also take a
 * particle_order class. Each time they are called, they will store the block
 * that the particle was sorted into, plus the position of the particle within
 * the block. The particle_order class can used by the c_loop_order class to
 * specifically loop over the particles that have their information stored
 * within it. */
class particle_order {
    public:
        /** A pointer to the array holding the ordering. */
        int *o;
        /** A pointer to the next position in the ordering array in which to
         * store an entry. */
        int *op;
        /** The current memory allocation for the class, set to the number of
         * entries which can be stored. */
        int size;
        /** The particle_order constructor allocates memory to store the
         * ordering information.
         * \param[in] init_size the initial amount of memory to
         *                      allocate. */
        particle_order(int init_size=init_ordering_size)
            : o(new int[init_size<<1]),op(o),size(init_size) {}
        /** The particle_order destructor frees the dynamically allocated
         * memory used to store the ordering information. */
        ~particle_order() {
            delete [] o;
        }
        /** Adds a record to the order, corresponding to the memory address of
         * where a particle was placed into the container.
         * \param[in] ijk the block into which the particle was placed.
         * \param[in] q the position within the block where the particle was
         *              placed. */
        inline void add(int ijk,int q) {
            if(op==o+size) add_ordering_memory();
            *(op++)=ijk;*(op++)=q;
        }
    private:
        inline void add_ordering_memory() {
            /*
            int *no=new int[size*4];
            memcpy(no,o,2*sizeof(int)*size);
            delete [] o;
            size<<=1;
            o=no;op=o+size;
            */

            int *no=new int[size<<2],*nop=no,*opp=o;
            while(opp<op) *(nop++)=*(opp++);
            delete [] o;
            size<<=1;o=no;op=nop;
        }

};

}

#endif
