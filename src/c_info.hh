// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file c_info.hh
 * \brief Header file for c_info class for pointing at particles in the
 * container. */

#ifndef VOROPP_C_INFO_HH
#define VOROPP_C_INFO_HH

// XXX CHR - I shifted this into its own header file because it is used by both
// the 2D and 3D iterators

namespace voro {

struct c_info {
    int ijk;
    int q;
    c_info() {}
    // XXX CHR - I've added a constructor here, which should simplify some of
    // the code you wrote. Rather than writing "c_info ci;ci.set(ijk,q);return
    // ci;" you can write "return c_info(ijk,q);".
    c_info(int ijk_,int q_) ijk(ijk_), q(q_) {}
    inline void set(int ijk_,int q_){
        ijk=ijk_; q=q_;
    }
};

// XXX CHR - I am thinking to change particle_order into std::vector<c_info>,
// which will accomplish the same thing and cut down on the classes that need
// to be specified. It will also simplify the iterator_order classes.

}

#endif
