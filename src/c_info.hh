// Voro++,a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file c_info.hh
 * \brief Header file for c_info class for pointing at particles in the
 * container. */

#ifndef VOROPP_C_INFO_HH
#define VOROPP_C_INFO_HH

namespace voro {

/** The types of geometrical region to loop over. */
enum subset_mode {
    sphere,
    circle,
    box,
    no_check
};

struct c_info {
    /** A particle block index. */
    int ijk;
    /** The index of a particle within the block. */
    int q;
    c_info() {}
    c_info(int ijk_,int q_) : ijk(ijk_), q(q_) {}
    inline void set(int ijk_,int q_) {
        ijk=ijk_;q=q_;
    }
};

}

#endif
