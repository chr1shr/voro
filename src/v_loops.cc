// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file v_loops.cc
 * \brief Function implementations for the loop classes. */

#include "v_loops.hh"

void voropp_order::add_ordering_memory() {
	int *no(new int[size<<2]),*nop(no),*opp(o);
	while(opp<op) *(nop++)=*(opp++);
	delete [] o;
	size<<=1;o=no;op=nop;
}
