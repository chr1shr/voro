// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "cell.cc"
#include "neighbor.cc"
#include "unitcell.cc"

int main() {
	unitcell uc(1,0,1,1,0,1);
	uc.draw_domain();
}
