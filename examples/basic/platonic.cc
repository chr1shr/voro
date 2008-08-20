// Platonic solids example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"

// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

int main() {
	voronoicell v;

	// Create a tetrahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(1,1,1);
	v.plane(1,-1,-1);
	v.plane(-1,1,-1);
	v.plane(-1,-1,1);
	
	v.draw_gnuplot("tetrahedron.gnu",0,0,0);

	// Create a cube. Since this is the default shape
	// we don't need to do any plane cutting. 
	v.init(-1,1,-1,1,-1,1);
	v.draw_gnuplot("cube.gnu",0,0,0);
	
	// Create an octahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(1,1,1);
	v.plane(-1,1,1);
	v.plane(1,-1,1);
	v.plane(-1,-1,1);
	v.plane(1,1,-1);
	v.plane(-1,1,-1);
	v.plane(1,-1,-1);
	v.plane(-1,-1,-1);

	v.draw_gnuplot("octahedron.gnu",0,0,0);

	// Create a dodecahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(0,Phi,1);
	v.plane(0,-Phi,1);
	v.plane(0,Phi,-1);
	v.plane(0,-Phi,-1);
	v.plane(1,0,Phi);
	v.plane(-1,0,Phi);
	v.plane(1,0,-Phi);
	v.plane(-1,0,-Phi);
	v.plane(Phi,1,0);
	v.plane(-Phi,1,0);
	v.plane(Phi,-1,0);
	v.plane(-Phi,-1,0);

	v.draw_gnuplot("dodecahedron.gnu",0,0,0);

	// Create an icosahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(1,1,1);
	v.plane(-1,1,1);
	v.plane(1,-1,1);
	v.plane(-1,-1,1);
	v.plane(1,1,-1);
	v.plane(-1,1,-1);
	v.plane(1,-1,-1);
	v.plane(-1,-1,-1);
	v.plane(0,phi,Phi);	
	v.plane(0,phi,-Phi);	
	v.plane(0,-phi,Phi);	
	v.plane(0,-phi,-Phi);
	v.plane(Phi,0,phi);
	v.plane(Phi,0,-phi);
	v.plane(-Phi,0,phi);
	v.plane(-Phi,0,-phi);
	v.plane(phi,Phi,0);
	v.plane(phi,-Phi,0);
	v.plane(-phi,Phi,0);
	v.plane(-phi,-Phi,0);

	v.draw_gnuplot("icosahedron.gnu",0,0,0);

}
