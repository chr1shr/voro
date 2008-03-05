// Platonic solids example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.cc"

// Golden ratio constants
const double Phi=0.5*(1+sqrt(5.0));
const double phi=0.5*(1-sqrt(5.0));

int main() {
	voronoicell v;
	ofstream file;
/*
	// Create a tetrahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(1,1,1);
	v.plane(1,-1,-1);
	v.plane(-1,1,-1);
	v.plane(-1,-1,1);
	
	file.open("tetrahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();

	// Create a cube. Since this is the default shape
	// we don't need to do any plane cutting 
	v.init(-1,1,-1,1,-1,1);
//	v.init_octahedron(1);
	file.open("cube",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();
	v.facets();*/

	// Create an octahedron
	v.init(-2,2,-2,2,-2,2);
	v.nplane(1.8,1.8,1.8,1);
/*	v.nplane(-1,1,1,2);
	v.nplane(1,-1,1,3);
	v.nplane(-1,-1,1,4);
	v.nplane(1,1,-1,5);
	v.nplane(-1,1,-1,6);
	v.nplane(1,-1,-1,7);
	v.nplane(-1,-1,-1,8);*/

	file.open("octahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();
	cout << endl;
	v.facets();
/*
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

	file.open("dodecahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();

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

	file.open("icosahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();	*/
}
