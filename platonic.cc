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

	// Create a tetrahedron
	v.init(-2,2,-2,2,-2,2);
	v.plane(1,1,1);
	v.plane(1,-1,-1);
	v.plane(-1,1,-1);
	v.plane(-1,-1,1);

	v.facets();
	
	file.open("tetrahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();

	// Create a cube. Since this is the default shape
	// we don't need to do any plane cutting 
	v.init(-1,1,-1,1,-1,1);
	file.open("cube",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();
	
	v.facets();

	// Create an octahedron
	v.init(-2,2,-2,2,-2,2);
	v.nplane(1,1,1,1);
	v.nplane(-1,1,1,2);
	v.nplane(1,-1,1,3);
	v.nplane(-1,-1,1,4);
	v.nplane(1,1,-1,5);
	v.nplane(-1,1,-1,6);
	v.nplane(1,-1,-1,7);
	v.nplane(-1,-1,-1,8);

	file.open("octahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();
	cout << endl;
	v.edgeprint();
	v.relcheck();
	cout << "Vol is " << v.volume() << endl;
	v.facets();

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

	v.facets();

	// Create an icosahedron
	v.init(-2,2,-2,2,-2,2);
	v.nplane(1,1,1,1);
	v.nplane(-1,1,1,2);
	v.nplane(1,-1,1,3);
	v.nplane(-1,-1,1,4);
	v.nplane(1,1,-1,5);
	v.nplane(-1,1,-1,6);
	v.nplane(1,-1,-1,7);
	v.nplane(-1,-1,-1,8);
	v.nplane(0,phi,Phi,9);	
	v.nplane(0,phi,-Phi,10);	
	v.nplane(0,-phi,Phi,11);	
	v.nplane(0,-phi,-Phi,12);
	v.nplane(Phi,0,phi,13);
	v.nplane(Phi,0,-phi,14);
	v.nplane(-Phi,0,phi,15);
	v.nplane(-Phi,0,-phi,16);
	v.nplane(phi,Phi,0,17);
	v.nplane(phi,-Phi,0,18);
	v.nplane(-phi,Phi,0,19);
	v.nplane(-phi,-Phi,0,20);

	file.open("icosahedron",ofstream::out|ofstream::trunc);
	v.dumpgnuplot(file,0,0,0);
	file.close();

	v.facets();
}
