// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_poly::put(int n,fpoint x,fpoint y,fpoint z,fpoint r) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][sz*co[i]]=x;p[i][sz*co[i]+1]=y;p[i][sz*co[i]+2]=z;p[i][sz*co[i]+3]=r;
			if (r>max_radius) max_radius=r;
			id[i][co[i]++]=n;
		}
	}
}

/** If the radius argument is not supplied to the polydisperse put() routine
 * then assume that we're using a particle with a unit diameter.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_poly::put(int n,fpoint x,fpoint y,fpoint z) {
	put(n,x,y,z,0.5);
}

/** Import a list of particles.
 * \param[in] &is an open input stream to read from. */
void container_poly::import(istream &is) {
	int n;fpoint x,y,z;
	fpoint r;
	is >> n >> x >> y >> z >> r;
	while(!is.eof()) {
		put(n,x,y,z,r);
		is >> n >> x >> y >> z >> r;
	}
}

/** An overloaded version of the import routine, that reads the standard input.
 */
inline void container_poly::import() {
	import(cin);
}

/** An overloaded version of the import routine, that reads in particles from
 * a particular file.
 * \param[in] filename the name of the file to open and import. */
inline void container_poly::import(char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	import(is);
	is.close();
}
