// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

/** Put a particle into the correct region of the container. */
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
};

void container_poly::put(int n,fpoint x,fpoint y,fpoint z) {
	put(n,x,y,z,0.5);
};

/** Import a list of particles from standard input. */
void container_poly::import(istream &is) {
	int n;fpoint x,y,z;
	fpoint r;
	is >> n >> x >> y >> z >> r;
	while(!is.eof()) {
		put(n,x,y,z,r);
		is >> n >> x >> y >> z >> r;
	}
};

/** Clears a container of particles. */
void container_poly::clear() {
	for(int ijk=0;ijk<nxyz;ijk++) co[ijk]=0;
	max_radius=0;
};

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output. */
void container_poly::vprintall(ostream &os) {
	fpoint x,y,z;
	voronoicell c;
	facets_loop l(this);
	int i,s;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			x=p[s][sz*i];y=p[s][sz*i+1];z=p[s][sz*i+2];
			compute_cell(c,s,i,x,y,z);
			os << id[s][i] << " " << x << " " << y << " " << z << " " << p[s][4*i+3];
			os << " " << c.volume();
#ifdef FACETS_NEIGHBOR
			c.neighbors(os);			
#endif
			os << endl;
		}
	}
};


/** Computes a single Voronoi cell in the container. This routine can be run by
 * the user, and it is also called multiple times by the functions vprintall,
 * vcomputeall and vdraw. */
inline void container_poly::compute_cell(voronoicell &c,int s,int i,fpoint x,fpoint y,fpoint z) {
	fpoint x1,y1,z1,x2,y2,z2,qx,qy,qz,lr=0,lrs=0,ur,urs,rs;
	int j,t;
	facets_loop l(this);
	fpoint crad=p[s][4*i+3];
	const fpoint mul=1+(crad*crad-max_radius*max_radius)/((max_radius+crad)*(max_radius+crad));
	crad*=crad;

	// Initialize the voronoi cell to be the entire container. For
	// non-periodic coordinates, this is set by the position of the walls.
	// For periodic coordinates, the space is equally divided in either
	// direction from the particle's initial position. That makes sense
	// since those boundaries would be made by the neighboring periodic
	// images of this particle. 
	if (xperiodic) x1=-(x2=0.5*(bx-ax));else {x1=ax-x;x2=bx-x;}
	if (yperiodic) y1=-(y2=0.5*(by-ay));else {y1=ay-y;y2=by-y;}
	if (zperiodic) z1=-(z2=0.5*(bz-az));else {z1=az-z;z2=bz-z;}
	c.init(x1,x2,y1,y2,z1,z2);

	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing. TODO: this can sometimes be inefficient.
	// For example, sometimes particles at the top of granular packings can
	// extend upwards by a long way, and the shells grow very big. It would
	// be better to use a box-by-box approach, but that's not
	// straightforward.
	while(lrs*mul<c.maxradsq()) {
		ur=lr+0.5;urs=ur*ur;
		t=l.init(x,y,z,ur,qx,qy,qz);
		do {
			for(j=0;j<co[t];j++) {
				x1=p[t][sz*j]+qx-x;y1=p[t][sz*j+1]+qy-y;z1=p[t][sz*j+2]+qz-z;
				rs=x1*x1+y1*y1+z1*z1;
				if (lrs-tolerance<rs&&rs<urs&&(j!=i||s!=t))
#ifdef FACETS_NEIGHBOR
					c.nplane(x1,y1,z1,rs+crad-p[t][4*j+3]*p[t][4*j+3],id[t][j]);
#else
					c.plane(x1,y1,z1,rs+crad-p[t][4*j+3]*p[t][4*j+3]);
#endif
			}
		} while ((t=l.inc(qx,qy,qz))!=-1);
		lr=ur;lrs=urs;
	}
};
