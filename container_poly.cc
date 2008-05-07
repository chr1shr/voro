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
	initialize_voronoicell(c);

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

inline void voronoicell_neighbor::neighbor_main_allocate() {
	mne=new int*[currentvertexorder];
	ne=new int*[currentvertices];
};

inline void voronoicell_neighbor::neighbor_allocate(int i,int m) {
	mne[i]=new int[m*i];
};	

inline void voronoicell_neighbor::neighbor_deallocate(int i) {
	delete [] mne[i];
};

inline void voronoicell_neighbor::neighbor_main_deallocate() {
	delete [] mne;
	delete [] ne;
};

inline void voronoicell_neighbor::neighbor_addmemory_vertices(int i) {
	int *pp;
	pp=new int*[i];
	for(int j=0;j<currentvertices;j++) pp[j]=ne[j];
	delete [] ne;ne=pp;
};

inline void voronoicell_neighbor::neighbor_addmemory_vorder(int i) {
	int **p2;
	p2=new int*[i];
	for(j=0;j<currentvertexorder;j++) p2[j]=mne[j];
	delete [] mne;mne=p2;
};

inline void voronoicell_neighbor::neighbor_init() {
	int *q;
	q=mne[3];
	q[0]=-5;q[1]=-3;q[2]=-1;
	q[3]=-5;q[4]=-2;q[5]=-3;
	q[6]=-5;q[7]=-1;q[8]=-4;
	q[9]=-5;q[10]=-4;q[11]=-2;
	q[12]=-6;q[13]=-1;q[14]=-3;
	q[15]=-6;q[16]=-3;q[17]=-2;
	q[18]=-6;q[19]=-4;q[20]=-1;
	q[21]=-6;q[22]=-2;q[23]=-4;
	ne[0]=q;ne[1]=q+3;ne[2]=q+6;ne[3]=q+9;
	ne[4]=q+12;ne[5]=q+15;ne[6]=q+18;ne[7]=q+21;
};

inline void voronoicell_neighbor::neighbor_init_octahedron() {
	int *q;
	q=mne[4];
	q[0]=-5;q[1]=-6;q[2]=-7;q[3]=-8;
	q[4]=-1;q[5]=-2;q[6]=-3;q[7]=-4;
	q[8]=-6;q[9]=-5;q[10]=-2;q[11]=-1;
	q[12]=-8;q[13]=-7;q[14]=-4;q[15]=-3;
	q[16]=-5;q[17]=-8;q[18]=-3;q[19]=-2;
	q[20]=-7;q[21]=-6;q[22]=-1;q[23]=-4;
	ne[0]=q;ne[1]=q+4;ne[2]=q+8;ne[3]=q+12;ne[4]=q+16;ne[5]=q+20;
};	
