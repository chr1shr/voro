// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

// Container constructor. The first six arguments set the corners of the box to
// be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
// grid of blocks, set by the following three arguments. The next three
// arguments are booleans, which set the periodicity in each direction. The
// final argument sets the amount of memory allocated to each block.
container::container(f_point xa,f_point xb,f_point ya,f_point yb,f_point za,f_point zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi)
	: ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),zsp(zn/(zb-za)),
	nx(xn),ny(yn),nz(zn),nxy(xn*yn),nxyz(xn*yn*zn),
	xperiodic(xper), yperiodic(yper), zperiodic(zper) {
	int l;
	co=new int[nxyz];
	for(l=0;l<nxyz;l++) co[l]=0;
	mem=new int[nxyz];
	for(l=0;l<nxyz;l++) mem[l]=memi;
	id=new int*[nxyz];
	for(l=0;l<nxyz;l++) id[l]=new int[memi];
	p=new f_point*[nxyz];
#ifdef FACETS_RADICAL
	max_radius=0;
	for(l=0;l<nxyz;l++) p[l]=new f_point[4*memi];
#else
	for(l=0;l<nxyz;l++) p[l]=new f_point[3*memi];
#endif
};

// Container destructor - free memory
container::~container() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] mem;
	delete [] co;	
};

// Dumps all the particle positions and identifies to a file
void container::dump(char *filename) {
	int c,l;
	ofstream file;
	file.open(filename,ofstream::out|ofstream::trunc);
	for(l=0;l<nxyz;l++) {
		for (c=0;c<co[l];c++)
#ifdef FACETS_RADICAL			
			file << id[l][c] << " " << p[l][4*c] << " " << p[l][4*c+1] << " " << p[l][4*c+2] << " " << p[l][4*c+3] << endl;
#else
			file << id[l][c] << " " << p[l][3*c] << " " << p[l][3*c+1] << " " << p[l][3*c+2] << endl;
#endif
	}
	file.close();
};

// Put a particle into the correct region of the container
#ifdef FACETS_RADICAL
void container::put(int n,f_point x,f_point y,f_point z,f_point r) {
#else
void container::put(int n,f_point x,f_point y,f_point z) {
#endif
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) addparticlemem(i);
#ifdef FACETS_RADICAL
			p[i][4*co[i]]=x;p[i][4*co[i]+1]=y;p[i][4*co[i]+2]=z;p[i][4*co[i]+3]=r;
			if (r>max_radius) max_radius=r;
#else
			p[i][3*co[i]]=x;p[i][3*co[i]+1]=y;p[i][3*co[i]+2]=z;
#endif
			id[i][co[i]++]=n;
		}
	}
};

// Increase memory for a particular region
void container::addparticlemem(int i) {
	int *idp;f_point *pp;
	int l,nmem=2*mem[i];
	if (nmem>maxparticlemem) throw fatal_error("Absolute maximum memory allocation exceeded");
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
#ifdef FACETS_RADICAL
	pp=new f_point[4*nmem];
	for(l=0;l<4*co[i];l++) pp[l]=p[i][l];
#else
	pp=new f_point[3*nmem];
	for(l=0;l<3*co[i];l++) pp[l]=p[i][l];
#endif
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
};

// Import a list of particles from standard input
void container::import(istream &is) {
	int n;f_point x,y,z;
#ifdef FACETS_RADICAL
	f_point r;
	is >> n >> x >> y >> z >> r;
	while(!is.eof()) {
		put(n,x,y,z,r);
		is >> n >> x >> y >> z >> r;
#else
	is >> n >> x >> y >> z;
	while(!is.eof()) {
		put(n,x,y,z);
		is >> n >> x >> y >> z;
#endif
	}
};

// An overloaded version of the import routine, that reads the
// standard input
inline void container::import() {
	import(cin);
};

// An overloaded version of the import routine, that reads
// in particles from <filename>
inline void container::import(char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	import(is);
	is.close();
};

// Outputs the number of particles within each region
void container::regioncount() {
	int i,j,k,ijk=0;
	for(k=0;k<nz;k++) {
		for(j=0;j<ny;j++) {
			for(i=0;i<nx;i++) cout << "Region (" << i << "," << j << "," << k << "): " << co[ijk++] << " particles" << endl;
		}
	}
};

// Clears a container of particles
void container::clear() {
	for(int ijk=0;ijk<nxyz;ijk++) co[ijk]=0;
#ifdef FACETS_RADICAL
	max_radius=0;
#endif
};

// Computes the Voronoi cells for all particles within a box with corners
// (xmin,ymin,zmin) and (xmax,ymax,zmax), and saves the output in a format
// that can be read by gnuplot
void container::vdraw_gnuplot(char *filename,f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax) {
	f_point x,y,z,px,py,pz;
	loop l1(this);
	int i,s;
	voronoicell c;
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(i=0;i<co[s];i++) {
#ifdef FACETS_RADICAL
			x=p[s][4*i]+px;y=p[s][4*i+1]+py;z=p[s][4*i+2]+pz;
#else
			x=p[s][3*i]+px;y=p[s][3*i+1]+py;z=p[s][3*i+2]+pz;
#endif
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				compute_cell(c,s,i,x,y,z);
				c.dumpgnuplot(of,x,y,z);
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	of.close();
};

// If only a filename is supplied to vdraw_gnuplot, then assume that we are
// calculating the entire simulation region
void container::vdraw_gnuplot(char *filename) {
	vdraw_gnuplot(filename,ax,bx,ay,by,az,bz);
};

// Computes the Voronoi cells for all particles within a box with corners
// (xmin,ymin,zmin) and (xmax,ymax,zmax), and saves the output in a format
// that can be read by gnuplot
void container::vdraw_pov(char *filename,f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax) {
	f_point x,y,z,px,py,pz;
	loop l1(this);
	int i,s;
	voronoicell c;
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	of << "#declare voronoi=union{\n";
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(i=0;i<co[s];i++) {
#ifdef FACETS_RADICAL
			x=p[s][4*i]+px;y=p[s][4*i+1]+py;z=p[s][4*i+2]+pz;
#else
			x=p[s][3*i]+px;y=p[s][3*i+1]+py;z=p[s][3*i+2]+pz;
#endif
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				compute_cell(c,s,i,x,y,z);
				c.dumppov(of,x,y,z);break;
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	of << "}\n";
	of.close();
};

// If only a filename is supplied to vdraw_pov, then assume that we are
// calculating the entire simulation region
void container::vdraw_pov(char *filename) {
	vdraw_pov(filename,ax,bx,ay,by,az,bz);
};


// Computes the Voronoi volumes for all the particles, and stores the
// results according to the particle label in the f_point array bb
void container::vcomputeall(f_point *bb) {
	voronoicell c;
	loop l(this);
	int i,s;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			compute_cell(c,s,i);
			bb[id[s][i]]=c.volume();
		}
	}
};

// Prints a list of all particle labels, positions, and Voronoi volumes to the
// standard output
void container::vprintall(ostream &of) {
	f_point x,y,z;
	voronoicell c;
	loop l(this);
	int i,s;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
#ifdef FACETS_RADICAL
			x=p[s][4*i];y=p[s][4*i+1];z=p[s][4*i+2];
#else
			x=p[s][3*i];y=p[s][3*i+1];z=p[s][3*i+2];
#endif
			compute_cell(c,s,i,x,y,z);
			of << id[s][i] << " " << x << " " << y << " " << z;
#ifdef FACETS_RADICAL
			of << " " << p[s][4*i+3];
#endif
			of << " " << c.volume();
#ifdef FACETS_NEIGHBOR
			c.neighbors(of);			
#endif
			of << endl;
		}
	}
};

// An overloaded version of vprintall, which just prints to standard output
inline void container::vprintall() {
	vprintall(cout);
};

// An overloaded version of vprintall, which outputs the result to <filename>
inline void container::vprintall(char* filename) {
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	vprintall(of);
	of.close();
};

// Computes a single Voronoi cell in the container. This routine can be run by
// the user, and it is also called multiple times by the functions vprintall,
// vcomputeall and vdraw.
inline void container::compute_cell(voronoicell &c,int s,int i,f_point x,f_point y,f_point z) {
	f_point x1,y1,z1,x2,y2,z2,qx,qy,qz,lr=0,lrs=0,ur,urs,rs;
	int j,t;
	loop l(this);
#ifdef FACETS_RADICAL
	f_point crad=p[s][4*i+3];
	const f_point mul=1+(crad*crad-max_radius*max_radius)/((max_radius+crad)*(max_radius+crad));
	crad*=crad;
#endif

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
#ifdef FACETS_RADICAL
	while(lrs*mul<c.maxradsq()) {
#else
	while(lrs<c.maxradsq()) {
#endif
		ur=lr+0.5;urs=ur*ur;
		t=l.init(x,y,z,ur,qx,qy,qz);
		do {
			for(j=0;j<co[t];j++) {
#ifdef FACETS_RADICAL
				x1=p[t][4*j]+qx-x;y1=p[t][4*j+1]+qy-y;z1=p[t][4*j+2]+qz-z;
#else
				x1=p[t][3*j]+qx-x;y1=p[t][3*j+1]+qy-y;z1=p[t][3*j+2]+qz-z;
#endif
				rs=x1*x1+y1*y1+z1*z1;
				if (lrs-tolerance<rs&&rs<urs&&(j!=i||s!=t))
#ifdef FACETS_RADICAL
#ifdef FACETS_NEIGHBOR
					c.nplane(x1,y1,z1,rs+crad-p[t][4*j+3]*p[t][4*j+3],id[t][j]);
#else
					c.plane(x1,y1,z1,rs+crad-p[t][4*j+3]*p[t][4*j+3]);
#endif
#else
#ifdef FACETS_NEIGHBOR					
					c.nplane(x1,y1,z1,rs,id[t][j]);
#else
					c.plane(x1,y1,z1,rs);
#endif
#endif
			}
		} while ((t=l.inc(qx,qy,qz))!=-1);
		lr=ur;lrs=urs;
	}
};

// A overloaded version of compute_cell, that sets up the x, y, and z variables.
inline void container::compute_cell(voronoicell &c,int s,int i) {
#ifdef FACETS_RADICAL
	double x=p[s][4*i],y=p[s][4*i+1],z=p[s][4*i+2];
#else
	double x=p[s][3*i],y=p[s][3*i+1],z=p[s][3*i+2];
#endif
	compute_cell(c,s,i,x,y,z);
}

// Creates a loop object, by pulling the necesssary constants about the container
// geometry from a pointer to the current container class
loop::loop(container *q) : sx(q->bx-q->ax), sy(q->by-q->ay), sz(q->bz-q->az),
	xsp(q->xsp),ysp(q->ysp),zsp(q->zsp),
	ax(q->ax),ay(q->ay),az(q->az),
	nx(q->nx),ny(q->ny),nz(q->nz),nxy(q->nxy),nxyz(q->nxyz),
	xperiodic(q->xperiodic),yperiodic(q->yperiodic),zperiodic(q->zperiodic) {};

// Initializes a loop object, by finding all blocks which are within a distance
// r of the vector (vx,vy,vz). It returns the first block which is to be
// tested, and sets the periodic displacement vector (px,py,pz) accordingly.
inline int loop::init(f_point vx,f_point vy,f_point vz,f_point r,f_point &px,f_point &py,f_point &pz) {
	ai=myint((vx-ax-r)*xsp);
	bi=myint((vx-ax+r)*xsp);
	if (!xperiodic) {
		if(ai<0) {ai=0;if (bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if (ai>=nx) ai=nx-1;} 
	}
	aj=myint((vy-ay-r)*ysp);
	bj=myint((vy-ay+r)*ysp);
	if (!yperiodic) {
		if(aj<0) {aj=0;if (bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if (aj>=ny) aj=ny-1;} 
	}
	ak=myint((vz-az-r)*zsp);
	bk=myint((vz-az+r)*zsp);
	if (!zperiodic) {
		if(ak<0) {ak=0;if (bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if (ak>=nz) ak=nz-1;} 
	}
	i=ai;j=aj;k=ak;
	aip=ip=mymod(i,nx);apx=px=mydiv(i,nx)*sx;
	ajp=jp=mymod(j,ny);apy=py=mydiv(j,ny)*sy;
	akp=kp=mymod(k,nz);apz=pz=mydiv(k,nz)*sz;
	inc1=aip-mymod(bi,nx);
	inc2=nx*(ny+ajp-mymod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
};

// Initializes a loop object, by finding all blocks which overlap the box with
// corners (xmin,ymin,zmin) and (xmax,ymax,zmax). It returns the first block
// which is to be tested, and sets the periodic displacement vector (px,py,pz)
// accordingly.
inline int loop::init(f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax,f_point &px,f_point &py,f_point &pz) {
	ai=myint((xmin-ax)*xsp);
	bi=myint((xmax-ax)*xsp);
	if (!xperiodic) {
		if(ai<0) {ai=0;if (bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if (ai>=nx) ai=nx-1;} 
	}
	aj=myint((ymin-ay)*ysp);
	bj=myint((ymax-ay)*ysp);
	if (!yperiodic) {
		if(aj<0) {aj=0;if (bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if (aj>=ny) aj=ny-1;} 
	}
	ak=myint((zmin-az)*zsp);
	bk=myint((zmax-az)*zsp);
	if (!zperiodic) {
		if(ak<0) {ak=0;if (bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if (ak>=nz) ak=nz-1;} 
	}
	i=ai;j=aj;k=ak;
	aip=ip=mymod(i,nx);apx=px=mydiv(i,nx)*sx;
	ajp=jp=mymod(j,ny);apy=py=mydiv(j,ny)*sy;
	akp=kp=mymod(k,nz);apz=pz=mydiv(k,nz)*sz;
	inc1=aip-mymod(bi,nx);
	inc2=nx*(ny+ajp-mymod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
};

// Returns the next block to be tested in a loop, and updates the periodicity
// vector if necessary.
inline int loop::inc(f_point &px,f_point &py,f_point &pz) {
	if (i<bi) {
		i++;
		if (ip<nx-1) {ip++;s++;} else {ip=0;s+=1-nx;px+=sx;}
		return s;
	} else if (j<bj) {
		i=ai;ip=aip;px=apx;j++;
		if (jp<ny-1) {jp++;s+=inc1;} else {jp=0;s+=inc1-nxy;py+=sy;}
		return s;
	} else if (k<bk) {
		i=ai;ip=aip;j=aj;jp=ajp;px=apx;py=apy;k++;
		if (kp<nz-1) {kp++;s+=inc2;} else {kp=0;s+=inc2-nxyz;pz+=sz;}
		return s;
	} else return -1;
};

// Custom int function, that gives consistent stepping for negative numbers.
// With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
// With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).
template <class T>
inline int loop::myint(T a) {
	return a<0?int(a)-1:int(a);
};

// Custom mod function, that gives consistent stepping for negative numbers
inline int loop::mymod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
};

// Custom div function, that gives consistent stepping for negative numbers
inline int loop::mydiv(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
};
