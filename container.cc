// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

// If this is defined, then extra checks are built into the code to
// notify the user of memory overflows in the container grid, and in
// the Voronoi calculation
#define OVERFLOW_CHECKING

// This constant sets the maximum number of vertices for a single Voronoi cell
const int initvertices=2560;
const int initvertexorder=640;
const int init3vertices=2560;
const int initnvertices=80;
const int initdubious=2560;
const int initdeletesize=2560;
const int initdeletesize2=2560;

const int maxvertices=1048576;
const int maxvertexorder=2048;
const int maxnvertices=1048576;
const int maxdubious=1048576;
const int maxdeletesize=1048576;
const int maxdeletesize2=1048576;
const int maxparticlemem=1048576;

const double tolerance=1e-9;
const double tolerance2=2e-9;

#include "container.hh"

// Container constructor. The first six arguments set the corners of the box to
// be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
// grid of blocks, set by the following three arguments. The next three
// arguments are booleans, which set the periodicity in each direction. The
// final argument sets the amount of memory allocated to each block.
container::container(double xa,double xb,double ya,double yb,double za,double zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi)
	: ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),zsp(zn/(zb-za)),
	nx(xn),ny(yn),nz(zn),nxy(xn*yn),nxyz(xn*yn*zn),
	xperiodic(xper), yperiodic(yper), zperiodic(zper) {
	co=new int[nxyz];
	for(int l=0;l<nxyz;l++) co[l]=0;
	mem=new int[nxyz];
	for(int l=0;l<nxyz;l++) mem[l]=memi;
	id=new int*[nxyz];
	for(int l=0;l<nxyz;l++) id[l]=new int[memi];
	p=new double*[nxyz];
	for(int l=0;l<nxyz;l++) p[l]=new double[3*memi];
};

// Dumps all the particle positions and identifies to a file
void container::dump(char *filename) {
	int c,l;
	ofstream file;
	file.open(filename,ofstream::out|ofstream::trunc);
	for(l=0;l<nxyz;l++) {
		for (c=0;c<co[l];c++) file << id[l][c] << " " << p[l][3*c] << " " << p[l][3*c+1] << " " << p[l][3*c+2] << endl;
	}
	file.close();
};

// Put a particle into the correct region of the container
void container::put(int n,double x,double y,double z) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
#ifdef OVERFLOW_CHECKING
			if(co[i]==mem[i]) addparticlemem(i);
#endif
			p[i][3*co[i]]=x;
			p[i][3*co[i]+1]=y;
			p[i][3*co[i]+2]=z;
			id[i][co[i]++]=n;
		}
	}
};

// Increase memory for a particular region
void container::addparticlemem(int i) {
	cout << "memscaleup" << i << endl;
	int *idp;double *pp;
	int l,nmem=2*mem[i];
	if (nmem>maxparticlemem) throw overflow("Absolute maximum memory allocation exceeded");
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new double[3*nmem];
	for(l=0;l<3*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete id[i];id[i]=idp;
	delete p[i];p[i]=pp;
}


// Import a list of particles from standard input
void container::import() {
	int n;double x,y,z;
	cin >> n >> x >> y >> z;
	while(!cin.eof()) {
		put(n,x,y,z);
		cin >> n >> x >> y >> z;
	}
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
};

// Computes the Voronoi cells for all particles within a box with corners
// (xmin,ymin,zmin) and (xmax,ymax,zmax). The cell edges are then saved to a
// file. The final argument sets the output format, which can be either "pov"
// for the POV-Ray raytracer, or "gnuplot" for gnuplot. 
void container::vdraw(char *filename,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,out_type ot) {
	double x,y,z,x1,y1,z1,x2,y2,z2,px,py,pz,qx,qy,qz,lr,lrs,ur,urs,rs;
	loop l1(this),l2(this);
	int i,j,s,t;
	voronoicell c;
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	if (ot==pov) of << "#declare voronoi=union{\n";
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(i=0;i<co[s];i++) {
			x=p[s][3*i]+px;
			y=p[s][3*i+1]+py;
			z=p[s][3*i+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				if (xperiodic) x1=-(x2=0.5*(bx-ax));
				else {x1=ax-x;x2=bx-x;}
				if (yperiodic) y1=-(y2=0.5*(by-ay));
				else {y1=ay-y;y2=by-y;}
				if (zperiodic) z1=-(z2=0.5*(bz-az));
				else {z1=az-z;z2=bz-z;}
				c.init(x1,x2,y1,y2,z1,z2);
				lr=lrs=0;
				while(lrs<c.maxradsq()) {
					ur=lr+0.5;urs=ur*ur;
					t=l2.init(x,y,z,ur,qx,qy,qz);
					do {
						for(j=0;j<co[t];j++) {
							x1=p[t][3*j]+qx-x;
							y1=p[t][3*j+1]+qy-y;
							z1=p[t][3*j+2]+qz-z;
							rs=x1*x1+y1*y1+z1*z1;
							if ((j!=i||s!=t)&&lrs-tolerance<rs&&rs<urs+tolerance) c.plane(x1,y1,z1,rs);
						}
					} while ((t=l2.inc(qx,qy,qz))!=-1);
					lr=ur;lrs=urs;
				}
				switch(ot) {
					case pov: c.dumppov(of,x,y,z);break;
					case gnuplot: c.dumpgnuplot(of,x,y,z);
				}
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	if (ot==pov) of << "}\n";
	of.close();
};

// If only a filename and an output type are supplied to vdraw, then assume
// that we are calculating the entire simulation region
void container::vdraw(char *filename,out_type ot) {
	vdraw(filename,ax,bx,ay,by,az,bz,ot);
};

// Computes the Voronoi volumes for all the particles, and stores the
// results according to the particle label in the double array bb
void container::vcomputeall(double *bb) {
	double x,y,z,x1,y1,z1,x2,y2,z2,lr,lrs,ur,urs,rs,qx,qy,qz;
	voronoicell c;
	loop l(this);
	int i,j,s,t;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			x=p[s][3*i];y=p[s][3*i+1];z=p[s][3*i+2];
			if (xperiodic) x1=-(x2=0.5*(bx-ax));
			else {x1=ax-x;x2=bx-x;}
			if (yperiodic) y1=-(y2=0.5*(by-ay));
			else {y1=ay-y;y2=by-y;}
			if (zperiodic) z1=-(z2=0.5*(bz-az));
			else {z1=az-z;z2=bz-z;}
			c.init(x1,x2,y1,y2,z1,z2);
			lr=lrs=0;
			while(lrs<c.maxradsq()) {
				ur=lr+0.5;urs=ur*ur;
				t=l.init(x,y,z,ur,qx,qy,qz);
				do {
					for(j=0;j<co[t];j++) {
						x1=p[t][3*j]+qx-x;
						y1=p[t][3*j+1]+qy-y;
						z1=p[t][3*j+2]+qz-z;
						rs=x1*x1+y1*y1+z1*z1;
						if ((j!=i||s!=t)&&lrs-tolerance<rs&&rs<urs+tolerance) c.plane(x1,y1,z1,rs);
					}
				} while ((t=l.inc(qx,qy,qz))!=-1);
				lr=ur;lrs=urs;
			}
			bb[id[s][i]]=c.volume();
		}
	}
};

// Prints a list of all particle labels, positions, and Voronoi volumes to the
// standard output
void container::vprintall() {
	double x,y,z,x1,y1,z1,x2,y2,z2,lr,lrs,ur,urs,rs,qx,qy,qz;
	voronoicell c;
	loop l(this);
	int i,j,s,t;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			x=p[s][3*i];y=p[s][3*i+1];z=p[s][3*i+2];
			if (xperiodic) x1=-(x2=0.5*(bx-ax));
			else {x1=ax-x;x2=bx-x;}
			if (yperiodic) y1=-(y2=0.5*(by-ay));
			else {y1=ay-y;y2=by-y;}
			if (zperiodic) z1=-(z2=0.5*(bz-az));
			else {z1=az-z;z2=bz-z;}
			c.init(x1,x2,y1,y2,z1,z2);
			lr=lrs=0;
			while(lrs<c.maxradsq()) {
				ur=lr+0.5;urs=ur*ur;
				t=l.init(x,y,z,ur,qx,qy,qz);
				do {
					for(j=0;j<co[t];j++) {
						x1=p[t][3*j]+qx-x;
						y1=p[t][3*j+1]+qy-y;
						z1=p[t][3*j+2]+qz-z;
						rs=x1*x1+y1*y1+z1*z1;
						if ((j!=i||s!=t)&&lrs-tolerance<rs&&rs<urs+tolerance) c.plane(x1,y1,z1,rs);
					}
				} while ((t=l.inc(qx,qy,qz))!=-1);
				lr=ur;lrs=urs;
			}
			cout << id[s][i] << " " << x << " " << y << " " << z << " " << c.volume() << endl;
		}
	}
};


// Prints a list of all particle labels, positions, and Voronoi volumes to the
// standard output
void container::vprintall(char *filename) {
	double x,y,z,x1,y1,z1,x2,y2,z2,lr,lrs,ur,urs,rs,qx,qy,qz;
	voronoicell c;
	loop l(this);
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	int i,j,s,t;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			x=p[s][3*i];y=p[s][3*i+1];z=p[s][3*i+2];
			if (xperiodic) x1=-(x2=0.5*(bx-ax));
			else {x1=ax-x;x2=bx-x;}
			if (yperiodic) y1=-(y2=0.5*(by-ay));
			else {y1=ay-y;y2=by-y;}
			if (zperiodic) z1=-(z2=0.5*(bz-az));
			else {z1=az-z;z2=bz-z;}
			c.init(x1,x2,y1,y2,z1,z2);
			lr=lrs=0;
			while(lrs<c.maxradsq()) {
				ur=lr+0.5;urs=ur*ur;
				t=l.init(x,y,z,ur,qx,qy,qz);
				do {
					for(j=0;j<co[t];j++) {
						x1=p[t][3*j]+qx-x;
						y1=p[t][3*j+1]+qy-y;
						z1=p[t][3*j+2]+qz-z;
						rs=x1*x1+y1*y1+z1*z1;
						if ((j!=i||s!=t)&&lrs-tolerance<rs&&rs<urs+tolerance) c.plane(x1,y1,z1,rs);
					}
				} while ((t=l.inc(qx,qy,qz))!=-1);
				lr=ur;lrs=urs;
			}
			of << id[s][i] << " " << x << " " << y << " " << z << " " << c.volume() << endl;
		}
	}
	of.close();
};


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
inline int loop::init(double vx,double vy,double vz,double r,double &px,double &py,double &pz) {
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
inline int loop::init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double &px,double &py,double &pz) {
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
inline int loop::inc(double &px,double &py,double &pz) {
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

// Constructs a Voronoi cell and sets up the initial memory
voronoicell::voronoicell() :
	currentvertices(initvertices), currentvertexorder(initvertexorder),
	currentdeletesize(initdeletesize), currentdeletesize2(initdeletesize2) {
	int i;
	ds=new int[currentdeletesize];
	ds2=new int[currentdeletesize2];
	mem=new int[currentvertexorder];
	mec=new int[currentvertexorder];
	mep=new int*[currentvertexorder];
	ed=new int*[currentvertices];
	nu=new int[currentvertices];
	pts=new double[3*currentvertices];
	sure.p=pts;
	for(i=0;i<3;i++) {
		mem[i]=initnvertices;
		mep[i]=new int[initnvertices*(2*i+1)];
		mec[i]=0;
	}
	mem[3]=init3vertices;
	mep[3]=new int[init3vertices*7];
	mec[3]=0;
	for(i=4;i<currentvertexorder;i++) {
		mem[i]=initnvertices;
		mep[i]=new int[initnvertices*(2*i+1)];
		mec[i]=0;
	}
};

voronoicell::~voronoicell() {
	delete ds;
	delete ds2;
	for(int i=0;i<currentvertexorder;i++) if (mem[i]>0) delete mep[i];
	delete mem;
	delete mec;
	delete mep;
	delete ed;
	delete nu;
	delete pts;
}

// Increases the memory storage for a particular vertex order
void voronoicell::addmemory(int i) {
	int s=2*i+1;
	if(mem[i]==0) {
		mep[i]=new int[initnvertices*s];
		mem[i]=initnvertices;
		cerr << "Order " << i << " vertex memory created " << endl;
	} else {
		int j,k,*l;
		mem[i]*=2;
		if (mem[i]>maxnvertices) throw overflow("Point memory allocation exceeded absolute maximum");
		cerr << "Order " << i << " vertex memory scaled up to " << mem[i] << endl;
		l=new int[s*mem[i]];
		j=0;
		while(j<s*mec[i]) {
			k=mep[i][j+2*i];
			if(k>=0) {
				ed[k]=l+j;
			} else {
				int m;
				for(m=0;m<stack2;m++) {
					if(ed[ds2[m]]==mep[i]+j) {
						ed[ds2[m]]=l+j;
						break;
					}
				}
				if(m==stack2) throw overflow("Couldn't relocate dangling pointer");
				cerr << "Relocated dangling pointer" << endl;
			}
			for(k=0;k<=2*i;k++) {
				l[j]=mep[i][j];
				j++;
			}
		}
		delete mep[i];
		mep[i]=l;
	}
};

void voronoicell::addmemory_vertices() {
	int i=2*currentvertices,j,**ped,*pnu;
	if (i>maxvertices) throw overflow("Vertex memory allocation exceeded absolute maximum");
	cerr << "Vertex memory scaled up to " << i << endl;
	double *ppts;
	ped=new int*[i];
	for(j=0;j<currentvertices;j++) ped[j]=ed[j];
	delete ed;ed=ped;
	pnu=new int[i];
	for(j=0;j<currentvertices;j++) pnu[j]=nu[j];
	delete nu;nu=pnu;
	ppts=new double[3*i];
	for(j=0;j<3*currentvertices;j++) ppts[j]=pts[j];
	delete pts;sure.p=pts=ppts;
	currentvertices=i;
};

void voronoicell::addmemory_vorder() {
	int i=2*currentvertexorder,j,*pmem,**pmep,*pmec;
	if (i>maxvertexorder) throw overflow("Vertex order memory allocation exceeded absolute maximum");
	cerr << "Vertex order memory scaled up to " << i << endl;
	pmem=new int[i];
	for(j=0;j<currentvertexorder;j++) pmem[j]=mem[j];while(j<i) pmem[j++]=0;
	delete mem;mem=pmem;
	pmep=new int*[i];
	for(j=0;j<currentvertexorder;j++) pmep[j]=mep[j];
	delete mep;mep=pmep;
	pmec=new int[i];
	for(j=0;j<currentvertexorder;j++) pmec[j]=mec[j];while(j<i) pmec[j++]=0;
	delete mec;mec=pmec;
	currentvertexorder=i;
};

void voronoicell::addmemory_ds() {
	int i=2*currentdeletesize,j,*pds;
	if (i>maxdeletesize) throw overflow("Delete stack 1 memory allocation exceeded absolute maximum");
	cerr << "Delete stack 1 memory scaled up to " << i << endl;
	pds=new int[i];
	for(j=0;j<currentdeletesize;j++) pds[j]=ds[j];
	delete ds;ds=pds;
	currentdeletesize=i;
};

void voronoicell::addmemory_ds2() {
	int i=2*currentdeletesize2,j,*pds2;
	if (i>maxdeletesize2) throw overflow("Delete stack 2 memory allocation exceeded absolute maximum");
	cerr << "Delete stack 2 memory scaled up to " << i << endl;
	pds2=new int[i];
	for(j=0;j<currentdeletesize2;j++) pds2[j]=ds2[j];
	delete ds2;ds2=pds2;
	currentdeletesize2=i;
};

// Initializes a Voronoi cell as a rectangular box with the given dimensions
inline void voronoicell::init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;
	mec[3]=p=8;xmin*=2;xmax*=2;ymin*=2;ymax*=2;zmin*=2;zmax*=2;
	pts[0]=xmin;pts[1]=ymin;pts[2]=zmin;
	pts[3]=xmax;pts[4]=ymin;pts[5]=zmin;
	pts[6]=xmin;pts[7]=ymax;pts[8]=zmin;
	pts[9]=xmax;pts[10]=ymax;pts[11]=zmin;
	pts[12]=xmin;pts[13]=ymin;pts[14]=zmax;
	pts[15]=xmax;pts[16]=ymin;pts[17]=zmax;
	pts[18]=xmin;pts[19]=ymax;pts[20]=zmax;
	pts[21]=xmax;pts[22]=ymax;pts[23]=zmax;
	int *q=mep[3];
	q[0]=1;q[1]=4;q[2]=2;q[3]=2;q[4]=1;q[5]=0;q[6]=0;
	q[7]=3;q[8]=5;q[9]=0;q[10]=2;q[11]=1;q[12]=0;q[13]=1;
	q[14]=0;q[15]=6;q[16]=3;q[17]=2;q[18]=1;q[19]=0;q[20]=2;
	q[21]=2;q[22]=7;q[23]=1;q[24]=2;q[25]=1;q[26]=0;q[27]=3;
	q[28]=6;q[29]=0;q[30]=5;q[31]=2;q[32]=1;q[33]=0;q[34]=4;
	q[35]=4;q[36]=1;q[37]=7;q[38]=2;q[39]=1;q[40]=0;q[41]=5;
	q[42]=7;q[43]=2;q[44]=4;q[45]=2;q[46]=1;q[47]=0;q[48]=6;
	q[49]=5;q[50]=3;q[51]=6;q[52]=2;q[53]=1;q[54]=0;q[55]=7;
	ed[0]=q;ed[1]=q+7;ed[2]=q+14;ed[3]=q+21;
	ed[4]=q+28;ed[5]=q+35;ed[6]=q+42;ed[7]=q+49;
	nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=nu[6]=nu[7]=3;
};

// Initializes a Voroni cell as a regular octahedron
inline void voronoicell::init_octahedron(double l) {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;
	mec[4]=p=6;l*=2;
	pts[0]=-l;pts[1]=0;pts[2]=0;
	pts[3]=l;pts[4]=0;pts[5]=0;
	pts[6]=0;pts[7]=-l;pts[8]=0;
	pts[9]=0;pts[10]=l;pts[11]=0;
	pts[12]=0;pts[13]=0;pts[14]=-l;
	pts[15]=0;pts[16]=0;pts[17]=l;
	int *q=mep[4];
	q[0]=2;q[1]=5;q[2]=3;q[3]=4;q[4]=0;q[5]=0;q[6]=0;q[7]=0;q[8]=0;
	q[9]=2;q[10]=4;q[11]=3;q[12]=5;q[13]=2;q[14]=2;q[15]=2;q[16]=2;q[17]=1;
	q[18]=0;q[19]=4;q[20]=1;q[21]=5;q[22]=0;q[23]=3;q[24]=0;q[25]=1;q[26]=2;
	q[27]=0;q[28]=5;q[29]=1;q[30]=4;q[31]=2;q[32]=3;q[33]=2;q[34]=1;q[35]=3;
	q[36]=0;q[37]=3;q[38]=1;q[39]=2;q[40]=3;q[41]=3;q[42]=1;q[43]=1;q[44]=4;
	q[45]=0;q[46]=2;q[47]=1;q[48]=3;q[49]=1;q[50]=3;q[51]=3;q[52]=1;q[53]=5;
	ed[0]=q;ed[1]=q+9;ed[2]=q+18;ed[3]=q+27;ed[4]=q+36;ed[5]=q+45;
	nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=4;
};

// Initializes an arbitrary test object using the add_vertex() and
// relconstruct() routines
inline void voronoicell::init_test() {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;p=0;

	/*add_vertex(1,-2,-1,5,1,3);
	add_vertex(0,-1,1,2,0,5);
	add_vertex(0,1,0,4,6,3,1);
	add_vertex(1,2,-1,0,2,6);
	add_vertex(-1,2,-1,6,2,5);
	add_vertex(-1,-2,-1,1,0,4);
	add_vertex(0,3,0,3,2,4);*/

	/*add_vertex(-2,2,-1,1,4,3);
	add_vertex(2,2,-1,2,5,0);
	add_vertex(2,-2,-1,3,6,1);
	add_vertex(-2,-2,-1,0,7,2);
	add_vertex(-1,1,0,5,7,0);
	add_vertex(1,1,1,6,4,1);
	add_vertex(1,-1,1,5,2,7);
	add_vertex(-1,-1,1,6,3,4);*/

	/*add_vertex(1,-2,-1,4,3,1);
	add_vertex(-1,-2,-1,5,4,0,2);
	add_vertex(-1,2,-1,3,6,1);
	add_vertex(1,2,-1,0,5,6,2);
	add_vertex(0,-1,1,0,1,5);
	add_vertex(0,0,0,6,3,4,1);
	add_vertex(0,1,1,3,5,2);*/

	/*add_vertex(-1,-1,-1,1,3,4);
	add_vertex(1,-1,-1,5,2,0);
	add_vertex(1,1,-1,3,1,6);
	add_vertex(-1,1,-1,7,0,2);
	add_vertex(-1,-1,1,8,5,0,7);
	add_vertex(1,-1,1,8,6,1,4);
	add_vertex(1,1,1,8,7,2,5);
	add_vertex(-1,1,1,8,4,3,6);
	add_vertex(0,0,2,6,5,4,7);*/

	/*add_vertex(1,-3,-1,1,6,5);
	add_vertex(-1,-3,-1,2,6,0);
	add_vertex(-3,0,-1,3,8,7,1);
	add_vertex(-1,3,-1,4,9,2);
	add_vertex(1,3,-1,5,9,3);
	add_vertex(3,0,-1,0,7,8,4);
	add_vertex(0,-2,1,0,1,7);
	add_vertex(0,-1,0,5,6,2,8);
	add_vertex(0,1,0,5,7,2,9);
	add_vertex(0,2,1,4,8,3);*/

	/*add_vertex(-1,-3,-1,12,8,1,7);
	add_vertex(1,-3,-1,0,8,12,2);
	add_vertex(3,-1,-1,1,13,9,3);
	add_vertex(3,1,-1,2,9,13,4);
	add_vertex(1,3,-1,3,14,10,5);
	add_vertex(-1,3,-1,4,10,14,6);
	add_vertex(-3,1,-1,5,15,11,7);
	add_vertex(-3,-1,-1,6,11,15,0);
	add_vertex(0,-2,1,12,1,0);
	add_vertex(2,0,1,13,3,2);
	add_vertex(0,2,1,5,4,14);
	add_vertex(-2,0,1,6,15,7);
	add_vertex(0,-1,0.5,16,1,8,0);
	add_vertex(1,0,0.5,16,3,9,2);
	add_vertex(0,1,0.5,16,5,10,4);
	add_vertex(-1,0,0.5,16,7,11,6);
	add_vertex(0,0,0,14,13,12,15);*/

	/*add_vertex(2,-3,-1,1,4,3);
	add_vertex(-2,-3,-1,2,4,0);
	add_vertex(-2,3,-1,3,7,1);
	add_vertex(2,3,-1,0,6,2);
	add_vertex(0,-2,0,1,5,0);
	add_vertex(0,1,0,7,6,4);
	add_vertex(1,2,1,3,5,7);
	add_vertex(-1,2,1,2,6,5);*/

	/*add_vertex(3,-2,-1,1,3,2);
	add_vertex(-3,-2,-1,2,4,0);
	add_vertex(0,4,-1,0,5,1);
	add_vertex(1.5,-1,0,6,0,7);
	add_vertex(-1.5,-1,0,8,7,1);
	add_vertex(0,2,0,2,6,8);
	add_vertex(0.75,0.5,0,5,3,9);
	add_vertex(0,-1,0,9,3,4);
	add_vertex(-0.75,0.5,0,5,9,4);
	add_vertex(0,0,1,6,7,8);*/

	add_vertex(0,0,0,2,1,3);
	add_vertex(1,0,1,0,2,3);
	add_vertex(1,1,0,1,0,3);
	add_vertex(2,0,0,0,1,2,4,6);
	add_vertex(3,1,0,5,8,6,3);
	add_vertex(3,2,0,4);
	add_vertex(4,0,0,8,7,3,4);
	add_vertex(5,0,0,6);
	add_vertex(4,1,0,4,6);

	relconstruct();
};

// Adds an order 1 vertex to the memory structure, and specifies its edge
void voronoicell::add_vertex(double x,double y,double z,int a) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=1;
	if (mem[1]=mec[1]) addmemory(1);
	int *q=mep[1]+3*mec[1]++;ed[p]=q;
	q[0]=a;q[2]=p++;
};

// Adds an order 2 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(double x,double y,double z,int a,int b) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=2;
	if (mem[2]=mec[2]) addmemory(2);
	int *q=mep[2]+5*mec[2]++;ed[p]=q;
	q[0]=a;q[1]=b;q[4]=p++;
};

// Adds an order 3 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(double x,double y,double z,int a,int b,int c) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=3;
	if (mem[3]=mec[3]) addmemory(3);
	int *q=mep[3]+7*mec[3]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[6]=p++;
};

// Adds an order 4 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(double x,double y,double z,int a,int b,int c,int d) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=4;
	if (mem[4]=mec[4]) addmemory(4);
	int *q=mep[4]+9*mec[4]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[8]=p++;
};

// Adds an order 5 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(double x,double y,double z,int a,int b,int c,int d,int e) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=5;
	if (mem[5]=mec[5]) addmemory(5);
	int *q=mep[5]+11*mec[5]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[4]=e;q[10]=p++;
};

// Checks that the relational table of the Voronoi cell is accurate, and prints
// out any errors. This algorithm is O(p), so running it every time the plane
// routine is called will result in a significant slowdown.
inline void voronoicell::relcheck() {
	int i,j;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if (ed[ed[i][j]][ed[i][nu[i]+j]]!=i) cout << "Relational error at point " << i << ", edge " << j << "." << endl;
		}
	}
};

// This routine checks for any two vertices that are connected by more than one
// edge. The plane algorithm is designed so that this should not happen, so any
// occurrences are most likely errors. Note that the routine is O(p), so
// running it every time the plane routine is called will result in a significant
// slowdown. 
inline void voronoicell::duplicatecheck() {
	int i,j,k;
	for(i=0;i<p;i++) {
		for(j=1;j<nu[i];j++) {
			for(k=0;k<j;k++) {
				if (ed[i][j]==ed[i][k]) cout << "Duplicate edges: (" << i << "," << j << ") and (" << i << "," << k << ") [" << ed[i][j] << "]" << endl;
			}
		}
	}
};

// Constructs the relational table if the edges have been specified
inline void voronoicell::relconstruct() {
	int i,j,k,l;
	for(i=0;i<p;i++) for(j=0;j<nu[i];j++) {
		k=ed[i][j];
		l=0;
		while(ed[k][l]!=i) {
			l++;
			if (l==nu[k]) throw overflow("Relation table construction failed");
		}
		ed[i][nu[i]+j]=l;
	}
};

// Cuts the Voronoi cell by a particle whose center is at a separation of
// (x,y,z) from the cell center. The value of rsq should be initially set to
// x*x+y*y+z*z.
bool voronoicell::plane(double x,double y,double z,double rsq) {
	int count=0,i,j,k,up=0,lp=0,tp,cp,qp=1,rp,stack=0;stack2=0;
	int us,ls,qs,iqs,cs,uw,qw,lw,tw;
	int *edp,*emp;
	double u,l,t,r,q;bool complicatedsetup=false,newdoubleedge=false,doubleedge=false;

	//Initialize the safe testing routine
	sure.init(x,y,z,rsq);

	//Test approximately sqrt(n)/4 points for their proximity to the plane
	//and keep the one which is closest
	uw=sure.test(up,u);t=abs(u);
	tw=qp=1;rp=p>>3;
	while(tw<rp) {
		qw=sure.test(qp,q);
		r=abs(q);
		if(r<t) {up=qp;u=q;t=r;uw=qw;}
		tw+=qp++;
	}
	lp=up;l=u;

	// Starting from an initial guess, we now move from vertex to vertex,
	// to try and find an edge which intersects the cutting plane,
	// or a vertex which is on the plane
	try {
		if(uw==1) {
			// The test point is within the cutting space
			do {
				// If we have been around this loop more times
				// than there are points, there's a floating
				// point problem, so we'll bail out 
				if (++count>=p) throw true;
				
				// Test all the neighbors of the current point
				// and find the one which is closest to the
				// plane
				u=l;up=lp;uw=lw;
				for(i=0;i<nu[up];i++) {
					tp=ed[up][i];
					tw=sure.test(tp,t);
					if(t<l) {l=t;lw=tw;lp=tp;us=i;}
				}

				// If we couldn't find a point and the object
				// is convex, then the whole cell must be
				// within the cutting space, so there's nothing
				// left
				if (lp==up) {
					cerr << "Failed to find intersection" << endl;
					return false;
				}
			} while (lw==1);
			ls=ed[up][nu[up]+us];

			// If the last point in the iteration is within the
			// plane, we need to do the complicated setup
			// routine. Otherwise, we use the regular iteration.
			if (lw==0) {
				up=lp;
				complicatedsetup=true;
			} else complicatedsetup=false;
		} else if (uw==-1) {
			// The test point is outside of the cutting space
			do {
				// If we have been around this loop more times
				// than there are points, there's a floating
				// point problem, so we'll bail out 
				if (++count>=p) throw true;
				
				// Test all the neighbors of the current point
				// and find the one which is closest to the
				// plane
				l=u;lp=up;lw=uw;
				for(i=0;i<nu[lp];i++) {
					tp=ed[lp][i];
					tw=sure.test(tp,t);
					if(t>u) {u=t;uw=tw;up=tp;ls=i;}
				}

				// If we couldn't find a point and the object
				// is convex, then the whole cell must be
				// outside the cutting space, so it's not
				// intersected at all 
				if (up==lp) return true;
			} while (uw==-1);
			us=ed[lp][nu[lp]+ls];
			complicatedsetup=(uw!=1);
		} else {
			// Our original test point was on the plane, so we
			// automatically head for the complicated setup
			// routine
			complicatedsetup=true;
		}
	}
	catch(bool except) {
		// This routine is a fall-back, in case floating point errors
		// cause the usual search routine to fail. In the fall-back
		// routine, we just test every edge to find one straddling
		// the plane.
		cerr << "Bailed out of convex calculation\n";
		for(qp=0;qp<p;qp++) {
			qw=sure.test(qp,q);
			if (qw==1) {

				// The point is inside the cutting space. Now
				// see if we can find a neighbor which isn't.
				for(us=0;us<nu[qp];us++) {
					lp=ed[qp][us];
					if(lp<qp) {
						lw=sure.test(lp,l);
						if (lw!=1) break;
					}
				}
				if(us<nu[qp]) {
					up=qp;
					if(lw==0) {
						complicatedsetup=true;
					} else {
						complicatedsetup=false;
						u=q;
						ls=ed[up][nu[up]+us];
					}
					break;
				}
			} else if (qw==-1) {

				// The point is outside the cutting space. See
				// if we can find a neighbor which isn't.
				for(ls=0;ls<nu[qp];ls++) {
					up=ed[qp][ls];
					if(up<qp) {
						uw=sure.test(up,u);
						if (uw!=-1) break;
					}
				}
				if(ls<nu[qp]) {
					if(uw==0) {
						up=qp;
						complicatedsetup=true;
					} else {
						complicatedsetup=false;
						lp=qp;l=q;
						us=ed[lp][nu[lp]+ls];
					}
					break;
				}
			} else {
				
				// The point is in the plane, so we just
				// proceed with the complicated setup routine
				up=qp;
				complicatedsetup=true;
				break;
			}
		}
		if(qp==p) return qw==-1?true:false;
	}

	// We're about to add the first point of the new facet. In either
	// routine, we have to add a point, so first check there's space for
	// it.
	if(p==currentvertices) addmemory_vertices();

	if (complicatedsetup) {

		// The search algorithm found a point which is on the cutting
		// plane. We leave that point in place, and create a new one at
		// the same location.
		pts[3*p]=pts[3*up];
		pts[3*p+1]=pts[3*up+1];
		pts[3*p+2]=pts[3*up+2];
		
		// Search for a collection of edges of the test vertex which
		// are outside of the cutting space. Begin by testing the
		// zeroth edge. 
		i=0;
		lp=ed[up][0];
		lw=sure.test(lp,l);
		if(lw!=-1) {

			// The first edge is either inside the cutting space,
			// or lies within the cutting plane. Test the edges
			// sequentially until we find one that is outside.
			rp=lw;
			do {
				i++;

				// If we reached the last edge with no luck
				// then all of the vertices are inside
				// or on the plane, so the cell is completely
				// deleted
				if (i==nu[up]) return false;
				lp=ed[up][i];
				lw=sure.test(lp,l);
			} while (lw!=-1);
			j=i+1;

			// We found an edge outside the cutting space. Keep
			// moving through these edges until we find one that's
			// inside or on the plane.
			while(j<nu[up]) {
				lp=ed[up][j];
				lw=sure.test(lp,l);
				if (lw!=-1) break;
				j++;
			}

			// Compute the number of edges for the new vertex. In
			// general it will be the number of outside edges
			// found, plus two. But we need to recognize the
			// special case when all but one edge is outside, and
			// the remaining one is on the plane. For that case we
			// have to reduce the edge count by one to prevent
			// doubling up.
			if(j==nu[up]&&i==1&&rp==0) {
				nu[p]=nu[up];
				doubleedge=true;
			} else nu[p]=j-i+2;
			k=1;

			// Add memory for the new vertex if needed, and
			// initialize
			while (nu[p]>=currentvertexorder) addmemory_vorder();
			if (mec[nu[p]]==mem[nu[p]]) addmemory(nu[p]);
			ed[p]=mep[nu[p]]+(2*nu[p]+1)*mec[nu[p]]++;
			ed[p][2*nu[p]]=p;

			// Copy the edges of the original vertex into the new
			// one. Delete the edges of the original vertex, and
			// update the relational table.
			us=vor_down(i,up);
			while(i<j) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			qs=i==nu[up]?0:i;

		} else {

			// In this case, the zeroth edge is outside the cutting
			// plane. Begin by searching backwards from the last
			// edge until we find an edge which isn't outside.
			i=nu[up]-1;
			lp=ed[up][i];
			lw=sure.test(lp,l);
			while(lw==-1) {
				i--;

				// If i reaches zero, then we have a point in
				// the plane all of whose edges are outside
				// the cutting space, so we just exit
				if (i==0) return true;
				lp=ed[up][i];
				lw=sure.test(lp,l);
			}

			// Now search forwards from zero
			j=1;
			qp=ed[up][j];
			qw=sure.test(qp,q);
			while(qw==-1) {
				j++;
				qp=ed[up][j];
				qw=sure.test(qp,l);
			}

			// Compute the number of edges for the new vertex. In
			// general it will be the number of outside edges
			// found, plus two. But we need to recognize the
			// special case when all but one edge is outside, and
			// the remaining one is on the plane. For that case we
			// have to reduce the edge count by one to prevent
			// doubling up.
			if (i==j&&qw==0) {
				doubleedge=true;
				nu[p]=nu[up];
			} else {
				nu[p]=nu[up]-i+j+1;
			}

			// Add memory to store the vertex if it doesn't exist
			// already
			k=1;
			while(nu[p]>=currentvertexorder) addmemory_vorder();
			if (mec[nu[p]]==mem[nu[p]]) addmemory(nu[p]);

			// Copy the edges of the original vertex into the new
			// one. Delete the edges of the original vertex, and
			// update the relational table.
			ed[p]=mep[nu[p]]+(2*nu[p]+1)*mec[nu[p]]++;
			ed[p][2*nu[p]]=p;
			us=i++;
			while(i<nu[up]) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			i=0;
			while(i<j) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			qs=j;
		}
		
		// Add this point to the auxiliary delete stack
		if (stack2==currentdeletesize2) addmemory_ds2();
		ds2[stack2++]=up;

		// Look at the edges on either side of the group that was
		// detected. We're going to commence facet computation by
		// moving along one of them. We are going to end up coming back
		// along the other one.
		cs=k;
		qp=up;q=u;
		i=ed[up][us];
		us=ed[up][nu[up]+us];
		up=i;
		ed[qp][2*nu[qp]]=-p;

	} else {
		// The search algorithm found an intersected edge between the
		// points lp and up. Create a new vertex between them which
		// lies on the cutting plane. Since u and l differ by at least
		// the tolerance, this division should never screw up.
		if (stack==currentdeletesize) addmemory_ds();
		ds[stack++]=up;
		r=1/(u-l);
		pts[3*p]=(pts[3*lp]*u-pts[3*up]*l)*r;
		pts[3*p+1]=(pts[3*lp+1]*u-pts[3*up+1]*l)*r;
		pts[3*p+2]=(pts[3*lp+2]*u-pts[3*up+2]*l)*r;

		// This point will always have three edges. Connect one of them
		// to lp.
		nu[p]=3;
		if (mec[3]==mem[3]) addmemory(3);
		ed[p]=mep[3]+7*mec[3]++;
		ed[p][6]=p;
		ed[up][us]=-1;
		ed[lp][ls]=p;
		ed[lp][nu[lp]+ls]=1;
		ed[p][1]=lp;
		ed[p][nu[p]+1]=ls;
		cs=2;

		// Set the direction to move in
		qs=vor_up(us,up);
		qp=up;q=u;
	}

	// When the code reaches here, we have initialized the first point, and
	// we have a direction for moving it to construct the rest of the facet
	cp=p;rp=p;p++;
	while(qp!=up||qs!=us) {

		// We're currently tracing round an intersected facet. Keep
		// moving around it until we find a point or edge which
		// intersects the plane. 
		lp=ed[qp][qs];
		lw=sure.test(lp,l);
		
		if (lw==1) {

			// The point is still in the cutting space. Just add it
			// to the delete stack and keep moving.
			if (stack==currentdeletesize) addmemory_ds();
			qs=vor_up(ed[qp][nu[qp]+qs],lp);
			qp=lp;
			q=l;
			ds[stack++]=qp;
		
		} else if (lw==-1) {

			// The point is outside of the cutting space, so we've
			// found an intersected edge. Introduce a regular point
			// at the point of intersection. Connect it to the
			// point we just tested. Also connect it to the previous
			// new point in the facet we're constructing.
			if(p==currentvertices) addmemory_vertices();
			r=1/(q-l);
			pts[3*p]=(pts[3*lp]*q-pts[3*qp]*l)*r;
			pts[3*p+1]=(pts[3*lp+1]*q-pts[3*qp+1]*l)*r;
			pts[3*p+2]=(pts[3*lp+2]*q-pts[3*qp+2]*l)*r;
			nu[p]=3;
			if (mec[3]==mem[3]) addmemory(3);
			ed[p]=mep[3]+7*mec[3]++;
			ed[p][6]=p;
			ls=ed[qp][qs+nu[qp]];
			ed[lp][ls]=p;
			ed[lp][nu[lp]+ls]=1;
			ed[p][1]=lp;
			ed[p][0]=cp;
			ed[p][nu[p]+1]=ls;
			ed[p][nu[p]]=cs;
			ed[cp][cs]=p;
			ed[cp][nu[cp]+cs]=0;
			ed[qp][qs]=-1;
			qs=vor_up(qs,qp);
			cp=p++;
			cs=2;
		} else {

			// We've found a point which is on the cutting plane.
			// We're going to introduce a new point right here, but
			// first we need to figure out the number of edges it
			// has.
			if(p==currentvertices) addmemory_vertices();
			
			// If the previous vertex detected a double edge, our
			// new vertex will have one less edge.
			k=doubleedge?0:1;
			qs=ed[qp][nu[qp]+qs];
			qp=lp;
			iqs=qs;

			// Start testing the edges of the current point until
			// we find one which isn't outside the cutting space
			do {
				k++;
				qs=vor_up(qs,qp);
				lp=ed[qp][qs];
				lw=sure.test(lp,l);
			} while (lw==-1);
			
			// Now we need to find out whether this marginal vertex
			// we are on has been visited before, because if that's
			// the case, we need to add vertices to the existing
			// new vertex, rather than creating a fresh one. We also
			// need to figure out whether we're in a case where we
			// might be creating a duplicate edge.
			j=-ed[qp][2*nu[qp]];
	 		if(qp==up&&qs==us) {

				// If we're heading into the final part of the
				// new facet, then we never worry about the
				// duplicate edge calculation.
				newdoubleedge=false;
				if(j>0) k+=nu[j];
			} else {
				if (j>0) {

					// This vertex was visited before, so
					// count those vertices to the ones we
					// already have.
					k+=nu[j];
					
					// The only time when we might make a
					// duplicate edge is if the point we're
					// going to move to next is also a
					// marginal point, so test for that
					// first.
					if(lw==0) {

						// Now see whether this marginal point
						// has been visited before.
						i=-ed[lp][2*nu[lp]];
						if(i>0) {

							// Now see if the last edge of that other
							// marginal point actually ends up here.
							if(ed[i][nu[i]-1]==j) {
								newdoubleedge=true;
								k-=1;
							} else newdoubleedge=false;
						} else {

							// That marginal point hasn't been visited
							// before, so we probably don't have to worry
							// about duplicate edges, except in the
							// case when that's the way into the end
							// of the facet, because that way always creates
							// an edge.
							if (j==rp&&lp==up&&ed[qp][nu[qp]+qs]==us) {
								newdoubleedge=true;
								k-=1;
							} else newdoubleedge=false;
						}
					} else newdoubleedge=false;
				} else {

					// The vertex hasn't been visited
					// before, but let's see if it's
					// marginal
					if(lw==0) {

						// If it is, we need to check
						// for the case that it's a
						// small branch, and that we're
						// heading right back to where
						// we came from
						i=-ed[lp][2*nu[lp]];
						if(i==cp) {
							newdoubleedge=true;
							k-=1;
						} else newdoubleedge=false;
					} else newdoubleedge=false;
				}
			}
			
			// k now holds the number of edges of the new vertex
			// we are forming. Add memory for it if it doesn't exist
			// already.
			while(k>=currentvertexorder) addmemory_vorder();
			if (mec[k]==mem[k]) addmemory(k);
			
			// Now create a new vertex with order k, or augment
			// the existing one.
			if(j>0) {

				// If we're augmenting a vertex but we don't
				// actually need any more edges, just skip this
				// routine to avoid memory confusion
				if(nu[j]!=k) {

					// Allocate memory and copy the edges
					// of the previous instance into it
					edp=mep[k]+(2*k+1)*mec[k]++;
					i=0;
					while(i<nu[j]) {
						edp[i]=ed[j][i];
						edp[k+i]=ed[j][nu[j]+i];
						i++;
					}
					edp[2*k]=j;

					// Remove the previous instance with
					// fewer vertices from the memory
					// structure
					mec[nu[j]]--;
					emp=mep[nu[j]]+(2*nu[j]+1)*mec[nu[j]];
					if(emp!=ed[j]) {
						for(lw=0;lw<=2*nu[j];lw++) ed[j][lw]=emp[lw];
						ed[emp[2*nu[j]]]=ed[j];
					}
					ed[j]=edp;
				} else i=nu[j];
			} else {

				// Allocate a new vertex of order k
				ed[p]=mep[k]+(2*k+1)*mec[k]++;
				ed[p][2*k]=p;
				if (stack2==currentdeletesize2) addmemory_ds2();
				ds2[stack2++]=qp;
				pts[3*p]=pts[3*qp];
				pts[3*p+1]=pts[3*qp+1];
				pts[3*p+2]=pts[3*qp+2];
				ed[qp][2*nu[qp]]=-p;
				j=p++;
				i=0;
			}
			nu[j]=k;

			// Unless the previous case was a double edge, connect
			// the first available edge of the new vertex to the 
			// last one in the facet
			if (!doubleedge) {
				ed[j][i]=cp;
				ed[j][nu[j]+i]=cs;
				ed[cp][cs]=j;
				ed[cp][nu[cp]+cs]=i;
				i++;
			}

			// Copy in the edges of the underlying vertex,
			// and do one less if this was a double edge 
			qs=iqs;
			while(i<(newdoubleedge?k:k-1)) {
				qs=vor_up(qs,qp);
				lp=ed[qp][qs];ls=ed[qp][nu[qp]+qs];
				ed[j][i]=lp;
				ed[j][nu[j]+i]=ls;
				ed[lp][ls]=j;
				ed[lp][nu[lp]+ls]=i;
				ed[qp][qs]=-1;
				i++;
			}
			qs=vor_up(qs,qp);
			cs=i;
			cp=j;

			// Update the doubleedge flag, to pass it
			// to the next instance of this routine
			doubleedge=newdoubleedge;
		}
	}

	// Connect the final created vertex to the initial one
	ed[cp][cs]=rp;
	ed[rp][0]=cp;
	ed[cp][nu[cp]+cs]=0;
	ed[rp][nu[rp]+0]=cs;

	// Delete points: first, remove any duplicates
	i=0;
	while(i<stack) {
		j=ds[i];
		if(ed[j][nu[j]]!=-1) {
			ed[j][nu[j]]=-1;
			i++;
		} else ds[i]=ds[--stack];
	}
	
	// Add the points in the auxiliary delete stack,
	// and reset their back pointers
	for(i=0;i<stack2;i++) {
		j=ds2[i];
		ed[j][2*nu[j]]=j;
		if(ed[j][nu[j]]!=-1) {
			ed[j][nu[j]]=-1;
			if (stack==currentdeletesize) addmemory_ds();
			ds[stack++]=j;
		}
	}
	
	// Scan connections and add in extras
	for(i=0;i<stack;i++) {
		cp=ds[i];
		for(j=0;j<nu[cp];j++) {
			qp=ed[cp][j];
			if(qp!=-1) {
				if (ed[qp][nu[qp]]!=-1) {
					if (stack==currentdeletesize) addmemory_ds();
					ds[stack++]=qp;
					ed[qp][nu[qp]]=-1;
				}
			}
		}
	}

	// Delete them from the array structure
	while(stack>0) {
		while(ed[--p][nu[p]]==-1) {
			j=nu[p];
			mec[j]--;
			for(i=0;i<=2*j;i++) ed[p][i]=(mep[j]+(2*j+1)*mec[j])[i];
			ed[ed[p][2*j]]=ed[p];
		}
		qp=ds[--stack];
		if (qp<p) {
			// Vertex management
			pts[3*qp]=pts[3*p];
			pts[3*qp+1]=pts[3*p+1];
			pts[3*qp+2]=pts[3*p+2];

			// Memory management
			j=nu[qp];
			mec[j]--;
			for(i=0;i<=2*j;i++) ed[qp][i]=(mep[j]+(2*j+1)*mec[j])[i];
			ed[ed[qp][2*j]]=ed[qp];

			// Edge management
			ed[qp]=ed[p];
			nu[qp]=nu[p];
			for(i=0;i<nu[qp];i++) {
				if (ed[qp][i]==-1) throw overflow("fishy");
				ed[ed[qp][i]][ed[qp][nu[qp]+i]]=qp;
			}
			ed[qp][2*nu[qp]]=qp;
		} else p++;
	}

	// Check for any vertices of zero order
	if (mec[0]>0) throw overflow("Zero order vertex formed");

	// Collapse any order 2 vertices and exit	
	return collapseorder2();
};

// During the creation of a new facet in the plane routine, it is possible
// that some order 2 vertices may arise. This routine removes them.
// Suppose an order 2 vertex joins c and d. If there's a edge between
// c and d already, then the order 2 vertex is just removed; otherwise,
// the order 2 vertex is removed and c and d are joined together directly.
// It is possible this process will create order 2 or order 1 vertices,
// and the routine is continually run until all of them are removed.
inline bool voronoicell::collapseorder2() {
	if(!collapseorder1()) return false;
	int a,b,i,j,k,l;
	while(mec[2]>0) {

		// Pick a order 2 vertex and read in its edges
		i=--mec[2];
		j=mep[2][5*i];k=mep[2][5*i+1];
		if (j==k) {
			cerr << "Order two vertex joins itself" << endl;
			return false;
		}

		// Scan the edges of j to see if joins k
		for(l=0;l<nu[j];l++) {
			if(ed[j][l]==k) break;
		}

		// If j doesn't already join k, join them together.
		// Otherwise delete the connection to the current
		// vertex from j and k.
		a=mep[2][5*i+2];b=mep[2][5*i+3];i=mep[2][5*i+4];
		if(l==nu[j]) {
			ed[j][a]=k;
			ed[k][b]=j;
			ed[j][nu[j]+a]=b;
			ed[k][nu[k]+b]=a;
		} else {
			if (!delete_connection(j,a)) return false;
			if (!delete_connection(k,b)) return false;
		}

		// Compact the memory
		--p;
		if(p!=i) {
			pts[3*i]=pts[3*p];
			pts[3*i+1]=pts[3*p+1];
			pts[3*i+2]=pts[3*p+2];
			for(k=0;k<nu[p];k++) ed[ed[p][k]][ed[p][nu[p]+k]]=i;
			ed[i]=ed[p];
			nu[i]=nu[p];
			ed[i][2*nu[i]]=i;
		}

		// Collapse any order 1 vertices if they were created
		if(!collapseorder1()) return false;
	}
	return true;
};

// Order 1 vertices can potentially be created during the order 2 collapse
// routine. This routine removes them.
inline bool voronoicell::collapseorder1() {
	int i,j,k;
	while(mec[1]>0) {
		cerr << "Order one collapse" << endl;
		i=--mec[1];
		j=mep[1][3*i];k=mep[1][3*i+1];
		i=mep[1][3*i+2];
		if(!delete_connection(j,k)) return false;
		--p;
		if(p!=i) {
			pts[3*i]=pts[3*p];
			pts[3*i+1]=pts[3*p+1];
			pts[3*i+2]=pts[3*p+2];
			for(k=0;k<nu[p];k++) ed[ed[p][k]][ed[p][nu[p]+k]]=i;
			ed[i]=ed[p];
			nu[i]=nu[p];
			ed[i][2*nu[i]]=i;
		}
	}
	return true;
};

// This routine deletes the kth edge of vertex j and reorganizes the memory
inline bool voronoicell::delete_connection(int j,int k) {
	int i=nu[j]-1,l,*edp,*edd,m;
	if(i<1) {
		cout << "Zero order vertex formed" << endl;
		return false;
	}
	if(mec[i]==mem[i]) addmemory(i);
	edp=mep[i]+(2*i+1)*mec[i]++;
	edp[2*i]=j;
	for(l=0;l<k;l++) {
		edp[l]=ed[j][l];
		edp[l+i]=ed[j][l+nu[j]];
	}
	while(l<i) {
		m=ed[j][l+1];
		edp[l]=m;
		k=ed[j][l+nu[j]+1];
		edp[l+i]=k;
		ed[m][nu[m]+k]--;
		l++;
	}
	edd=mep[nu[j]]+(2*nu[j]+1)*--mec[nu[j]];
	for(l=0;l<=2*nu[j];l++) ed[j][l]=edd[l];
	ed[edd[2*nu[j]]]=edd;
	ed[j]=edp;
	nu[j]=i;
	return true;
};

// Cuts a Voronoi cell using the influence of a particle at (x,y,z), first
// calculating the modulus squared of this vector before passing it to the
// routine above
inline bool voronoicell::plane(double x,double y,double z) {
	double rsq=x*x+y*y+z*z;
	return plane(x,y,z,rsq);
};

// Simple functions for moving around the edges of a given Voronoi vertex
inline int voronoicell::vor_up(int a,int p) {
	return a==nu[p]-1?0:a+1;
};
inline int voronoicell::vor_down(int a,int p) {
	return a==0?nu[p]-1:a-1;
};

// Calculates the volume of a Voronoi cell
inline double voronoicell::volume() {
	const double fe=1/48.0;
	// TODO - This memory allocation needs to be dynamic
	static unsigned int b[maxvertices];double vol=0;
	int i,j,k,l,m,n;
	double ux,uy,uz,vx,vy,vz,wx,wy,wz;
	for(i=0;i<p;i++) b[i]=0;
	for(i=1;i<p;i++) {

		// TODO - this routine will fail if there are vertices
		// whose order is bigger than the number of bits in an integer.
		if (nu[i]>30) throw overflow("Volume routine doesn't support really high order vertices");
		
		ux=pts[0]-pts[3*i];
		uy=pts[1]-pts[3*i+1];
		uz=pts[2]-pts[3*i+2];
		for(j=0;j<nu[i];j++) {
			if ((b[i]&(1<<j))==0) {
				b[i]|=1<<j;
				k=ed[i][j];
				l=vor_up(ed[i][nu[i]+j],k);
				vx=pts[3*k]-pts[0];
				vy=pts[3*k+1]-pts[1];
				vz=pts[3*k+2]-pts[2];
				b[k]|=1<<l;
				m=ed[k][l];
				while(m!=i) {
					n=vor_up(ed[k][nu[k]+l],m);
					wx=pts[3*m]-pts[0];
					wy=pts[3*m+1]-pts[1];
					wz=pts[3*m+2]-pts[2];
					vol+=ux*vy*wz+uy*vz*wx+uz*vx*wy-uz*vy*wx-uy*vx*wz-ux*vz*wy;
					b[m]|=1<<n;
					k=m;l=n;vx=wx;vy=wy;vz=wz;
					m=ed[k][l];
				}
			}
		}
	}
	return vol*fe;
};

// Computes the maximum radius squared
inline double voronoicell::maxradsq() {
	int i;double r,s;
	r=pts[0]*pts[0]+pts[1]*pts[1]+pts[2]*pts[2];
	for(i=3;i<3*p;i+=3) {
		s=pts[i]*pts[i]+pts[i+1]*pts[i+1]+pts[i+2]*pts[i+2];
		if(s>r) r=s;
	}
	return r;
};

// Outputs the edges of the Voronoi cell (in POV-Ray format) to an open file
// stream, displacing the cell by an amount (x,y,z)
inline void voronoicell::dumppov(ofstream &of,double x,double y,double z) {
	int i,j,k;double ux,uy,uz;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];uz=z+0.5*pts[3*i+2];
		of << "sphere{<" << ux << "," << uy << "," << uz << ">,r}" << endl;
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k<i) of << "cylinder{<" << ux << "," << uy << "," << uz << ">,<" << x+0.5*pts[3*k] << "," << y+0.5*pts[3*k+1] << "," << z+0.5*pts[3*k+2] << ">,r}" << endl;
		}
	}
};

// Outputs the edges of the Voronoi cell (in gnuplot format) to an open file
// stream, displacing the cell by an amount (x,y,z)
inline void voronoicell::dumpgnuplot(ofstream &of,double x,double y,double z) {
	int i,j,k;double ux,uy,uz;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];uz=z+0.5*pts[3*i+2];
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (ed[i][j]<i) of << ux << " " << uy << " " << uz << endl << x+0.5*pts[3*k] << " " << y+0.5*pts[3*k+1] << " " << z+0.5*pts[3*k+2] << endl << endl << endl;
		}
	}
};

// Randomly perturbs the points in the Voronoi cell by an amount r
inline bool voronoicell::perturb(double r) {
	for(int i=0;i<3*p;i++) {
		pts[i]+=(2*double(rand())/RAND_MAX-1)*r;
	}
};

//Initialises the suretest class and creates a buffer for dubious points
suretest::suretest() : currentdubious(initdubious) {
	sn=new int[2*currentdubious];
};

// Sets up the suretest class with a particular test plane, and removes
// any special cases from the table
inline void suretest::init(double x,double y,double z,double rsq) {
	sc=0;px=x;py=y;pz=z;prsq=rsq;
};

inline int suretest::test(int n,double &ans) {
	ans=px*p[3*n]+py*p[3*n+1]+pz*p[3*n+2]-prsq;
	if(ans>tolerance2) {
		return 1;
	} else if(ans<-tolerance2) {
		return -1;
	} else {
		int i;
		for(i=0;i<sc;i+=2) if(sn[i]==n) return sn[i+1];
		if (sc==2*currentdubious) {
			i=2*currentdubious;
			if (i>maxdubious) throw overflow("Dubious case buffer allocation exceeded absolute maximum");
			cerr << "Dubious cases buffer scaled up to " << i << endl;
			int *psn=new int[2*i];
			for(int j=0;j<2*currentdubious;j++) psn[j]=sn[j];
			delete sn;sn=psn;
		}
		sn[sc++]=n;
		sn[sc++]=ans>tolerance?1:(ans<-tolerance?-1:0);
		return sn[sc-1];
	}
};

// Prints the vertices, their edges, the relation table,
// and also notifies if any glaring memory errors are visible.
void voronoicell::edgeprint() {
	int j;
	for(int i=0;i<p;i++) {
		cout << i << " " << nu[i] << "  ";
		for(j=0;j<nu[i];j++) cout << " " << ed[i][j];
		cout << "    ";
		while(j<2*nu[i]) cout << " " << ed[i][j++];
		cout << "     " << ed[i][j];
		cout << " " << pts[3*i] << " " << pts[3*i+1] << " " << pts[3*i+2];
		cout << " " << ed[i];
		if (ed[i]>=mep[nu[i]]+mec[nu[i]]*(2*nu[i]+1)) cout << " Memory error";
		cout << endl;
	}
};
