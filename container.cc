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

const double tolerance=1e-7;
const double tolerance2=2e-7;// was -4

#include "container.hh"

// Container constructor. The first six arguments set the corners of the box to
// be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
// grid of blocks, set by the following three arguments. The next three
// arguments are booleans, which set the periodicity in each direction. The
// final argument sets the amount of memory allocated to each block.
container::container(double xa,double xb,double ya,double yb,double za,double zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi)
	: ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),zsp(zn/(zb-za)),
	nx(xn),ny(yn),nz(zn),nxy(xn*yn),nxyz(xn*yn*zn),mem(memi),
	xperiodic(xper), yperiodic(yper), zperiodic(zper) {
	co=new int[nxyz];
	for(int l=0;l<nxyz;l++) co[l]=0;
	id=new int*[nxyz];
	for(int l=0;l<nxyz;l++) id[l]=new int[mem];
	p=new double*[nxyz];
	for(int l=0;l<nxyz;l++) p[l]=new double[3*mem];
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
}

// Put a particle into the correct region of the container
void container::put(int n,double x,double y,double z) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
#ifdef OVERFLOW_CHECKING
			if(co[i]==mem) throw overflow("Not enough space in the grid");
#endif
			p[i][3*co[i]]=x;
			p[i][3*co[i]+1]=y;
			p[i][3*co[i]+2]=z;
			id[i][co[i]++]=n;
		}
	}
}

// Import a list of particles from standard input
void container::import() {
	int n;double x,y,z;
	cin >> n >> x >> y >> z;
	while(!cin.eof()) {
		put(n,x,y,z);
		cin >> n >> x >> y >> z;
	}
}

// Outputs the number of particles within each region
void container::regioncount() {
	int i,j,k,ijk=0;
	for(k=0;k<nz;k++) {
		for(j=0;j<ny;j++) {
			for(i=0;i<nx;i++) cout << "Region (" << i << "," << j << "," << k << "): " << co[ijk++] << " particles" << endl;
		}
	}
}

// Clears a container of particles
void container::clear() {
	for(int ijk=0;ijk<nxyz;ijk++) co[ijk]=0;
}

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
							if ((j!=i||s!=t)&&lrs<rs&&rs<urs) c.plane(x1,y1,z1,rs);
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
}

// If only a filename and an output type are supplied to vdraw, then assume
// that we are calculating the entire simulation region
void container::vdraw(char *filename,out_type ot) {
	vdraw(filename,ax,bx,ay,by,az,bz,ot);
}

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
						if ((j!=i||s!=t)&&lrs<rs&&rs<urs) c.plane(x1,y1,z1,rs);
					}
				} while ((t=l.inc(qx,qy,qz))!=-1);
				lr=ur;lrs=urs;
			}
			bb[id[s][i]]=c.volume();
		}
	}
}

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
						if ((j!=i||s!=t)&&lrs<rs&&rs<urs) c.plane(x1,y1,z1,rs);
					}
				} while ((t=l.inc(qx,qy,qz))!=-1);
				lr=ur;lrs=urs;
			}
			cout << id[s][i] << " " << x << " " << y << " " << z << " " << c.volume() << endl;
		}
	}
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
}

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
}

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
}

// Custom int function, that gives consistent stepping for negative numbers.
// With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
// With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).
template <class T>
inline int loop::myint(T a) {
	return a<0?int(a)-1:int(a);
}

// Custom mod function, that gives consistent stepping for negative numbers
inline int loop::mymod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

// Custom div function, that gives consistent stepping for negative numbers
inline int loop::mydiv(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}

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
}

// Increases the memory storage for a particular vertex order
void voronoicell::addmemory(int i) {
	int s=2*i+1;
	cout << "Extramem " << i << endl;
	if(mem[i]==0) {
		mep[i]=new int[initnvertices*s];
		mem[i]=initnvertices;
	} else {
		int j,k,*l;
		mem[i]*=2;
		l=new int[s*mem[i]];
		j=0;
		while(j<s*mec[i]) {
			ed[mep[i][j+2*i]]=l+j;
			for(k=0;k<=2*i;k++) {
				l[j]=mep[i][j];
				j++;
			}
		}
		delete mep[i];
		mep[i]=l;
	}
}

void voronoicell::addmemory_vertices() {
	int i=2*currentvertices,j,**ped,*pnu;
	cout << "Vertex memory scaled up to " << i << endl;
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
}

void voronoicell::addmemory_vorder() {
	int i=2*currentvertexorder,j,*pmem,**pmep,*pmec;
	cout << "Vertex order memory scaled up to " << i << endl;
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
}

void voronoicell::addmemory_ds() {
	int i=2*currentdeletesize,j,*pds;
	cout << "Delete stack 1 memory scaled up to " << i << endl;
	pds=new int[i];
	for(j=0;j<currentdeletesize;j++) pds[j]=ds[j];
	delete ds;ds=pds;
	currentdeletesize=i;
}

void voronoicell::addmemory_ds2() {
	int i=2*currentdeletesize2,j,*pds2;
	cout << "Delete stack 2 memory scaled up to " << i << endl;
	pds2=new int[i];
	for(j=0;j<currentdeletesize2;j++) pds2[j]=ds2[j];
	delete ds2;ds2=pds2;
	currentdeletesize2=i;
}

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

inline void voronoicell::init_test() {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;p=0;

	/*
	add_vertex(1,-2,-1,5,1,3);
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

/*	add_vertex(1,-2,-1,4,3,1);
	add_vertex(-1,-2,-1,5,4,0,2);
	add_vertex(-1,2,-1,3,6,1);
	add_vertex(1,2,-1,0,5,6,2);
	add_vertex(0,-1,1,0,1,5);
	add_vertex(0,0,0,6,3,4,1);
	add_vertex(0,1,1,3,5,2);*/

/*	add_vertex(-1,-1,-1,1,3,4);
	add_vertex(1,-1,-1,5,2,0);
	add_vertex(1,1,-1,3,1,6);
	add_vertex(-1,1,-1,7,0,2);
	add_vertex(-1,-1,1,8,5,0,7);
	add_vertex(1,-1,1,8,6,1,4);
	add_vertex(1,1,1,8,7,2,5);
	add_vertex(-1,1,1,8,4,3,6);
	add_vertex(0,0,2,6,5,4,7);*/

/*	add_vertex(1,-3,-1,1,6,5);
	add_vertex(-1,-3,-1,2,6,0);
	add_vertex(-3,0,-1,3,8,7,1);
	add_vertex(-1,3,-1,4,9,2);
	add_vertex(1,3,-1,5,9,3);
	add_vertex(3,0,-1,0,7,8,4);
	add_vertex(0,-2,1,0,1,7);
	add_vertex(0,-1,0,5,6,2,8);
	add_vertex(0,1,0,5,7,2,9);
	add_vertex(0,2,1,4,8,3);*/

/*	add_vertex(-1,-3,-1,12,8,1,7);
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

/*	add_vertex(2,-3,-1,1,4,3);
	add_vertex(-2,-3,-1,2,4,0);
	add_vertex(-2,3,-1,3,7,1);
	add_vertex(2,3,-1,0,6,2);
	add_vertex(0,-2,0,1,5,0);
	add_vertex(0,1,0,7,6,4);
	add_vertex(1,2,1,3,5,7);
	add_vertex(-1,2,1,2,6,5);*/

/*	add_vertex(3,-2,-1,1,3,2);
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
}

void voronoicell::add_vertex(double x,double y,double z,int a) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=1;
	if (mem[1]=mec[1]) addmemory(1);
	int *q=mep[1]+3*mec[1]++;ed[p]=q;
	q[0]=a;q[2]=p++;
}

void voronoicell::add_vertex(double x,double y,double z,int a,int b) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=2;
	if (mem[2]=mec[2]) addmemory(2);
	int *q=mep[2]+5*mec[2]++;ed[p]=q;
	q[0]=a;q[1]=b;q[4]=p++;
}

void voronoicell::add_vertex(double x,double y,double z,int a,int b,int c) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=3;
	if (mem[3]=mec[3]) addmemory(3);
	int *q=mep[3]+7*mec[3]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[6]=p++;
}

void voronoicell::add_vertex(double x,double y,double z,int a,int b,int c,int d) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=4;
	if (mem[4]=mec[4]) addmemory(4);
	int *q=mep[4]+9*mec[4]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[8]=p++;
}

void voronoicell::add_vertex(double x,double y,double z,int a,int b,int c,int d,int e) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=5;
	if (mem[5]=mec[5]) addmemory(5);
	int *q=mep[5]+11*mec[5]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[4]=e;q[10]=p++;
}

// Checks that the relational table of the Voronoi cell is accurate, and prints
// out any errors
inline void voronoicell::relcheck() {
	int i,j;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if (ed[ed[i][j]][ed[i][nu[i]+j]]!=i) cout << "Relational error at point " << i << ", edge " << j << "." << endl;
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
inline bool voronoicell::plane(double x,double y,double z,double rsq) {
	int count=0,i,j,k,up=6,lp=6,tp,cp,qp,rp,stack=0,stack2=0;
	int us,ls,qs,iqs,cs,uw,qw,lw,tw;
	int *edp,*emp;
	double u,l,t,r,q;bool complicatedsetup=false,newdoubleedge=false,doubleedge=false;
	sure.init(x,y,z,rsq);
	uw=sure.test(up,u);l=u;
	try {
		if(uw==1) {
			do {
				if (++count>=p) throw true;
				u=l;up=lp;uw=lw;
					for(i=0;i<nu[up];i++) {
					tp=ed[up][i];
					tw=sure.test(tp,t);
					if(t<l) {l=t;lw=tw;lp=tp;us=i;}
				}
				if (lp==up) return false;  // Cell no longer exists
			} while (lw==1);
			ls=ed[up][nu[up]+us];
			if (lw==-1) {
				complicatedsetup=false;
			} else {
				// lp is in the plane, up is inside
				up=lp;
				complicatedsetup=true;
			}
		} else if (uw==-1) {
			do {
				if (++count>=p) throw true;
				l=u;lp=up;lw=uw;
				for(i=0;i<nu[lp];i++) {
					tp=ed[lp][i];
					tw=sure.test(tp,t);
					if(t>u) {u=t;uw=tw;up=tp;ls=i;}
				}
				if (up==lp) return true;  // Cell isn't intersected at all
			} while (uw==-1);
			us=ed[lp][nu[lp]+ls];
			complicatedsetup=(uw!=1);
		} else {
			complicatedsetup=true;
		}
	}
	catch(bool except) {
		cout << "Bailed out of convex calculation\n";
		for(qp=0;qp<p;qp++) {
			qw=sure.test(qp,q);
			if (qw==1) {
				for(i=0;i<nu[qp];i++) {
					lp=ed[qp][i];
					if(lp<qp) {
						lw=sure.test(lp,u);
						if (lw!=1) break;
					}
				}
				if(lw!=1) {
					up=qp;
					if(lw==0) {
						complicatedsetup=true;
					} else {
						complicatedsetup=false;
						u=q;
					}
					break;
				}
			} else if (qw==-1) {
				for(i=0;i<nu[qp];i++) {
					up=ed[qp][i];
					if(up<qp) {
						uw=sure.test(up,u);
						if (uw!=-1) break;
					}
				}
				if(uw!=-1) {
					lp=qp;
					if(uw==0) {
						complicatedsetup=true;
					} else {
						complicatedsetup=false;
						u=q;
					}
					break;
				}
			} else {
				up=qp;
				complicatedsetup=true;
				break;
			}
		}
	}
	if(p==currentvertices) addmemory_vertices();
	if (complicatedsetup) {
		pts[3*p]=pts[3*up];
		pts[3*p+1]=pts[3*up+1];
		pts[3*p+2]=pts[3*up+2];
		i=0;
		lp=ed[up][0];
		lw=sure.test(lp,l);
		if(lw!=-1) {
			cout << "case 1\n";
			rp=lw;
			do {
				i++;
				if (i==nu[up]) return false;  // Marginal vertex with all neighbors inside/marginal
				lp=ed[up][i];
				lw=sure.test(lp,l);
			} while (lw!=-1);
			j=i+1;
			while(j<nu[up]) {
				lp=ed[up][j];
				lw=sure.test(lp,l);
				if (lw!=-1) break;
				j++;
			}
			cout << i << " " << j << " " << rp << endl;
			if(j==nu[up]&&i==1&&rp==0) {
				nu[p]=nu[up];
				doubleedge=true;
			} else {
				nu[p]=j-i+2;
			}
			k=1;
			if (nu[p]==currentvertexorder) addmemory_vorder();
			if (mec[nu[p]]==mem[nu[p]]) addmemory(nu[p]);
			ed[p]=mep[nu[p]]+(2*nu[p]+1)*mec[nu[p]]++;
			ed[p][2*nu[p]]=p;
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
			cs=k;
		//	if (k!=nu[p]-1) throw overflow("Doesn't add up");
			if (stack2==currentdeletesize2) addmemory_ds2();
			ds2[stack2++]=up;

			qp=up;q=u;
			i=ed[up][us];
			us=ed[up][nu[up]+us];
			up=i;
			uw=sure.test(up,u);
			if(uw==0) ed[up][2*nu[up]]=-1;
			ed[qp][2*nu[qp]]=-p;
		} else {
			cout << "case 2\n";
			i=nu[up]-1;
			lp=ed[up][i];
			lw=sure.test(lp,l);
			while(lw==-1) {
				i--;
				if (i==0) return true; // Marginal vertex with all neighbors outside
				lp=ed[up][i];
				lw=sure.test(lp,l);
			}
			j=1;
			qp=ed[up][j];
			qw=sure.test(qp,q);
			while(qw==-1) {
				j++;
				qp=ed[up][j];
				qw=sure.test(qp,l);
			}
			if (i==j&&qw==0) {
				doubleedge=true;
				nu[p]=nu[up];
			} else {
				nu[p]=nu[up]-i+j+1;
			}
			k=1;
			if(nu[p]==currentvertexorder) addmemory_vorder();
			if (mec[nu[p]]==mem[nu[p]]) addmemory(nu[p]);
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
			cs=k;
		//	if (k!=nu[p]-1) throw overflow("Doesn't add up");
			if (stack2==currentdeletesize2) addmemory_ds2();
			ds2[stack2++]=up;

			qp=up;q=u;
			i=ed[up][us];
			us=ed[up][nu[up]+us];
			up=i;
			uw=sure.test(up,u);
			if(uw==0) ed[up][2*nu[up]]=-1;
			ed[qp][2*nu[qp]]=-p;
		}
	} else {
		if (stack==currentdeletesize) addmemory_ds();
		ds[stack++]=up;
		r=1/(u-l);
		pts[3*p]=(pts[3*lp]*u-pts[3*up]*l)*r;
		pts[3*p+1]=(pts[3*lp+1]*u-pts[3*up+1]*l)*r;
		pts[3*p+2]=(pts[3*lp+2]*u-pts[3*up+2]*l)*r;
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
		qs=vor_up(us,up);
		qp=up;q=u;
	}

	// What do we want by this point?
	// Set up first point in facet
	// Either create it, or modify existing vertex
	// Have a direction to move in
	// Know whether that direction is inside, or boundary
	// Know enough stuff to join back the end
	cp=p;rp=p;p++;
	while(qp!=up||qs!=us) {
		lp=ed[qp][qs];
		lw=sure.test(lp,l);
		if (lw==1) {
			if (stack==currentdeletesize) addmemory_ds();
			qs=vor_up(ed[qp][nu[qp]+qs],lp);qp=lp;q=l;ds[stack++]=qp;
		} else if (lw==-1) {
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
			if(p==currentvertices) addmemory_vertices();
			k=doubleedge?0:1;
			qs=ed[qp][nu[qp]+qs];
			qp=lp;iqs=qs;
			do {
				k++;
				qs=vor_up(qs,qp);
				lp=ed[qp][qs];
				lw=sure.test(lp,l);
			} while (lw==-1);
			j=-ed[qp][2*nu[qp]];
			if(j>0) {
				if (j!=1) {
					k+=nu[j];
					if(lw==0) {
						i=-ed[lp][2*nu[lp]];
						if(i>1) {
							if(ed[i][nu[i]-1]==j) {
								// Then we need to do one less edge
								newdoubleedge=true;
								k-=1;
							} else newdoubleedge=false;
						} else {
							if (i==1&&j==rp) {
							       k-=1;
						       	       newdoubleedge=true;
							} else newdoubleedge=false;
						}
					} else newdoubleedge=false;
				} else {
					newdoubleedge=false;
					j=0;
				}
			} else newdoubleedge=false;
			
			if(k==currentvertexorder) addmemory_vorder();
			if (mec[k]==mem[k]) addmemory(k);
			if(j>0) {
				edp=mep[k]+(2*k+1)*mec[k]++;
				i=0;
				while(i<nu[j]) {
					edp[i]=ed[j][i];
					edp[k+i]=ed[j][nu[j]+i];
					i++;
				}
				edp[2*k]=j;
				mec[nu[j]]--;
				emp=mep[nu[j]]+(2*nu[j]+1)*mec[nu[j]];
				if(emp!=ed[j]) {
					for(lw=0;lw<=2*nu[j];lw++) ed[j][lw]=emp[lw];
					ed[emp[2*nu[j]]]=ed[j];
				}
				ed[j]=edp;
			} else {
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

			if (!doubleedge) {
				ed[j][i]=cp;
				ed[j][nu[j]+i]=cs;
				ed[cp][cs]=j;
				ed[cp][nu[cp]+cs]=i;
				i++;
			}

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
			doubleedge=newdoubleedge;
		} 
	}

	// Connect the final created vertex to the initial one
	ed[cp][cs]=rp;
	ed[rp][0]=cp;
	ed[cp][nu[cp]+cs]=0;
	ed[rp][nu[rp]+0]=cs;

	// Delete points
	// First remove any duplicates
	i=0;
	while(i<stack) {
		j=ds[i];
		if(ed[j][nu[j]]!=-1) {
			ed[j][nu[j]]=-1;
			i++;
		} else ds[i]=ds[--stack];
	}
	
	// Add in the degenerate ones, and reset their back pointers
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
			ed[qp]=ed[p];nu[qp]=nu[p];
			for(i=0;i<nu[qp];i++) {
				ed[ed[qp][i]][ed[qp][nu[qp]+i]]=qp;
			}
			ed[qp][2*nu[qp]]=qp;
		} else p++;
	}

	// Check for any vertices with zero or one connection
	if (mec[0]>0) throw overflow("Zero order vertex formed");
	
	return collapseorder2();
}

// Collapse any order one vertices
inline bool voronoicell::collapseorder1() {
	int i,j,k;
	while(mec[1]>0) {
		cout << "Order one collapse" << endl;
		i=--mec[1];
		j=mep[1][3*i];k=mep[1][3*i+1];
		i=mep[1][3*i+2];
		delete_connection(j,k);
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
}

inline bool voronoicell::collapseorder2() {
	if(!collapseorder1()) return false;
	int a,b,i,j,k,l;
	while(mec[2]>0) {
		cout << "Order two collapse" << endl;
		i=--mec[2];
		j=mep[2][5*i];k=mep[2][5*i+1];
		if (j==k) return false;
		for(l=0;l<nu[j];l++) {
			if(ed[j][l]==k) break;
		}
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
		if(!collapseorder1()) return false;
	}
	return true;
}

inline bool voronoicell::delete_connection(int j,int k) {
	int i=nu[j]-1,l,*edp,*edd,m;
	if(i<1) return false;
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
}


// Cuts a Voronoi cell using the influence of a particle at (x,y,z), first
// calculating the modulus squared of this vector before passing it to the
// routine above
inline bool voronoicell::plane(double x,double y,double z) {
	double rsq=x*x+y*y+z*z;
	return plane(x,y,z,rsq);
}

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
	static unsigned int b[maxvertices];double vol=0;
	int i,j,k,l,m,n;
	double ux,uy,uz,vx,vy,vz,wx,wy,wz;
	for(i=0;i<p;i++) b[i]=0;
	for(i=1;i<p;i++) {
		for(j=0;j<3;j++) {
			ux=pts[0]-pts[i];
			uy=pts[1]-pts[i+1];
			uz=pts[2]-pts[i+2];
			if (b[i]&(1<<j)==0) {
				b[i]|=1<<j;
				k=ed[i][j];
				l=vor_up(ed[i][nu[i]+j],i);
				vx=pts[k]-pts[0];
				vy=pts[k+1]-pts[1];
				vz=pts[k+2]-pts[2];
				b[k]|=1<<l;
				m=ed[k][l];
				while(m!=i) {
					n=vor_up(ed[k][nu[k]+l],k);
					wx=pts[m]-pts[0];
					wy=pts[m+1]-pts[1];
					wz=pts[m+2]-pts[2];
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
	for(i=3;i<p;i+=3) {
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

// Sets up the suretest class with a particular test plane, and removes
// any special cases from the table
inline void suretest::init(double x,double y,double z,double rsq) {
	currentdubious=initdubious;
	sn=new int[2*currentdubious];
	sc=0;px=x;py=y;pz=z;prsq=rsq;
};

/* 
inline void suretest::test(int n,double &ans) {
	ans=px*p[3*n]+py*p[3*n+1]+pz*p[3*n+2]-prsq;
	if(ans>tolerance2||ans<-tolerance2) {
		return;
	} else {
		for(int i=0;i<sc;i+=2) if(sn[i]==n) return;
		sn[sc]=ans>tolerance?1:(ans<-tolerance?-1:0);
		sc+=2;
	}
};*/

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
			cout << "Dubious cases buffer scaled up to " << i << endl;
			int *psn=new int[2*i];
			for(int j=0;j<2*currentdubious;j++) psn[j]=sn[j];
			delete sn;sn=psn;
		}
		sn[sc++]=n;
		sn[sc++]=ans>tolerance?1:(ans<-tolerance?-1:0);
		return sn[sc-1];
	}
};

void voronoicell::edgeprint(bool extend) {
	int j;
	for(int i=0;i<p;i++) {
		cout << i << " " << nu[i] << "  ";
		for(j=0;j<nu[i];j++) cout << " " << ed[i][j];
		if (extend) {
			cout << "  ";
			while(j<2*nu[i]) cout << " " << ed[i][j++];
			cout << "   " << ed[i][j];
			cout << " " << pts[3*i] << " " << pts[3*i+1] << " " << pts[3*i+2];
		}
		cout << " " << ed[i];
		if (ed[i]>=mep[nu[i]]+mec[nu[i]]*(2*nu[i]+1)) cout << " crazzzy";
		cout << endl;
	}
};
