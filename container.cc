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
const int maxvertices=400;
const int maxvertexorder=20;
const int initial3vertices=256;
const int initialnvertices=16;
const int maxdubiouscases=64;
const double tolerance=1e-8;
const double tolerance2=2e-8;

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
voronoicell::voronoicell() : sure(pts) {
	mem[0]=initial3vertices;
	mep[0]=new int[initial3vertices*7];
	mec[0]=0;
	for(int i=1;i<maxvertexorder-3;i++) {
		mem[i]=initialnvertices;
		mep[i]=new int[initialnvertices*(2*i+1)];
		mec[i]=0;
	}
}

// Increases the memory storage for a particular vertex order
void voronoicell::addmemory(int i) {
	cout << "Extramem " << i << endl;
	if(mem[i]==0) {
		mep[i]=new int[initialnvertices*(2*i+1)];
		mem[i]=initialnvertices;
	} else {
		int j=0,k;int *l;
		mem[i]*=2;
		l=new int[mem[i]];
		while(j<(2*i+1)*mec[i]) {
			for(k=0;k<2*i;k++) l[j++]=mep[i][j];
			ed[mep[i][j]]=l+j;
			l[j++]=mep[i][j];
		}
		delete mep[i];
		mep[i]=l;
	}
}

// Initializes a Voronoi cell as a rectangular box with the given dimensions
inline void voronoicell::init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
	mec[0]=p=8;xmin*=2;xmax*=2;ymin*=2;ymax*=2;zmin*=2;zmax*=2;
	pts[0]=xmin;pts[1]=ymin;pts[2]=zmin;
	pts[3]=xmax;pts[4]=ymin;pts[5]=zmin;
	pts[6]=xmin;pts[7]=ymax;pts[8]=zmin;
	pts[9]=xmax;pts[10]=ymax;pts[11]=zmin;
	pts[12]=xmin;pts[13]=ymin;pts[14]=zmax;
	pts[15]=xmax;pts[16]=ymin;pts[17]=zmax;
	pts[18]=xmin;pts[19]=ymax;pts[20]=zmax;
	pts[21]=xmax;pts[22]=ymax;pts[23]=zmax;
	int *q=mep[0];
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
	for(int i=1;i<maxvertexorder-3;i++) mec[i]=0;
};

// Initializes a Voroni cell as a regular octahedron
inline void voronoicell::init_octahedron(double l) {
	mec[0]=0;mec[1]=p=6;l*=2;
	pts[0]=-l;pts[1]=0;pts[2]=0;
	pts[3]=l;pts[4]=0;pts[5]=0;
	pts[6]=0;pts[7]=-l;pts[8]=0;
	pts[9]=0;pts[10]=l;pts[11]=0;
	pts[12]=0;pts[13]=0;pts[14]=-l;
	pts[15]=0;pts[16]=0;pts[17]=l;
	int *q=mep[1];
	q[0]=2;q[1]=5;q[2]=3;q[3]=4;q[4]=0;q[5]=0;q[6]=0;q[7]=0;q[8]=0;
	q[9]=2;q[10]=4;q[11]=3;q[12]=5;q[13]=2;q[14]=2;q[15]=2;q[16]=2;q[17]=1;
	q[18]=0;q[19]=4;q[20]=1;q[21]=5;q[22]=0;q[23]=3;q[24]=0;q[25]=1;q[26]=2;
	q[27]=0;q[28]=5;q[29]=1;q[30]=4;q[31]=2;q[32]=3;q[33]=2;q[34]=1;q[35]=3;
	q[36]=0;q[37]=3;q[38]=1;q[39]=2;q[40]=3;q[41]=3;q[42]=1;q[43]=1;q[44]=4;
	q[45]=0;q[46]=2;q[47]=1;q[48]=3;q[49]=1;q[50]=3;q[51]=3;q[52]=1;q[53]=5;
	ed[0]=q;ed[1]=q+9;ed[2]=q+18;ed[3]=q+27;ed[4]=q+36;ed[5]=q+45;
	nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=4;
	for(int i=2;i<maxvertexorder-3;i++) mec[i]=0;
};

// Checks that the relational table of the Voronoi cell is accurate, and prints
// out any errors
inline void voronoicell::relcheck() {
	int i,j,ij=0;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if (ed[ed[i][j]][ed[i][nu[i]+j]]!=i) cout << "Relational error at point " << i << ", edge " << j << "." << endl;
		}
	}
};

// Cuts the Voronoi cell by a particle whose center is at a separation of
// (x,y,z) from the cell center. The value of rsq should be initially set to
// x*x+y*y+z*z.
inline bool voronoicell::plane(double x,double y,double z,double rsq) {
	int i,j,up=0,lp=0,tp,cp,qp,rp,stack;
	int us,ls,ts,qs,uw,lw,tw;
	double u,l,t,r,q;
	cout << "Start" << endl;
	edgeprint(true);cout << endl;
	sure.init(x,y,z,rsq);
	uw=sure.test(up,u);l=u;
	static int ds[3*maxvertices];
	if(uw==1) {
		do {
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
			cout << "top\n";
			// lp is outside, up is inside
			// Proceed regular iteration by making new point
			ds[stack++]=up;
			r=1/(u-l);
			pts[3*p]=(pts[3*lp]*u-pts[3*up]*l)*r;
			pts[3*p+1]=(pts[3*lp+1]*u-pts[3*up+1]*l)*r;
			pts[3*p+2]=(pts[3*lp+2]*u-pts[3*up+2]*l)*r;
			nu[p]=3;
			if (mec[0]==mem[0]) addmemory(0);
			ed[p]=mep[0]+7*mec[0]++;
			ed[p][6]=p;
			ed[up][us]=-1;
			ed[lp][ls]=p;
			ed[lp][nu[lp]+ls]=0;
			ed[p][0]=lp;
		} else {
			// lp is in the plane, up is inside
			throw overflow("Not supported yet");
		}
	} else if (uw==-1) {
		do {
			l=u;lp=up;lw=uw;
			for(i=0;i<nu[lp];i++) {
				tp=ed[lp][i];
				tw=sure.test(tp,t);
				if(t>u) {u=t;uw=tw;up=tp;ls=i;}
			}
			if (up==lp) return true;  // Cell isn't intersected at all
		} while (uw==-1);
		us=ed[lp][nu[lp]+ls];cout << up << " " << lp << " " << p << endl;
		if (uw==1) {
			// lp is outside, up is inside
			// Proceed regular iteration
			ds[stack++]=up;
			r=1/(u-l);
			pts[3*p]=(pts[3*lp]*u-pts[3*up]*l)*r;
			pts[3*p+1]=(pts[3*lp+1]*u-pts[3*up+1]*l)*r;
			pts[3*p+2]=(pts[3*lp+2]*u-pts[3*up+2]*l)*r;
			nu[p]=3;
			if (mec[0]==mem[0]) addmemory(0);
			ed[p]=mep[0]+7*mec[0]++;
			ed[p][6]=p;
			ed[up][us]=-1;
			ed[lp][ls]=p;
			ed[lp][nu[lp]+ls]=0;
			ed[p][0]=lp;
		} else {
			throw overflow("Not supported yet");
			// up is in the plane, lp is outside
		}
	} else {
		// The point is in the cutting plane and that's all we know
		// Scan connections
		// A points inside, B points outside, C on boundary
		// C==2, B==0, then bail out with false
		// Put all As on delete stack
		throw overflow("Not supported yet");
	}
	// What do we want by this point?
	// Set up first point in facet
	// Either create it, or modify existing vertex
	// Have a direction to move in
	// Know whether that direction is inside, or boundary
	// Know enough stuff to join back the end
	edgeprint(true);
	cout << "smd" << endl;
	cp=p;
	rp=p;
	ed[p][nu[p]]=ls;
	qp=up;q=u;qs=vor_up(us,qp);p++;
	edgeprint(true);
	while(qp!=up||qs!=us) {
		cout << "st" << qp << " " << up << " " << qs << " " << us << endl;
		lp=ed[qp][qs];
		lw=sure.test(lp,l);
		if (lw==1) {
			qs=vor_up(ed[qp][nu[qp]+qs],qp);qp=lp;q=l;ds[stack++]=qp;
		} else if (lw==-1) {
			r=1/(q-l);
			pts[3*p]=(pts[3*lp]*q-pts[3*qp]*l)*r;
			pts[3*p+1]=(pts[3*lp+1]*q-pts[3*qp+1]*l)*r;
			pts[3*p+2]=(pts[3*lp+2]*q-pts[3*qp+2]*l)*r;
			nu[p]=3;
			if (mec[0]==mem[0]) addmemory(0);
			ed[p]=mep[0]+7*mec[0]++;
			ed[p][6]=p;
			ls=ed[qp][qs+nu[qp]];
			ed[lp][ls]=p;
			ed[lp][nu[lp]+ls]=0;
			ed[p][0]=lp;
			ed[p][2]=cp;
			ed[p][nu[p]]=ls;
			ed[p][nu[p]+2]=1;
			ed[cp][1]=p;
			ed[cp][nu[cp]+1]=2;
			ed[qp][qs]=-1;
			qs=vor_up(qs,qp);
			cp=p++;
		} else {
			throw overflow("Not supported yet");
		}
		cout << "en" << qp << " " << up << " " << qs << " " << us << endl;
	}
	ed[cp][1]=rp;
	ed[rp][2]=cp;
	ed[cp][nu[cp]+1]=2;
	ed[rp][nu[rp]+2]=1;

	cout << "Planed" << endl;
	edgeprint(true);
	
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

	// Scan connections and add in extras
	for(i=0;i<stack;i++) {
		cp=ds[i];
		for(j=0;j<nu[cp];j++) {
			qp=ed[cp][j];
			if(qp!=-1) {
				if (ed[qp][nu[qp]]!=-1) {
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
			mec[j-3]--;
			for(i=0;i<=2*j;i++) ed[p][i]=(mep[j-3]+(2*j+1)*mec[j-3])[i];
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
			mec[j-3]--;
			for(i=0;i<=2*j;i++) ed[qp][i]=(mep[j-3]+(2*j+1)*mec[j-3])[i];
			ed[ed[qp][2*j]]=ed[qp];

			// Edge management
			ed[qp]=ed[p];nu[qp]=nu[p];
			for(i=0;i<nu[qp];i++) {
				ed[ed[qp][i]][ed[qp][nu[qp]+i]]=qp;
			}
			ed[qp][2*nu[qp]]=qp;
		} else p++;
	}
	cout << endl;
	cout << "Deleted" << endl;
	edgeprint(true);
	return false;
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
		cout << "dubious\n";
		for(int i=0;i<sc;i+=2) if(sn[i]==n) return sn[i+1];
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
		cout << endl;
	}
};
