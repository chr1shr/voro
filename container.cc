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

// Initializes a Voronoi cell as a rectangular box with the given dimensions
inline void voronoicell::init(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
	p=24;xmin*=2;xmax*=2;ymin*=2;ymax*=2;zmin*=2;zmax*=2;
	pts[0]=xmin;pts[1]=ymin;pts[2]=zmin;
	pts[3]=xmax;pts[4]=ymin;pts[5]=zmin;
	pts[6]=xmin;pts[7]=ymax;pts[8]=zmin;
	pts[9]=xmax;pts[10]=ymax;pts[11]=zmin;
	pts[12]=xmin;pts[13]=ymin;pts[14]=zmax;
	pts[15]=xmax;pts[16]=ymin;pts[17]=zmax;
	pts[18]=xmin;pts[19]=ymax;pts[20]=zmax;
	pts[21]=xmax;pts[22]=ymax;pts[23]=zmax;
	ed[0]=3;ed[1]=12;ed[2]=6;ed[3]=9;ed[4]=15;ed[5]=0;
	ed[6]=0;ed[7]=18;ed[8]=9;ed[9]=6;ed[10]=21;ed[11]=3;
	ed[12]=18;ed[13]=0;ed[14]=15;ed[15]=12;ed[16]=3;ed[17]=21;
	ed[18]=21;ed[19]=6;ed[20]=12;ed[21]=15;ed[22]=9;ed[23]=18;
	rl[0]=2;rl[1]=1;rl[2]=0;rl[3]=2;rl[4]=1;rl[5]=0;
	rl[6]=2;rl[7]=1;rl[8]=0;rl[9]=2;rl[10]=1;rl[11]=0;
	rl[12]=2;rl[13]=1;rl[14]=0;rl[15]=2;rl[16]=1;rl[17]=0;
	rl[18]=2;rl[19]=1;rl[20]=0;rl[21]=2;rl[22]=1;rl[23]=0;
};

// Checks that the relational table of the Voronoi cell is accurate, and prints
// out any errors
inline void voronoicell::relcheck() {
	int i,j,ij=0;
	for(i=0;i<p;i+=3) {
		for(j=0;j<3;j++,ij++) {
			if (ed[ed[ij]+rl[ij]]!=i) cout << "Relational error at point " << i/3 << ", edge " << j << "." << endl;
		}
	}
};

// Cuts the Voronoi cell by a particle whose center is at a separation of
// (x,y,z) from the cell center. The value of rsq should be initially set to
// x*x+y*y+z*z
inline bool voronoicell::plane(double x,double y,double z,double rsq) {
	int up=0,lp=0,qp,rp,cp,us,cs,qs,rs,ls,stack=1,i,j;
	static int ds[3*maxvertices];
	double u=x*pts[0]+y*pts[1]+z*pts[2],l=u,q,r;
	if (u>rsq) {
		do {
			u=l;up=lp;
			lp=ed[up];qp=ed[up+1];rp=ed[up+2];
			l=x*pts[lp]+y*pts[lp+1]+z*pts[lp+2];
			q=x*pts[qp]+y*pts[qp+1]+z*pts[qp+2];
			r=x*pts[rp]+y*pts[rp+1]+z*pts[rp+2];
			if (r<u) {
				if (q<r) {if (l<q) us=0;else {lp=qp;l=q;us=1;}}
				else {if (l<r) us=0;else {lp=rp;l=r;us=2;}}
			} else {
				if (q<u) {if (l<q) us=0;else {lp=qp;l=q;us=1;}}
				else {if (l<u) us=0;else return true;}
			}
		} while(l>rsq);
		ls=rl[up+us];
	} else {
		do {
			lp=up;l=u;
			up=ed[lp];qp=ed[lp+1];rp=ed[lp+2];
			u=x*pts[up]+y*pts[up+1]+z*pts[up+2];
			q=x*pts[qp]+y*pts[qp+1]+z*pts[qp+2];
			r=x*pts[rp]+y*pts[rp+1]+z*pts[rp+2];
			if (r>l) {
				if (q>r) {if (u>q) ls=0;else {up=qp;u=q;ls=1;}}
				else {if (u>r) ls=0;else {up=rp;u=r;ls=2;}}
			} else {
				if (q>l) {if (u>q) ls=0;else {up=qp;u=q;ls=1;}}
				else {if (u>l) ls=0;else return false;}
			}
		} while(u<rsq);
		us=rl[lp+ls];
	}
	ds[0]=up;
	ed[up+us]=-1;
	ed[lp+ls]=p;
	ed[p+us]=lp;
	cp=p;cs=vor_up(us); 
	rp=p;rs=vor_down(us);
	rl[p+us]=ls;
	r=1/(u-l);
	pts[p++]=(pts[up]*(rsq-l)+pts[lp]*(u-rsq))*r;
	pts[p++]=(pts[up+1]*(rsq-l)+pts[lp+1]*(u-rsq))*r;
	pts[p++]=(pts[up+2]*(rsq-l)+pts[lp+2]*(u-rsq))*r;
	qp=up;q=u;qs=vor_up(us);
	while(qp!=up||qs!=us) {
		lp=ed[qp+qs];
		l=x*pts[lp]+y*pts[lp+1]+z*pts[lp+2];
		if (l>rsq) {
			qs=vor_up(rl[qp+qs]);qp=lp;q=l;ds[stack++]=qp;
		} else {
#ifdef OVERFLOW_CHECKING
			if (p==3*maxvertices) throw overflow("Not enough memory for points");
#endif
			ls=rl[qp+qs];
			ed[lp+ls]=p;
			ed[p+qs]=lp;
			ed[p+vor_down(qs)]=cp;
			rl[p+qs]=ls;
			rl[p+vor_down(qs)]=cs;
			ed[cp+cs]=p;
			rl[cp+cs]=vor_down(qs);
			cp=p;
			r=1/(q-l);
			pts[p++]=(pts[qp]*(rsq-l)+pts[lp]*(q-rsq))*r;
			pts[p++]=(pts[qp+1]*(rsq-l)+pts[lp+1]*(q-rsq))*r;
			pts[p++]=(pts[qp+2]*(rsq-l)+pts[lp+2]*(q-rsq))*r;
			ed[qp+qs]=-1;
			qs=vor_up(qs);
			cs=qs;
		}
	}
	ed[cp+cs]=rp;
	ed[rp+rs]=cp;
	rl[cp+cs]=rs;
	rl[rp+rs]=cs;
	i=0;
	while(i<stack) {
		if (rl[ds[i]]!=3) rl[ds[i++]]=3;
		else ds[i]=ds[--stack];
	}
	for(i=0;i<stack;i++) {
		cp=ds[i];
		for(j=0;j<3;j++) {
			qp=ed[cp+j];
			if(qp!=-1) {
				if (rl[qp]!=3) {ds[stack++]=qp;rl[qp]=3;}
			}
		}
	}
	while(stack>0) {
		while(rl[p-=3]==3);
		qp=ds[--stack];
		if (qp<p) {
			pts[qp]=pts[p];
			pts[qp+1]=pts[p+1];
			pts[qp+2]=pts[p+2];
			for(j=0;j<3;j++) {
				ed[ed[p+j]+rl[p+j]]=qp;
				ed[qp+j]=ed[p+j];
				rl[qp+j]=rl[p+j];
			}
		} else p+=3;
	}
	return false;
};

// Cuts a Voronoi cell using the influence of a particle at (x,y,z), first
// calculating the modulus squared of this vector before passing it to the
// routine above
inline bool voronoicell::plane(double x,double y,double z) {
	double rsq=x*x+y*y+z*z;
	return plane(x,y,z,rsq);
}

// Simple functions for moving around the edges of a given Voronoi vertex
inline int voronoicell::vor_up(int a) {
	return a==2?0:a+1;
};
inline int voronoicell::vor_down(int a) {
	return a==0?2:a-1;
};

// Calculates the volume of a Voronoi cell
inline double voronoicell::volume() {
	const double fe=1/48.0;
	static bool b[3*maxvertices];double vol=0;
	int i,j,k,l,m,n;
	double ux,uy,uz,vx,vy,vz,wx,wy,wz;
	for(i=0;i<p;i++) b[i]=true;
	for(i=3;i<p;i+=3) {
		for(j=0;j<3;j++) {
			ux=pts[0]-pts[i];
			uy=pts[1]-pts[i+1];
			uz=pts[2]-pts[i+2];
			if (b[i+j]) {
				b[i+j]=false;
				k=ed[i+j];
				l=vor_up(rl[i+j]);
				vx=pts[k]-pts[0];
				vy=pts[k+1]-pts[1];
				vz=pts[k+2]-pts[2];
				b[k+l]=false;
				m=ed[k+l];
				while(m!=i) {
					n=vor_up(rl[k+l]);
					wx=pts[m]-pts[0];
					wy=pts[m+1]-pts[1];
					wz=pts[m+2]-pts[2];
					vol+=ux*vy*wz+uy*vz*wx+uz*vx*wy-uz*vy*wx-uy*vx*wz-ux*vz*wy;
					b[m+n]=false;
					k=m;l=n;vx=wx;vy=wy;vz=wz;
					m=ed[k+l];
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
	for(i=0;i<p;i+=3) {
		ux=x+0.5*pts[i];uy=y+0.5*pts[i+1];uz=z+0.5*pts[i+2];
		of << "sphere{<" << ux << "," << uy << "," << uz << ">,r}" << endl;
		for(j=0;j<3;j++) {
			k=ed[i+j];
			if (k<i) of << "cylinder{<" << ux << "," << uy << "," << uz << ">,<" << x+0.5*pts[k] << "," << y+0.5*pts[k+1] << "," << z+0.5*pts[k+2] << ">,r}" << endl;
		}
	}
};

// Outputs the edges of the Voronoi cell (in gnuplot format) to an open file
// stream, displacing the cell by an amount (x,y,z)
inline void voronoicell::dumpgnuplot(ofstream &of,double x,double y,double z) {
	int i,j,k;double ux,uy,uz;
	for(i=0;i<p;i+=3) {
		ux=x+0.5*pts[i];uy=y+0.5*pts[i+1];uz=z+0.5*pts[i+2];
		for(j=0;j<3;j++) {
			k=ed[i+j];
			if (ed[i+j]<i) of << ux << " " << uy << " " << uz << endl << x+0.5*pts[k] << " " << y+0.5*pts[k+1] << " " << z+0.5*pts[k+2] << endl << endl << endl;
		}
	}
};
