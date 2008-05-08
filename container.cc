// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

/** Container constructor. The first six arguments set the corners of the box to
 * be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
 * grid of blocks, set by the following three arguments. The next three
 * arguments are booleans, which set the periodicity in each direction. The
 * final argument sets the amount of memory allocated to each block. */
container::container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi,int isz)
	: sz(isz), ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
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
	p=new fpoint*[nxyz];
	for(l=0;l<nxyz;l++) p[l]=new fpoint[sz*memi];
};

/** Container constructor. The first six arguments set the corners of the box to
 * be (xa,ya,za) and (xb,yb,zb). The box is then divided into an nx by ny by nz
 * grid of blocks, set by the following three arguments. The next three
 * arguments are booleans, which set the periodicity in each direction. The
 * final argument sets the amount of memory allocated to each block. */
container::container(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi)
	: sz(3), ax(xa),bx(xb),ay(ya),by(yb),az(za),bz(zb),
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
	p=new fpoint*[nxyz];
	for(l=0;l<nxyz;l++) p[l]=new fpoint[sz*memi];
};

/** Container destructor - free memory. */
container::~container() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] mem;
	delete [] co;	
};

/** Dumps all the particle positions and identifies to a file. */
void container::dump(char *filename) {
	int c,l,i;
	ofstream file;
	file.open(filename,ofstream::out|ofstream::trunc);
	for(l=0;l<nxyz;l++) {
		for (c=0;c<co[l];c++)
			file << id[l][c];
			for(i=sz*c;i<sz*(c+1);i++) file << " " << p[l][i];
			file << endl;
	}
	file.close();
};

/** Put a particle into the correct region of the container. */
void container::put(int n,fpoint x,fpoint y,fpoint z) {
	if(x>ax&&y>ay&&z>az) {
		int i,j,k;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);k=int((z-az)*zsp);
		if(i<nx&&j<ny&&k<nz) {
			i+=nx*j+nxy*k;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][3*co[i]]=x;p[i][3*co[i]+1]=y;p[i][3*co[i]+2]=z;
			id[i][co[i]++]=n;
		}
	}
};

/** Increase memory for a particular region. */
void container::add_particle_memory(int i) {
	int *idp;fpoint *pp;
	int l,nmem=2*mem[i];
	if (nmem>maxparticlemem) throw fatal_error("Absolute maximum memory allocation exceeded");
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new fpoint[sz*nmem];
	for(l=0;l<sz*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
};

/** Import a list of particles from standard input. */
void container::import(istream &is) {
	int n;fpoint x,y,z;
	is >> n >> x >> y >> z;
	while(!is.eof()) {
		put(n,x,y,z);
		is >> n >> x >> y >> z;
	}
};

/** An overloaded version of the import routine, that reads the standard input.
 */
inline void container::import() {
	import(cin);
};

/** An overloaded version of the import routine, that reads in particles from
 * <filename>. */
inline void container::import(char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	import(is);
	is.close();
};

/** Outputs the number of particles within each region. */
void container::region_count() {
	int i,j,k,ijk=0;
	for(k=0;k<nz;k++) {
		for(j=0;j<ny;j++) {
			for(i=0;i<nx;i++) cout << "Region (" << i << "," << j << "," << k << "): " << co[ijk++] << " particles" << endl;
		}
	}
};

/** Clears a container of particles. */
void container::clear() {
	for(int ijk=0;ijk<nxyz;ijk++) co[ijk]=0;
	poly_clear_radius();
};

/** Computes the Voronoi cells for all particles within a box with corners
 * (xmin,ymin,zmin) and (xmax,ymax,zmax), and saves the output in a format
 * that can be read by gnuplot. */
void container::draw_gnuplot(char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	facets_loop l1(this);
	int i,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(i=0;i<co[s];i++) {
			x=p[s][sz*i]+px;y=p[s][sz*i+1]+py;z=p[s][sz*i+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				compute_cell(c,s,i,x,y,z);
				c.dump_gnuplot(os,x,y,z);
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	os.close();
};

/** If only a filename is supplied to draw_gnuplot(), then assume that we are
 * calculating the entire simulation region. */
void container::draw_gnuplot(char *filename) {
	draw_gnuplot(filename,ax,bx,ay,by,az,bz);
};

/** Computes the Voronoi cells for all particles within a box with corners
 * (xmin,ymin,zmin) and (xmax,ymax,zmax), and saves the output in a format
 * that can be read by gnuplot.*/
void container::draw_pov(char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px,py,pz;
	facets_loop l1(this);
	int i,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	os << "#declare voronoi=union{\n";
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz);
	do {
		for(i=0;i<co[s];i++) {
			x=p[s][sz*i]+px;y=p[s][sz*i+1]+py;z=p[s][sz*i+2]+pz;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				compute_cell(c,s,i,x,y,z);
				c.dump_pov(os,x,y,z);break;
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
	os << "}\n";
	os.close();
};

/** If only a filename is supplied to draw_pov(), then assume that we are
 * calculating the entire simulation region.*/
void container::draw_pov(char *filename) {
	draw_pov(filename,ax,bx,ay,by,az,bz);
};


/** Computes the Voronoi volumes for all the particles, and stores the
 * results according to the particle label in the fpoint array bb.*/
void container::store_cell_volumes(fpoint *bb) {
	voronoicell c;
	facets_loop l(this);
	int i,s;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			compute_cell(c,s,i);
			bb[id[s][i]]=c.volume();
		}
	}
};

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output. */
inline void container::print_all(ostream &os,voronoicell &c) {
	fpoint x,y,z;
	facets_loop l(this);
	int i,s;
	for(s=0;s<nxyz;s++) {
		for(i=0;i<co[s];i++) {
			x=p[s][sz*i];y=p[s][sz*i+1];z=p[s][sz*i+2];
			compute_cell(c,s,i,x,y,z);
			os << id[s][i] << " " << x << " " << y << " " << z;
			if (sz==4) cout << " " << p[s][4*i+3];
			os << " " << c.volume();
			c.neighbors(os);			
			os << endl;
		}
	}
};

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output. */
void container::print_all(ostream &os) {
	voronoicell c;
	print_all(os,c);
};

/** An overloaded version of print_all(), which just prints to standard output. */
void container::print_all() {
	voronoicell c;
	print_all(cout);
};

/** An overloaded version of print_all(), which outputs the result to <filename>. */
inline void container::print_all(char* filename) {
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all(os,c);
	os.close();
};

/** Prints a list of all particle labels, positions, Voronoi volumes, and a list
 * of neighboring particles to an output stream.
 * \param[in] os The output stream to print to.*/
void container::print_all_neighbor(ostream &os) {
	voronoicell_neighbor c;
	print_all(os,c);
};

/** An overloaded version of print_all_neighbor(), which just prints to standard output. */
void container::print_all_neighbor() {
	voronoicell_neighbor c;
	print_all(cout,c);
};

/** An overloaded version of print_all_neighbor(), which outputs the result to <filename>. */
inline void container::print_all_neighbor(char* filename) {
	voronoicell_neighbor c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all(os,c);
	os.close();
};

/** Initialize the Voronoi cell to be the entire container. For non-periodic
 * coordinates, this is set by the position of the walls. For periodic
 * coordinates, the space is equally divided in either direction from the
 * particle's initial position. That makes sense since those boundaries would
 * be made by the neighboring periodic images of this particle. */
inline void container::initialize_voronoicell(voronoicell &c,fpoint x,fpoint y,fpoint z) {
	float x1,x2,y1,y2,z1,z2;
	if (xperiodic) x1=-(x2=0.5*(bx-ax));else {x1=ax-x;x2=bx-x;}
	if (yperiodic) y1=-(y2=0.5*(by-ay));else {y1=ay-y;y2=by-y;}
	if (zperiodic) z1=-(z2=0.5*(bz-az));else {z1=az-z;z2=bz-z;}
	c.init(x1,x2,y1,y2,z1,z2);
};

/** Computes a single Voronoi cell in the container. This routine can be run by
 * the user, and it is also called multiple times by the functions vprintall,
 * store_cell_volumes() and draw(). */
inline void container::compute_cell(voronoicell &c,int s,int i,fpoint x,fpoint y,fpoint z) {
	fpoint x1,y1,z1,qx,qy,qz,lr=0,lrs=0,ur,urs,rs;
	int j,t;
	facets_loop l(this);
	initialize_voronoicell(c,x,y,z);

	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing. TODO: this can sometimes be inefficient.
	// For example, sometimes particles at the top of granular packings can
	// extend upwards by a long way, and the shells grow very big. It would
	// be better to use a box-by-box approach, but that's not
	// straightforward.
	while(lrs<c.maxradsq()) {
		ur=lr+0.5;urs=ur*ur;
		t=l.init(x,y,z,ur,qx,qy,qz);
		do {
			for(j=0;j<co[t];j++) {
				x1=p[t][sz*j]+qx-x;y1=p[t][sz*j+1]+qy-y;z1=p[t][sz*j+2]+qz-z;
				rs=x1*x1+y1*y1+z1*z1;
				if (lrs-tolerance<rs&&rs<urs&&(j!=i||s!=t))
					c.nplane(x1,y1,z1,rs,id[t][j]);
			}
		} while ((t=l.inc(qx,qy,qz))!=-1);
		lr=ur;lrs=urs;
	}
};

/** A overloaded version of compute_cell, that sets up the x, y, and z variables. */
inline void container::compute_cell(voronoicell &c,int s,int i) {
	double x=p[s][sz*i],y=p[s][sz*i+1],z=p[s][sz*i+2];
	compute_cell(c,s,i,x,y,z);
}

/** Creates a facets_loop object, by pulling the necesssary constants about the container
 * geometry from a pointer to the current container class. */
facets_loop::facets_loop(container *q) : sx(q->bx-q->ax), sy(q->by-q->ay), sz(q->bz-q->az),
	xsp(q->xsp),ysp(q->ysp),zsp(q->zsp),
	ax(q->ax),ay(q->ay),az(q->az),
	nx(q->nx),ny(q->ny),nz(q->nz),nxy(q->nxy),nxyz(q->nxyz),
	xperiodic(q->xperiodic),yperiodic(q->yperiodic),zperiodic(q->zperiodic) {};

/** Initializes a facets_loop object, by finding all blocks which are within a distance
 * r of the vector (vx,vy,vz). It returns the first block which is to be
 * tested, and sets the periodic displacement vector (px,py,pz) accordingly. */
inline int facets_loop::init(fpoint vx,fpoint vy,fpoint vz,fpoint r,fpoint &px,fpoint &py,fpoint &pz) {
	ai=step_int((vx-ax-r)*xsp);
	bi=step_int((vx-ax+r)*xsp);
	if (!xperiodic) {
		if(ai<0) {ai=0;if (bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if (ai>=nx) ai=nx-1;} 
	}
	aj=step_int((vy-ay-r)*ysp);
	bj=step_int((vy-ay+r)*ysp);
	if (!yperiodic) {
		if(aj<0) {aj=0;if (bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if (aj>=ny) aj=ny-1;} 
	}
	ak=step_int((vz-az-r)*zsp);
	bk=step_int((vz-az+r)*zsp);
	if (!zperiodic) {
		if(ak<0) {ak=0;if (bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if (ak>=nz) ak=nz-1;} 
	}
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	akp=kp=step_mod(k,nz);apz=pz=step_div(k,nz)*sz;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(ny+ajp-step_mod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
};

/** Initializes a facets_loop object, by finding all blocks which overlap the box with
 * corners (xmin,ymin,zmin) and (xmax,ymax,zmax). It returns the first block
 * which is to be tested, and sets the periodic displacement vector (px,py,pz)
 * accordingly. */
inline int facets_loop::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax,fpoint &px,fpoint &py,fpoint &pz) {
	ai=step_int((xmin-ax)*xsp);
	bi=step_int((xmax-ax)*xsp);
	if (!xperiodic) {
		if(ai<0) {ai=0;if (bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if (ai>=nx) ai=nx-1;} 
	}
	aj=step_int((ymin-ay)*ysp);
	bj=step_int((ymax-ay)*ysp);
	if (!yperiodic) {
		if(aj<0) {aj=0;if (bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if (aj>=ny) aj=ny-1;} 
	}
	ak=step_int((zmin-az)*zsp);
	bk=step_int((zmax-az)*zsp);
	if (!zperiodic) {
		if(ak<0) {ak=0;if (bk<0) bk=0;}
		if(bk>=nz) {bk=nz-1;if (ak>=nz) ak=nz-1;} 
	}
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	akp=kp=step_mod(k,nz);apz=pz=step_div(k,nz)*sz;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(ny+ajp-step_mod(bj,ny))+inc1;
	inc1+=nx;
	s=aip+nx*(ajp+ny*akp);
	return s;
};

/** Returns the next block to be tested in a loop, and updates the periodicity
 * vector if necessary. */
inline int facets_loop::inc(fpoint &px,fpoint &py,fpoint &pz) {
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

/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).*/
inline int facets_loop::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
};

/** Custom mod function, that gives consistent stepping for negative numbers. */
inline int facets_loop::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
};

/** Custom div function, that gives consistent stepping for negative numbers. */
inline int facets_loop::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
};
