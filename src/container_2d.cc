/** \file container_2d.cc
 * \brief Function implementations for the container_2d class. */

#include "container_2d.hh"

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (xa,xb) the minimum and maximum x coordinates.
 * \param[in] (ya,yb) the minimum and maximum y coordinates.
 * \param[in] (xn,yn) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] (xper,yper) flags setting whether the container is periodic
 *                        in each coordinate direction.
 * \param[in] memi the initial memory allocation for each block. */
container_2d::container_2d(fpoint xa,fpoint xb,fpoint ya,
		fpoint yb,int xn,int yn,bool xper,bool yper,int memi)
	: ax(xa),bx(xb),ay(ya),by(yb),xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),
	nx(xn),ny(yn),nxy(xn*yn),xperiodic(xper),yperiodic(yper),
	co(new int[nxy]),mem(new int[nxy]),id(new int*[nxy]),p(new fpoint*[nxy]) {
	int l;
	for(l=0;l<nxy;l++) co[l]=0;
	for(l=0;l<nxy;l++) mem[l]=memi;
	for(l=0;l<nxy;l++) id[l]=new int[memi];
	for(l=0;l<nxy;l++) p[l]=new fpoint[2*memi];
}

/** The container destructor frees the dynamically allocated memory. */
container_2d::~container_2d() {
	int l;
	for(l=0;l<nxy;l++) delete [] p[l];
	for(l=0;l<nxy;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] mem;
	delete [] co;
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] os an output stream to write to. */
void container_2d::draw_particles(ostream &os) {
	int c,l,i;
	for(l=0;l<nxy;l++) for(c=0;c<co[l];c++) {
		os << id[l][c];
		for(i=2*c;i<2*(c+1);i++) os << " " << p[l][i];
		os << "\n";
	}
}

/** An overloaded version of the draw_particles() routine, that just prints
 * to standard output. */
void container_2d::draw_particles() {
	draw_particles(cout);
}

/** An overloaded version of the draw_particles() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
void container_2d::draw_particles(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles(os);
	os.close();
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle. */
void container_2d::put(int n,fpoint x,fpoint y) {
	if(x>ax&&y>ay) {
		int i,j;
		i=int((x-ax)*xsp);j=int((y-ay)*ysp);
		if(i<nx&&j<ny) {
			i+=nx*j;
			if(co[i]==mem[i]) add_particle_memory(i);
			p[i][2*co[i]]=x;p[i][2*co[i]+1]=y;
			id[i][co[i]++]=n;
		}
	}
}


/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_2d::add_particle_memory(int i) {
	int *idp;fpoint *pp;
	int l,nmem(2*mem[i]);
#if VOROPP_VERBOSE >=3
	cerr << "Particle memory in region " << i << " scaled up to " << nmem << endl;
#endif
	if(nmem>max_particle_memory)
		voropp_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new fpoint[2*nmem];
	for(l=0;l<2*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}

/** An overloaded version of the import routine, that reads the standard input.
 */
inline void container_2d::import() {
	import(cin);
}

/** An overloaded version of the import routine, that reads in particles from
 * a particular file.
 * \param[in] filename the name of the file to read from. */
inline void container_2d::import(const char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	if(is.fail()) voropp_fatal_error("Unable to open file for import",VOROPP_FILE_ERROR);
	import(is);
	is.close();
}

/** Imports a list of particles from an input stream.
 * \param[in] is an input stream to read from. */
inline void container_2d::import(istream &is) {
	int n;fpoint x,y;
	is >> n >> x >> y;
	while(!is.eof()) {
		put(n,x,y);
		is >> n >> x >> y;
	}
}

/** Clears a container of particles. */
void container_2d::clear() {
	for(int ij=0;ij<nxy;ij++) co[ij]=0;
}


/** Computes the Voronoi cells for all particles within a rectangular box, and
 * saves the output in gnuplot format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box. */
void container_2d::draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax) {
	fpoint x,y,px,py;
	voropp_loop_2d l1(this);
	int q,s;
	voronoicell_2d c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,px,py);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][2*q]+px;y=p[s][2*q+1]+py;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax) {
				if(compute_cell_sphere(c,l1.ip,l1.jp,s,q,x,y)) {
					c.draw_gnuplot(os,x,y);os << endl;
				}
			}
		}
	} while((s=l1.inc(px,py))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_gnuplot() that computes the Voronoi
 * cells for the entire simulation region and saves the output in gnuplot
 * format.
 * \param[in] filename the name of the file to write to. */
void container_2d::draw_cells_gnuplot(const char *filename) {
	draw_cells_gnuplot(filename,ax,bx,ay,by);
}

/** Computes the Voronoi cells for all particles within a rectangular box, and
 * saves the output in POV-Ray format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box. */
void container_2d::draw_cells_pov(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax) {
	fpoint x,y,px,py;
	voropp_loop_2d l1(this);
	int q,s;
	voronoicell_2d c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,px,py);
	do {
		for(q=0;q<co[s];q++) {
			os << "// cell " << id[s][q] << "\n";
			x=p[s][2*q]+px;y=p[s][2*q+1]+py;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax) {
				if(compute_cell_sphere(c,l1.ip,l1.jp,s,q,x,y)) {
					c.draw_pov(os,x,y);os << endl;
				}
			}
		}
	} while((s=l1.inc(px,py))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_pov() that computes the Voronoi
 * cells for the entire simulation region and saves the output in POV-Ray
 * format.
 * \param[in] filename the name of the file to write to. */
void container_2d::draw_cells_pov(const char *filename) {
	draw_cells_pov(filename,ax,bx,ay,by);
}

/** Dumps all the particle positions in POV-Ray format.
 * \param[in] os an output stream to write to. */
void container_2d::draw_particles_pov(ostream &os) {
	int l,c;
	for(l=0;l<nxy;l++) for(c=0;c<co[l];c++) {
		os << "// id " << id[l][c] << "\n";
		os << "sphere{<" << p[l][2*c] << "," << p[l][2*c+1] << ",0>,s}\n";
	}
}

/** An overloaded version of the draw_particles_pov() routine, that just prints
 * to standard output. */
void container_2d::draw_particles_pov() {
	draw_particles_pov(cout);
}

/** An overloaded version of the draw_particles_pov() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
void container_2d::draw_particles_pov(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles_pov(os);
	os.close();
}

/** Computes the Voronoi cells for all particles in the container, and for each
 * cell, outputs a line containing custom information about the cell structure.
 * The output format is specified using an input string with control sequences
 * similar to the standard C printf() routine.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] os an open output stream to write to. */
void container_2d::print_all_custom(const char *format,ostream &os) {
	int i,j,ij=0,fp,q;
	fpoint x,y;
	voronoicell_2d c;
	for(j=0;j<ny;j++) for(i=0;i<nx;i++,ij++) for(q=0;q<co[ij];q++) {
		x=p[ij][2*q];y=p[ij][2*q+1];
		if(!compute_cell_sphere(c,i,j,ij,q,x,y)) continue;
		fp=0;
		while(format[fp]!=0) {
			if(format[fp]=='%') {
				fp++;
				switch(format[fp]) {

					// Particle-related output
					case 'i': os << id[ij][q];break;
					case 'x': os << x;break;
					case 'y': os << y;break;
					case 'q': os << x << " " << y;break;

					// Vertex-related output
					case 'w': os << c.p;break;
					case 'm': os << 0.25*c.max_radius_squared();break;

					// Edge-related output
					case 'p': os << c.perimeter();break;

					// Area-related output
					case 'v': os << c.area();break;
					case 'c': {
							  fpoint cx,cy;
							  c.centroid(cx,cy);
							  os << cx << " " << cy;
						  } break;
					case 'C': {
							  fpoint cx,cy;
							  c.centroid(cx,cy);
							  os << x+cx << " " << y+cy;
						  } break;

					// End-of-string reached
					case 0: fp--;break;

					// The percent sign is not part of a
					// control sequence
					default: os << '%' << format[fp];
				}
			} else os << format[fp];
			fp++;
		}
		os << "\n";
	}
}

/** An overloaded version of the print_all_custom routine, that prints to the
 * standard output.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics. */
inline void container_2d::print_all_custom(const char *format) {
	print_all_custom(format,cout);
}

/** An overloaded version of the print_all_custom routine that writes directly to a
 * file.
 * \param[in] filename The name of the file to write to.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics. */
inline void container_2d::print_all_custom(const char *format,const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_custom(format,os);
	os.close();
}

/** Initializes a voronoicell_2d class to fill the entire container.
 * \param[in] c a reference to a voronoicell_2d class.
 * \param[in] (x,y) the position of the particle that . */
inline bool container_2d::initialize_voronoicell(voronoicell_2d &c,fpoint x,fpoint y) {
	fpoint x1,x2,y1,y2;
	if(xperiodic) x1=-(x2=0.5*(bx-ax));else {x1=ax-x;x2=bx-x;}
	if(yperiodic) y1=-(y2=0.5*(by-ay));else {y1=ay-y;y2=by-y;}
	c.init(x1,x2,y1,y2);
	return true;
}

/** An overloaded version of the compute_cell_sphere routine, that sets up the x
 * and y variables.
 \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j) the coordinates of the block that the test particle is
 *                  in.
 * \param[in] ij the index of the block that the test particle is in, set to
 *               i+nx*j.
 * \param[in] s the index of the particle within the test block.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
inline bool container_2d::compute_cell_sphere(voronoicell_2d &c,int i,int j,int ij,int s) {
	fpoint x=p[s][2*ij],y=p[s][2*ij+1];
	return compute_cell_sphere(c,i,j,ij,s,x,y);
}

/** This routine computes the Voronoi cell for a give particle, by successively
 * testing over particles within larger and larger concentric circles. This
 * routine is simple and fast, although it may not work well for anisotropic
 * arrangements of particles.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j) the coordinates of the block that the test particle is
 *                  in.
 * \param[in] ij the index of the block that the test particle is in, set to
 *               i+nx*j.
 * \param[in] s the index of the particle within the test block.
 * \param[in] (x,y) the coordinates of the particle.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
bool container_2d::compute_cell_sphere(voronoicell_2d &c,int i,int j,int ij,int s,fpoint x,fpoint y) {

	// This length scale determines how large the spherical shells should
	// be, and it should be set to approximately the particle diameter
	const fpoint length_scale=0.5*sqrt((bx-ax)*(by-ay)/nx/ny);
	
	fpoint x1,y1,qx,qy,lr=0,lrs=0,ur,urs,rs;
	int q,t;
	voropp_loop_2d l(this);

	if(!initialize_voronoicell(c,x,y)) return false;

	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing.
	while(lrs<c.max_radius_squared()) {
		ur=lr+0.5*length_scale;urs=ur*ur;
		t=l.init(x,y,ur,qx,qy);
		do {
			for(q=0;q<co[t];q++) {
				x1=p[t][2*q]+qx-x;y1=p[t][2*q+1]+qy-y;
				rs=x1*x1+y1*y1;
				if(lrs-tolerance<rs&&rs<urs&&(q!=s||ij!=t)) {
					if(!c.plane(x1,y1,rs)) return false;
				}
			}
		} while((t=l.inc(qx,qy))!=-1);
		lr=ur;lrs=urs;
	}
	return true;
}


/** Creates a voropp_loop_2d object, by setting the necessary constants about the
 * container geometry from a pointer to the current container class.
 * \param[in] q a pointer to the current container class. */
voropp_loop_2d::voropp_loop_2d(container_2d *q) : sx(q->bx-q->ax), sy(q->by-q->ay),
	xsp(q->xsp),ysp(q->ysp),
	ax(q->ax),ay(q->ay),nx(q->nx),ny(q->ny),nxy(q->nxy),
	xperiodic(q->xperiodic),yperiodic(q->yperiodic) {}

/** Initializes a voropp_loop_2d object, by finding all blocks which are within a
 * given sphere. It calculates the index of the first block that needs to be
 * tested and sets the periodic displacement vector accordingly.
 * \param[in] (vx,vy) the position vector of the center of the sphere.
 * \param[in] r the radius of the sphere.
 * \param[out] (px,py) the periodic displacement vector for the first block to
 *                     be tested.
 * \return The index of the first block to be tested. */
inline int voropp_loop_2d::init(fpoint vx,fpoint vy,fpoint r,fpoint &px,fpoint &py) {
	ai=step_int((vx-ax-r)*xsp);
	bi=step_int((vx-ax+r)*xsp);
	if(!xperiodic) {
		if(ai<0) {ai=0;if(bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if(ai>=nx) ai=nx-1;}
	}
	aj=step_int((vy-ay-r)*ysp);
	bj=step_int((vy-ay+r)*ysp);
	if(!yperiodic) {
		if(aj<0) {aj=0;if(bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if(aj>=ny) aj=ny-1;}
	}
	i=ai;j=aj;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	inc1=aip-step_mod(bi,nx)+nx;
	s=aip+nx*ajp;
	return s;
}

/** Initializes a voropp_loop_2d object, by finding all blocks which overlap a given
 * rectangular box. It calculates the index of the first block that needs to be
 * tested and sets the periodic displacement vector (px,py,pz) accordingly.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[out] (px,py) the periodic displacement vector for the first block
 *                     to be tested.
 * \return The index of the first block to be tested. */
inline int voropp_loop_2d::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint &px,fpoint &py) {
	ai=step_int((xmin-ax)*xsp);
	bi=step_int((xmax-ax)*xsp);
	if(!xperiodic) {
		if(ai<0) {ai=0;if(bi<0) bi=0;}
		if(bi>=nx) {bi=nx-1;if(ai>=nx) ai=nx-1;}
	}
	aj=step_int((ymin-ay)*ysp);
	bj=step_int((ymax-ay)*ysp);
	if(!yperiodic) {
		if(aj<0) {aj=0;if(bj<0) bj=0;}
		if(bj>=ny) {bj=ny-1;if(aj>=ny) aj=ny-1;}
	}
	i=ai;j=aj;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	ajp=jp=step_mod(j,ny);apy=py=step_div(j,ny)*sy;
	inc1=aip-step_mod(bi,nx)+nx;
	s=aip+nx*ajp;
	return s;
}

/** Returns the next block to be tested in a loop, and updates the periodicity
 * vector if necessary.
 * \param[in,out] (px,py) the current block on entering the function, which is
 *                        updated to the next block on exiting the function. */
inline int voropp_loop_2d::inc(fpoint &px,fpoint &py) {
	if(i<bi) {
		i++;
		if(ip<nx-1) {ip++;s++;} else {ip=0;s+=1-nx;px+=sx;}
		return s;
	} else if(j<bj) {
		i=ai;ip=aip;px=apx;j++;
		if(jp<ny-1) {jp++;s+=inc1;} else {jp=0;s+=inc1-nxy;py+=sy;}
		return s;
	} else return -1;
}

/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1). */
inline int voropp_loop_2d::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
}

/** Custom modulo function, that gives consistent stepping for negative
 * numbers. */
inline int voropp_loop_2d::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

/** Custom integer division function, that gives consistent stepping for
 * negative numbers. */
inline int voropp_loop_2d::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}
