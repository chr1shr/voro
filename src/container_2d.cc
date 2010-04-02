

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (xa,xb) the minimum and maximum x coordinates.
 * \param[in] (ya,yb) the minimum and maximum y coordinates.
 * \param[in] (xn,yn) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] memi the initial memory allocation for each block. */
container_2d::container_base(fpoint xa,fpoint xb,fpoint ya,
		fpoint yb,int xn,int yn,int memi)
	: ax(xa),bx(xb),ay(ya),by(yb),
	xsp(xn/(xb-xa)),ysp(yn/(yb-ya)),nx(xn),ny(yn),
	nxy(xn*yn),co(new int[nxyz]),mem(new int[nxyz]),
	id(new int*[nxyz]),p(new fpoint*[nxyz]) {
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
template<class r_option>
void container_base<r_option>::add_particle_memory(int i) {
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

/** Clears a container of particles. */
void container_2d::clear() {
	for(int ij=0;ij<nxy;ij++) co[ij]=0;
}


/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in gnuplot format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box. */
template<class r_option>
void container_base<r_option>::draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax) {
	fpoint x,y,px,py;
	voropp_loop l1(this);
	int q,s;
	voronoicell_2d c;  //BIG FIXME
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,px,py);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][2*q]+px;y=p[s][2*q+1]+py;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax) {
				if(compute_cell_sphere(c,l1.ip,l1.jp,s,q,x,y)) c.draw_gnuplot(os,x,y);
			}
		}
	} while((s=l1.inc(px,py))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_gnuplot() that computes the Voronoi
 * cells for the entire simulation region and saves the output in gnuplot
 * format.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_base<r_option>::draw_cells_gnuplot(const char *filename) {
	draw_cells_gnuplot(filename,ax,bx,ay,by);
}

inline bool container_2d::initialize_voronoicell(voronoicell_2d &c,fpoint x,fpoint y) {
	fpoint x1,x2,y1,y2;
	if(xperiodic) x1=-(x2=0.5*(bx-ax));else {x1=ax-x;x2=bx-x;}
	if(yperiodic) y1=-(y2=0.5*(by-ay));else {y1=ay-y;y2=by-y;}
	c.init(x1,x2,y1,y2);
	return true;
}

bool container_base<r_option>::compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int ij,int s,fpoint x,fpoint y) {

	// This length scale determines how large the spherical shells should
	// be, and it should be set to approximately the particle diameter
	const fpoint length_scale=1;
	fpoint x1,y1,qx,qy,lr=0,lrs=0,ur,urs,rs;
	int q,t;
	voropp_loop_2d l(this);
	if(!initialize_voronoicell(c,x,y)) return false;

	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing.
	radius.init(ij,s);
	while(radius.cutoff(lrs)<c.max_radius_squared()) {
		ur=lr+0.5*length_scale;urs=ur*ur;
		t=l.init(x,y,ur,qx,qy);
		do {
			for(q=0;q<co[t];q++) {
				x1=p[t][sz*q]+qx-x;y1=p[t][sz*q+1]+qy-y;
				rs=x1*x1+y1*y1+z1*z1;
				if(lrs-tolerance<rs&&rs<urs&&(q!=s||ij!=t)) {
					if(!c.nplane(x1,y1,rs)) return false;
				}
			}
		} while((t=l.inc(qx,qy))!=-1);
		lr=ur;lrs=urs;
	}
	return true;
}




