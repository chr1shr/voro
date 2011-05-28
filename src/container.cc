// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file container.cc
 * \brief Function implementations for the container and related classes. */

#include "container.hh"

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
 *                                                container is periodic in each
 *                                                coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_base::container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
		int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem,int ps_)
	: voropp_base(nx_,ny_,nz_,(bx_-ax_)/nx_,(by_-ay_)/ny_,(bz_-az_)/nz_),
	ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
	xperiodic(xperiodic_), yperiodic(yperiodic_), zperiodic(zperiodic_),
	id(new int*[nxyz]), p(new double*[nxyz]), co(new int[nxyz]), mem(new int[nxyz]), ps(ps_) {
	int l;
	for(l=0;l<nxyz;l++) co[l]=0;
	for(l=0;l<nxyz;l++) mem[l]=init_mem;
	for(l=0;l<nxyz;l++) id[l]=new int[init_mem];
	for(l=0;l<nxyz;l++) p[l]=new double[ps*init_mem];
}

/** The container destructor frees the dynamically allocated memory. */
container_base::~container_base() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] mem;
	delete [] co;
}

/** The class constructor sets up the geometry of container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
 *                                                container is periodic in each
 *                                                coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container::container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
	int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem)
	: container_base(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,xperiodic_,yperiodic_,zperiodic_,init_mem,3),
	vc(*this,xperiodic_?2*nx_+1:nx_,yperiodic_?2*ny_+1:ny_,zperiodic_?2*nz_+1:nz_) {}

/** The class constructor sets up the geometry of container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
 *                                                container is periodic in each
 *                                                coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_poly::container_poly(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
	int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem)
	: container_base(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,xperiodic_,yperiodic_,zperiodic_,init_mem,4),
	max_radius(0), vc(*this,xperiodic_?2*nx_+1:nx_,yperiodic_?2*ny_+1:ny_,zperiodic_?2*nz_+1:nz_) {}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container::put(int n,double x,double y,double z) {
	int ijk;
	if(put_locate_block(ijk,x,y,z)) {
		id[ijk][co[ijk]]=n;
		double *pp(p[ijk]+3*co[ijk]++);
		*(pp++)=x;*(pp++)=y;*pp=z;
	}
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly::put(int n,double x,double y,double z,double r) {
	int ijk;
	if(put_locate_block(ijk,x,y,z)) {
		id[ijk][co[ijk]]=n;
		double *pp(p[ijk]+4*co[ijk]++);
		*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
		if(max_radius<r) max_radius=r;
	}
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container::put(voropp_order &vo,int n,double x,double y,double z) {
	int ijk;
	if(put_locate_block(ijk,x,y,z)) {
		id[ijk][co[ijk]]=n;
		vo.add(ijk,co[ijk]);
		double *pp(p[ijk]+3*co[ijk]++);
		*(pp++)=x;*(pp++)=y;*pp=z;
	}
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly::put(voropp_order &vo,int n,double x,double y,double z,double r) {
	int ijk;
	if(put_locate_block(ijk,x,y,z)) {
		id[ijk][co[ijk]]=n;
		vo.add(ijk,co[ijk]);
		double *pp(p[ijk]+4*co[ijk]++);
		*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
		if(max_radius<r) max_radius=r;
	}
}

/** This routine takes a particle position vector, tries to remap it into the
 * primary domain. If successful, it computes the region into which it can be
 * stored and checks that there is enough memory within this region to store
 * it.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base::put_locate_block(int &ijk,double &x,double &y,double &z) {
	if(put_remap(ijk,x,y,z)) {
		if(co[ijk]==mem[ijk]) add_particle_memory(ijk);
		return true;
	}
#if VOROPP_REPORT_OUT_OF_BOUNDS
	fprintf(stderr,"Out of bounds: [%d] (x,y,z)=(%g,%g,%g)\n",n,x,y,z);
#endif
	return false;
}

/** Takes a particle position vector and computes the region index into which
 * it should be stored. If the container is periodic, then the routine also
 * maps the particle position to ensure it is in the primary domain. If the
 * container is not periodic, the routine bails out. 
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base::put_remap(int &ijk,double &x,double &y,double &z) {
	int l;
	
	ijk=step_int((x-ax)*xsp);
	if(xperiodic) {l=step_mod(ijk,nx);x+=boxx*(l-ijk);ijk=l;}
	else if(ijk<0||ijk>=nx) return false;

	int j(step_int((y-ay)*ysp));
	if(yperiodic) {l=step_mod(j,ny);y+=boxy*(l-j);j=l;}
	else if(j<0||j>=ny) return false;

	int k(step_int((z-az)*zsp));
	if(xperiodic) {l=step_mod(k,nz);z+=boxz*(l-k);k=l;}
	else if(k<0||k>=nz) return false;

	ijk+=nx*j+nxy*k;
	return true;
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_base::add_particle_memory(int i) {
	int *idp,l,nmem(mem[i]<<1);
	double *pp;
#if VOROPP_VERBOSE >=3
	fprintf(stderr,"Particle memory in region %d scaled up to %d\n",i,nmem);
#endif
	if(nmem>max_particle_memory)
		voropp_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new double[ps*nmem];
	for(l=0;l<ps*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}

/** Import a list of particles from an open file stream into the container.
 * Entries of four numbers (Particle ID, x position, y position, z position)
 * are searched for. If the file cannot be successfully read, then the routine
 * causes a fatal error.
 * \param[fp] fp the file handle to read from. */
void container::import(FILE *fp) {
	int i,j;
	double x,y,z;
	while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(i,x,y,z);
	if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position) are searched for. If the file cannot be
 * successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container::import(voropp_order &vo,FILE *fp) {
	int i,j;
	double x,y,z;
	while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(vo,i,x,y,z);
	if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream into the container.
 * Entries of five numbers (Particle ID, x position, y position, z position,
 * radius) are searched for. If the file cannot be successfully read, then the
 * routine causes a fatal error.
 * \param[fp] fp the file handle to read from. */
void container_poly::import(FILE *fp) {
	int i,j;
	double x,y,z,r;
	while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(i,x,y,z,r);
	if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position, radius) are searched for. If the file
 * cannot be successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_poly::import(voropp_order &vo,FILE *fp) {
	int i,j;
	double x,y,z,r;
	while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(vo,i,x,y,z,r);
	if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Outputs the a list of all the container regions along with the number of
 * particles stored within each. */
void container_base::region_count() {
	int i,j,k,*cop(co);
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
		printf("Region (%d,%d,%d): %d particles\n",i,j,k,*(cop++));
}

/** Clears a container of particles. */
void container::clear() {
	for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
}

/** Clears a container of particles, also clearing resetting the maximum radius
 * to zero. */
void container_poly::clear() {
	for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
	max_radius=0;
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] fp a file handle to write to. */
template<class v_loop>
void container::draw_particles(v_loop &vl,FILE *fp) {
	double *pp;
	if(vl.start()) do {
		pp=p[vl.ijk]+3*vl.q;
		fprintf(fp,"%d %g %g %g\n",id[vl.ijk][vl.q],*pp,pp[1],pp[2]);
	} while(vl.inc());
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] fp a file handle to write to. */
template<class v_loop>
void container_poly::draw_particles(v_loop &vl,FILE *fp) {
	double *pp;
	if(vl.start()) do {
		pp=p[vl.ijk]+4*vl.q;
		fprintf(fp,"%d %g %g %g %g\n",id[vl.ijk][vl.q],*pp,pp[1],pp[2],pp[3]);
	} while(vl.inc());
}

/** Dumps all the particle positions in POV-Ray format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container::draw_particles_pov(v_loop &vl,FILE *fp) {
	double *pp;
	if(vl.start()) do {
		pp=p[vl.ijk]+3*vl.q;
		fprintf(fp,"// id %d\nsphere{<%g,%g,%g>,r}\n",id[vl.ijk][vl.q],*pp,pp[1],pp[2]);
	} while(vl.inc());
}

/** Dumps all the particle positions in POV-Ray format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container_poly::draw_particles_pov(v_loop &vl,FILE *fp) {
	double *pp;
	if(vl.start()) do {
		pp=p[vl.ijk]+4*vl.q;
		fprintf(fp,"// id %d\nsphere{<%g,%g,%g>,%g}\n",id[vl.ijk][vl.q],*pp,pp[1],pp[2],pp[3]);
	} while(vl.inc());
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in gnuplot format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container::draw_cells_gnuplot(v_loop &vl,FILE *fp) {
	voronoicell c;double *pp;
	if(vl.start()) do if(compute_cell(c,vl)) {
		pp=p[vl.ijk]+3*vl.q;
		c.draw_gnuplot(*pp,pp[1],pp[2],fp);
	} while(vl.inc());
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in gnuplot format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container_poly::draw_cells_gnuplot(v_loop &vl,FILE *fp) {
	voronoicell c;double *pp;
	if(vl.start()) do if(compute_cell(c,vl)) {
		pp=p[vl.ijk]+4*vl.q;
		c.draw_gnuplot(*pp,pp[1],pp[2],fp);
	} while(vl.inc());
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in POV-Ray format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container::draw_cells_pov(v_loop &vl,FILE *fp) {
	voronoicell c;double *pp;
	if(vl.start()) do if(compute_cell(c,vl)) {
		fprintf(fp,"// cell %d\n",id[vl.ijk][vl.q]);
		pp=p[vl.ijk]+3*vl.q;
		c.draw_pov(*pp,pp[1],pp[2],fp);
	} while(vl.inc());
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in POV-Ray format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container_poly::draw_cells_pov(v_loop &vl,FILE *fp) {
	voronoicell c;double *pp;
	if(vl.start()) do if(compute_cell(c,vl)) {
		fprintf(fp,"// cell %d\n",id[vl.ijk][vl.q]);
		pp=p[vl.ijk]+4*vl.q;
		c.draw_pov(*pp,pp[1],pp[2],fp);
	} while(vl.inc());
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in POV-Ray format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container::print_custom(v_loop &vl,const char *format,FILE *fp) {
	int ijk,q;double *pp;
	if(contains_neighbor(format)) {
		voronoicell_neighbor c;
		if(vl.start()) do if(compute_cell(c,vl)) {
			ijk=vl.ijk;q=vl.q;pp=p[ijk]+3*q;
			c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],default_radius,fp);
		} while(vl.inc());
	} else {
		voronoicell c;
		if(vl.start()) do if(compute_cell(c,vl)) {
			ijk=vl.ijk;q=vl.q;pp=p[ijk]+3*q;
			c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],default_radius,fp);
		} while(vl.inc());
	}
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in POV-Ray format.
 * \param[in] fp a file handle to write to. */ 
template<class v_loop>
void container_poly::print_custom(v_loop &vl,const char *format,FILE *fp) {
	int ijk,q;double *pp;
	if(contains_neighbor(format)) {
		voronoicell_neighbor c;
		if(vl.start()) do if(compute_cell(c,vl)) {
			ijk=vl.ijk;q=vl.q;pp=p[ijk]+4*q;
			c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],pp[3],fp);
		} while(vl.inc());
	} else {
		voronoicell c;
		if(vl.start()) do if(compute_cell(c,vl)) {
			ijk=vl.ijk;q=vl.q;pp=p[ijk]+4*q;
			c.output_custom(format,id[ijk][q],*pp,pp[1],pp[2],pp[3],fp);
		} while(vl.inc());
	}
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
void container::compute_all_cells() {
	voronoicell c;
	v_loop_all vl(*this);
	if(vl.start()) do compute_cell(c,vl);
	while(vl.inc());
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
void container_poly::compute_all_cells() {
	voronoicell c;
	v_loop_all vl(*this);
	if(vl.start()) do compute_cell(c,vl);while(vl.inc());
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container::sum_cell_volumes() {
	voronoicell c;
	double vol=0;
	v_loop_all vl(*this);
	if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
	return vol;	
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_poly::sum_cell_volumes() {
	voronoicell c;
	double vol=0;
	v_loop_all vl(*this);
	if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
	return vol;	
}

/** This function tests to see if a given vector lies within the container
 * bounds and any walls.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
bool container_base::point_inside(double x,double y,double z) {
	if(x<ax||x>bx||y<ay||y>by||z<az||z>bz) return false;
	return point_inside_walls(x,y,z);
}

/** The wall_list constructor sets up an array of pointers to wall classes. */
wall_list::wall_list() : walls(new wall*[init_wall_size]), wep(walls), wel(walls+init_wall_size),
	current_wall_size(init_wall_size) {}

/** The wall_list destructor frees the array of pointers to the wall classes.
 */
wall_list::~wall_list() {
	delete [] walls;
}

/** Adds all of the walls on another wall_list to this class. 
 * \param[in] wl a reference to the wall class. */
void wall_list::add_wall(wall_list &wl) {
	for(wall **wp=wl.walls;wp<wl.wep;wp++) add_wall(*wp);
}

/** Deallocates all of the wall classes pointed to by the wall_list. */
void wall_list::deallocate() {
	for(wall **wp=walls;wp<wep;wp++) delete *wp;
}

/** Increases the memory allocation for the walls array. */
void wall_list::increase_wall_memory() {
	current_wall_size<<=1;
	if(current_wall_size>max_wall_size)
		voropp_fatal_error("Wall memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
	wall **nwalls(new wall*[current_wall_size]),**nwp(nwalls),**wp(walls);
	while(wp<wep) *(nwp++)=*(wp++);
	delete [] walls;
	walls=nwalls;wel=walls+current_wall_size;wep=nwp;
}

// Explicit instantiation
template void container::draw_particles<v_loop_all>(v_loop_all&,FILE*);
template void container::draw_particles<v_loop_order>(v_loop_order&,FILE*);
template void container::draw_particles_pov<v_loop_all>(v_loop_all&,FILE*);
template void container::draw_particles_pov<v_loop_order>(v_loop_order&,FILE*);
template void container::draw_cells_pov<v_loop_all>(v_loop_all&,FILE*);
template void container::draw_cells_pov<v_loop_order>(v_loop_order&,FILE*);
template void container::draw_cells_gnuplot<v_loop_all>(v_loop_all&,FILE*);
template void container::draw_cells_gnuplot<v_loop_order>(v_loop_order&,FILE*);
template void container::print_custom<v_loop_all>(v_loop_all&,const char*,FILE*);
template void container::print_custom<v_loop_order>(v_loop_order&,const char*,FILE*);

template void container_poly::draw_particles<v_loop_all>(v_loop_all&,FILE*);
template void container_poly::draw_particles<v_loop_order>(v_loop_order&,FILE*);
template void container_poly::draw_particles_pov<v_loop_all>(v_loop_all&,FILE*);
template void container_poly::draw_particles_pov<v_loop_order>(v_loop_order&,FILE*);
template void container_poly::draw_cells_pov<v_loop_all>(v_loop_all&,FILE*);
template void container_poly::draw_cells_pov<v_loop_order>(v_loop_order&,FILE*);
template void container_poly::draw_cells_gnuplot<v_loop_all>(v_loop_all&,FILE*);
template void container_poly::draw_cells_gnuplot<v_loop_order>(v_loop_order&,FILE*);
template void container_poly::print_custom<v_loop_all>(v_loop_all&,const char*,FILE*);
template void container_poly::print_custom<v_loop_order>(v_loop_order&,const char*,FILE*);
