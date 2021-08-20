// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file container_3d.cc
 * \brief Function implementations for the container_3d and related classes. */

#include "container_3d.hh"

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                          coordinate directions.
 * \param[in] (x_prd_,y_prd_,z_prd_) flags setting whether the container is
 *                                   periodic in each coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_base_3d::container_base_3d(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
        int nx_,int ny_,int nz_,bool x_prd_,bool y_prd_,bool z_prd_,int init_mem,int ps_)
    : voro_base(nx_,ny_,nz_,(bx_-ax_)/nx_,(by_-ay_)/ny_,(bz_-az_)/nz_),
    ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
    max_len_sq((bx-ax)*(bx-ax)*(x_prd_?0.25:1)+(by-ay)*(by-ay)*(y_prd_?0.25:1)
          +(bz-az)*(bz-az)*(z_prd_?0.25:1)),
    x_prd(x_prd_), y_prd(y_prd_), z_prd(z_prd_),
    id(new int*[nxyz]), p(new double*[nxyz]), co(new int[nxyz]), mem(new int[nxyz]), ps(ps_) {

    int l;
    for(l=0;l<nxyz;l++) co[l]=0;
    for(l=0;l<nxyz;l++) mem[l]=init_mem;
    for(l=0;l<nxyz;l++) id[l]=new int[init_mem];
    for(l=0;l<nxyz;l++) p[l]=new double[ps*init_mem];
}

/** The container destructor frees the dynamically allocated memory. */
container_base_3d::~container_base_3d() {
    int l;
    for(l=0;l<nxyz;l++) delete [] p[l];
    for(l=0;l<nxyz;l++) delete [] id[l];
    delete [] id;
    delete [] p;
    delete [] co;
    delete [] mem;
}

/** The class constructor sets up the geometry of the container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                          coordinate directions.
 * \param[in] (x_prd_,y_prd_,z_prd_) flags setting whether the container is
 *                                   periodic in each coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block. */
container_3d::container_3d(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
    int nx_,int ny_,int nz_,bool x_prd_,bool y_prd_,bool z_prd_,int init_mem)
    : container_base_3d(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,x_prd_,y_prd_,z_prd_,init_mem,3),
    vc(*this,x_prd_?2*nx_+1:nx_,y_prd_?2*ny_+1:ny_,z_prd_?2*nz_+1:nz_) {}

/** The class constructor sets up the geometry of container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                          coordinate directions.
 * \param[in] (x_prd_,y_prd_,z_prd_) flags setting whether the container is
 *                                   periodic in each coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block. */
container_poly_3d::container_poly_3d(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
    int nx_,int ny_,int nz_,bool x_prd_,bool y_prd_,bool z_prd_,int init_mem)
    : container_base_3d(ax_,bx_,ay_,by_,az_,bz_,nx_,ny_,nz_,x_prd_,y_prd_,z_prd_,init_mem,4),
    vc(*this,x_prd_?2*nx_+1:nx_,y_prd_?2*ny_+1:ny_,z_prd_?2*nz_+1:nz_) {ppr=p;}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_3d::put(int n,double x,double y,double z) {
    int ijk;
    if(put_locate_block(ijk,x,y,z)) {
        id[ijk][co[ijk]]=n;
        double *pp=p[ijk]+3*co[ijk]++;
        *(pp++)=x;*(pp++)=y;*pp=z;
    }
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly_3d::put(int n,double x,double y,double z,double r) {
    int ijk;
    if(put_locate_block(ijk,x,y,z)) {
        id[ijk][co[ijk]]=n;
        double *pp=p[ijk]+4*co[ijk]++;
        *(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
        if(max_radius<r) max_radius=r;
    }
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_3d::put(particle_order &vo,int n,double x,double y,double z) {
    int ijk;
    if(put_locate_block(ijk,x,y,z)) {
        id[ijk][co[ijk]]=n;
        vo.add(ijk,co[ijk]);
        double *pp=p[ijk]+3*co[ijk]++;
        *(pp++)=x;*(pp++)=y;*pp=z;
    }
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly_3d::put(particle_order &vo,int n,double x,double y,double z,double r) {
    int ijk;
    if(put_locate_block(ijk,x,y,z)) {
        id[ijk][co[ijk]]=n;
        vo.add(ijk,co[ijk]);
        double *pp=p[ijk]+4*co[ijk]++;
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
bool container_base_3d::put_locate_block(int &ijk,double &x,double &y,double &z) {
    if(put_remap(ijk,x,y,z)) {
        if(co[ijk]==mem[ijk]) add_particle_memory(ijk);
        return true;
    }
#if VOROPP_REPORT_OUT_OF_BOUNDS ==1
    fprintf(stderr,"Out of bounds: (x,y,z)=(%g,%g,%g)\n",x,y,z);
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
inline bool container_base_3d::put_remap(int &ijk,double &x,double &y,double &z) {
    int l;

    ijk=step_int((x-ax)*xsp);
    if(x_prd) {l=step_mod(ijk,nx);x+=boxx*(l-ijk);ijk=l;}
    else if(ijk<0||ijk>=nx) return false;

    int j=step_int((y-ay)*ysp);
    if(y_prd) {l=step_mod(j,ny);y+=boxy*(l-j);j=l;}
    else if(j<0||j>=ny) return false;

    int k=step_int((z-az)*zsp);
    if(z_prd) {l=step_mod(k,nz);z+=boxz*(l-k);k=l;}
    else if(k<0||k>=nz) return false;

    ijk+=nx*j+nxy*k;
    return true;
}

/** Takes a position vector and attempts to remap it into the primary domain.
 * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
 *                       with (0,0,0) corresponding to the primary domain.
 * \param[out] (ci,cj,ck) the index of the block that the position vector is
 *                        within, once it has been remapped.
 * \param[in,out] (x,y,z) the position vector to consider, which is remapped
 *                        into the primary domain during the routine.
 * \param[out] ijk the block index that the vector is within.
 * \return True if the particle is within the container or can be remapped into
 * it, false if it lies outside of the container bounds. */
inline bool container_base_3d::remap(int &ai,int &aj,int &ak,int &ci,int &cj,int &ck,double &x,double &y,double &z,int &ijk) {
    ci=step_int((x-ax)*xsp);
    if(ci<0||ci>=nx) {
        if(x_prd) {ai=step_div(ci,nx);x-=ai*(bx-ax);ci-=ai*nx;}
        else return false;
    } else ai=0;

    cj=step_int((y-ay)*ysp);
    if(cj<0||cj>=ny) {
        if(y_prd) {aj=step_div(cj,ny);y-=aj*(by-ay);cj-=aj*ny;}
        else return false;
    } else aj=0;

    ck=step_int((z-az)*zsp);
    if(ck<0||ck>=nz) {
        if(z_prd) {ak=step_div(ck,nz);z-=ak*(bz-az);ck-=ak*nz;}
        else return false;
    } else ak=0;

    ijk=ci+nx*cj+nxy*ck;
    return true;
}

/** Takes a vector and finds the particle whose Voronoi cell contains that
 * vector. This is equivalent to finding the particle which is nearest to the
 * vector. Additional wall classes are not considered by this routine.
 * \param[in] (x,y,z) the vector to test.
 * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
 *                        contains the vector. If the container is periodic,
 *                        this may point to a particle in a periodic image of
 *                        the primary domain.
 * \param[out] pid the ID of the particle.
 * \return True if a particle was found. If the container has no particles,
 * then the search will not find a Voronoi cell and false is returned. */
bool container_3d::find_voronoi_cell(double x,double y,double z,double &rx,double &ry,double &rz,int &pid) {
    int ai,aj,ak,ci,cj,ck,ijk;
    particle_record w;
    double mrs;

    // If the given vector lies outside the domain, but the container
    // is periodic, then remap it back into the domain
    if(!remap(ai,aj,ak,ci,cj,ck,x,y,z,ijk)) return false;
    vc.find_voronoi_cell(x,y,z,ci,cj,ck,ijk,w,mrs);

    if(w.ijk!=-1) {

        // Assemble the position vector of the particle to be returned,
        // applying a periodic remapping if necessary
        if(x_prd) {ci+=w.di;if(ci<0||ci>=nx) ai+=step_div(ci,nx);}
        if(y_prd) {cj+=w.dj;if(cj<0||cj>=ny) aj+=step_div(cj,ny);}
        if(z_prd) {ck+=w.dk;if(ck<0||ck>=nz) ak+=step_div(ck,nz);}
        rx=p[w.ijk][3*w.l]+ai*(bx-ax);
        ry=p[w.ijk][3*w.l+1]+aj*(by-ay);
        rz=p[w.ijk][3*w.l+2]+ak*(bz-az);
        pid=id[w.ijk][w.l];
        return true;
    }

    // If no particle is found then just return false
    return false;
}

/** Takes a vector and finds the particle whose Voronoi cell contains that
 * vector. Additional wall classes are not considered by this routine.
 * \param[in] (x,y,z) the vector to test.
 * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
 *                        contains the vector. If the container is periodic,
 *                        this may point to a particle in a periodic image of
 *                        the primary domain.
 * \param[out] pid the ID of the particle.
 * \return True if a particle was found. If the container has no particles,
 * then the search will not find a Voronoi cell and false is returned. */
bool container_poly_3d::find_voronoi_cell(double x,double y,double z,double &rx,double &ry,double &rz,int &pid) {
    int ai,aj,ak,ci,cj,ck,ijk;
    particle_record w;
    double mrs;

    // If the given vector lies outside the domain, but the container is
    // periodic, then remap it back into the domain
    if(!remap(ai,aj,ak,ci,cj,ck,x,y,z,ijk)) return false;
    vc.find_voronoi_cell(x,y,z,ci,cj,ck,ijk,w,mrs);

    if(w.ijk!=-1) {

        // Assemble the position vector of the particle to be returned,
        // applying a periodic remapping if necessary
        if(x_prd) {ci+=w.di;if(ci<0||ci>=nx) ai+=step_div(ci,nx);}
        if(y_prd) {cj+=w.dj;if(cj<0||cj>=ny) aj+=step_div(cj,ny);}
        if(z_prd) {ck+=w.dk;if(ck<0||ck>=nz) ak+=step_div(ck,nz);}
        rx=p[w.ijk][4*w.l]+ai*(bx-ax);
        ry=p[w.ijk][4*w.l+1]+aj*(by-ay);
        rz=p[w.ijk][4*w.l+2]+ak*(bz-az);
        pid=id[w.ijk][w.l];
        return true;
    }

    // If no particle is found then just return false
    return false;
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_base_3d::add_particle_memory(int i) {
    int l,nmem=mem[i]<<1;

    // Carry out a check on the memory allocation size, and print a status
    // message if requested
    if(nmem>max_particle_memory)
        voro_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=3
    fprintf(stderr,"Particle memory in region %d scaled up to %d\n",i,nmem);
#endif

    // Allocate new memory and copy in the contents of the old arrays
    int *idp=new int[nmem];
    for(l=0;l<co[i];l++) idp[l]=id[i][l];
    double *pp=new double[ps*nmem];
    for(l=0;l<ps*co[i];l++) pp[l]=p[i][l];

    // Update pointers and delete old arrays
    mem[i]=nmem;
    delete [] id[i];id[i]=idp;
    delete [] p[i];p[i]=pp;
}

/** Import a list of particles from an open file stream into the container.
 * Entries of four numbers (Particle ID, x position, y position, z position)
 * are searched for. If the file cannot be successfully read, then the routine
 * causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_3d::import(FILE *fp) {
    int i,j;
    double x,y,z;
    while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(i,x,y,z);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position) are searched for. If the file cannot be
 * successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_3d::import(particle_order &vo,FILE *fp) {
    int i,j;
    double x,y,z;
    while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(vo,i,x,y,z);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream into the container.
 * Entries of five numbers (Particle ID, x position, y position, z position,
 * radius) are searched for. If the file cannot be successfully read, then the
 * routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_poly_3d::import(FILE *fp) {
    int i,j;
    double x,y,z,r;
    while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(i,x,y,z,r);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position, radius) are searched for. If the file
 * cannot be successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_poly_3d::import(particle_order &vo,FILE *fp) {
    int i,j;
    double x,y,z,r;
    while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(vo,i,x,y,z,r);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Outputs the a list of all the container regions along with the number of
 * particles stored within each. */
void container_base_3d::region_count() {
    int i,j,k,*cop=co;
    for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
        printf("Region (%d,%d,%d): %d particles\n",i,j,k,*(cop++));
}

/** Clears a container of particles. */
void container_3d::clear() {
    for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
}

/** Clears a container of particles, also clearing resetting the maximum radius
 * to zero. */
void container_poly_3d::clear() {
    for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
    max_radius=0;
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_3d::print_custom(const char *format,FILE *fp) {
    c_loop_all vl(*this);
    print_custom(vl,format,fp);
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_poly_3d::print_custom(const char *format,FILE *fp) {
    c_loop_all vl(*this);
    print_custom(vl,format,fp);
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] filename the name of the file to write to. */
void container_3d::print_custom(const char *format,const char *filename) {
    FILE *fp=safe_fopen(filename,"w");
    print_custom(format,fp);
    fclose(fp);
}

/** Computes all the Voronoi cells and saves customized information about them
 * \param[in] format the custom output string to use.
 * \param[in] filename the name of the file to write to. */
void container_poly_3d::print_custom(const char *format,const char *filename) {
    FILE *fp=safe_fopen(filename,"w");
    print_custom(format,fp);
    fclose(fp);
}

/** Computes all of the Voronoi cells in the container, but does nothing with
 * the output. It is useful for measuring the pure computation time of the
 * Voronoi algorithm, without any additional calculations such as volume
 * evaluation or cell output. */
void container_3d::compute_all_cells() {
    voronoicell c(*this);
    c_loop_all vl(*this);
    if(vl.start()) do compute_cell(c,vl);
    while(vl.inc());
}

/** Computes all of the Voronoi cells in the container, but does nothing with
 * the output. It is useful for measuring the pure computation time of the
 * Voronoi algorithm, without any additional calculations such as volume
 * evaluation or cell output. */
void container_poly_3d::compute_all_cells() {
    voronoicell c(*this);
    c_loop_all vl(*this);
    if(vl.start()) do compute_cell(c,vl);while(vl.inc());
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_3d::sum_cell_volumes() {
    voronoicell c(*this);
    double vol=0;
    c_loop_all vl(*this);
    if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
    return vol;
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_poly_3d::sum_cell_volumes() {
    voronoicell c(*this);
    double vol=0;
    c_loop_all vl(*this);
    if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
    return vol;
}

/** This function tests to see if a given vector lies within the container
 * bounds and any walls.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
bool container_base_3d::point_inside(double x,double y,double z) {
    if(x<ax||x>bx||y<ay||y>by||z<az||z>bz) return false;
    return point_inside_walls(x,y,z);
}

/** Draws an outline of the domain in Gnuplot format.
 * \param[in] fp the file handle to write to. */
void container_base_3d::draw_domain_gnuplot(FILE *fp) {
    fprintf(fp,"%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n"
               "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n"
               "%g %g %g\n\n%g %g %g\n%g %g %g\n\n"
               "%g %g %g\n%g %g %g\n\n%g %g %g\n%g %g %g\n\n",
               ax,ay,az,bx,ay,az,bx,by,az,ax,by,az,
               ax,by,bz,bx,by,bz,bx,ay,bz,ax,ay,bz,
               ax,by,bz,ax,ay,az,ax,ay,bz,
               bx,ay,az,bx,ay,bz,bx,by,az,bx,by,bz);
}

/** Draws an outline of the domain in POV-Ray format.
 * \param[in] fp the file handle to write to. */
void container_base_3d::draw_domain_pov(FILE *fp) {
    fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
               "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
               "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
               "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
               "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
               ax,ay,az,bx,ay,az,ax,by,az,bx,by,az,
               ax,by,bz,bx,by,bz,ax,ay,bz,bx,ay,bz,
               ax,ay,az,ax,by,az,bx,ay,az,bx,by,az,
               bx,ay,bz,bx,by,bz,ax,ay,bz,ax,by,bz,
               ax,ay,az,ax,ay,bz,bx,ay,az,bx,ay,bz,
               bx,by,az,bx,by,bz,ax,by,az,ax,by,bz,
               ax,ay,az,bx,ay,az,ax,by,az,bx,by,az,
               ax,ay,bz,bx,ay,bz,ax,by,bz,bx,by,bz);
}

}
