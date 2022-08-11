// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file container_2d.cc
 * \brief Function implementations for the container_2d and related classes. */

#include <cstring>

#include "container_2d.hh"
#include "iter_2d.hh"

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (nx_,ny_) the number of grid blocks in each of the three
 *                      coordinate directions.
 * \param[in] (x_prd_,y_prd_) flags setting whether the container is
 *                      periodic in each coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle.
 * \param[in] nt_ the maximum number of threads that will be used for Voronoi
 *                computations. */
container_base_2d::container_base_2d(double ax_,double bx_,double ay_,double by_,
        int nx_,int ny_,bool x_prd_,bool y_prd_,int init_mem,int ps_,int nt_)
    : voro_base_2d(nx_,ny_,(bx_-ax_)/nx_,(by_-ay_)/ny_),
    ax(ax_), bx(bx_), ay(ay_), by(by_), x_prd(x_prd_), y_prd(y_prd_),
    id(new int*[nxy]), p(new double*[nxy]), co(new int[nxy]), mem(new int[nxy]),
    ps(ps_), nt(nt_), oflow_co(0), oflow_mem(init_overflow_size),
    ij_m_id_oflow(new int[3*oflow_mem]), p_oflow(new double[ps*oflow_mem]) {
    int l;
    for(l=0;l<nxy;l++) co[l]=0;
    for(l=0;l<nxy;l++) mem[l]=init_mem;
    for(l=0;l<nxy;l++) id[l]=new int[init_mem];
    for(l=0;l<nxy;l++) p[l]=new double[ps*init_mem];
}

/** The class destructor frees the dynamically allocated memory. */
container_base_2d::~container_base_2d() {

    // Delete the per-block arrays
    int l;
    for(l=nxy-1;l>=0;l--) delete [] p[l];
    for(l=nxy-1;l>=0;l--) delete [] id[l];

    // Delete the overflow arrays
    delete [] p_oflow;
    delete [] ij_m_id_oflow;

    // Delete the block arrays
    delete [] id;
    delete [] p;
    delete [] co;
    delete [] mem;
}

/** The class constructor sets up the geometry of container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (nx_,ny_) the number of grid blocks in each of the three
 *                      coordinate directions.
 * \param[in] (x_prd_,y_prd_) flags setting whether the container is periodic
 *                            in each coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] nt_ the maximum number of threads that will be used for Voronoi
 *                computations. */
container_2d::container_2d(double ax_,double bx_,double ay_,double by_,
    int nx_,int ny_,bool x_prd_,bool y_prd_,int init_mem,int nt_)
    : container_base_2d(ax_,bx_,ay_,by_,nx_,ny_,x_prd_,y_prd_,init_mem,2,nt_),
    vc(new voro_compute_2d<container_2d>*[nt]) {

    // Allocate as many Voronoi computation objects as there are threads
#pragma omp parallel num_threads(nt)
    {
        vc[t_num()]= new voro_compute_2d<container_2d>(*this,x_prd_?2*nx_+1:nx_,y_prd_?2*ny_+1:ny_);
    }
}

/** The class destructor frees the dynamically allocated memory. */
container_2d::~container_2d(){

    // Delete the Voronoi computation objects
    for(int l=0;l<nt;l++) delete vc[l];

    // Delete the Voronoi computation pointer array
    delete [] vc;
}

/** Changes the maximum number of threads that can be used in multithreaded
 * computations.
 * \param[in] nt_ the new maximum number of threads. */
void container_2d::change_number_thread(int nt_){

    // Delete the previous Voronoi computation objects
    for(int l=0;l<nt;l++) delete vc[l];
    delete [] vc;

    // Allocate the new Voronoi computation objects
    nt=nt_;
    vc=new voro_compute_2d<container_2d>*[nt];
#pragma omp parallel num_threads(nt)
    {
        vc[t_num()]=new voro_compute_2d<container_2d>(*this,x_prd?2*nx+1:nx,y_prd?2*ny+1:ny);
    }
}

/** The class constructor sets up the geometry of container.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (nx_,ny_) the number of grid blocks in each of the three
 *                      coordinate directions.
 * \param[in] (x_prd_,y_prd_) flags setting whether the container is periodic
 *                            in each coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block. */
container_poly_2d::container_poly_2d(double ax_,double bx_,double ay_,double by_,
    int nx_,int ny_,bool x_prd_,bool y_prd_,int init_mem,int nt_)
    : container_base_2d(ax_,bx_,ay_,by_,nx_,ny_,x_prd_,y_prd_,init_mem,3,nt_),
    vc(new voro_compute_2d<container_poly_2d>*[nt]), max_r(new double[nt]) {
    for(int j=0;j<nt;j++) max_r[j]=0.;
#pragma omp parallel num_threads(nt)
    {
        vc[t_num()]=new voro_compute_2d<container_poly_2d>(*this,x_prd_?2*nx_+1:nx_,y_prd_?2*ny_+1:ny_);
    }
    ppr=p;
}

/** The class destructor frees the dynamically allocated memory. */
container_poly_2d::~container_poly_2d() {

    // Delete the Voronoi computation objects
    for(int l=0;l<nt;l++) delete vc[l];

    // Delete other arrays
    delete [] vc;
    delete [] max_r;
}

/** Changes the maximum number of threads that can be used in multithreaded
 * computations.
 * \param[in] nt_ the new maximum number of threads. */
void container_poly_2d::change_number_thread(int nt_) {

    // Delete the previous Voronoi computation objects and maximum radius array
    for(int l=0;l<nt;l++) delete vc[l];
    delete [] vc;
    delete [] max_r;

    // Allocate the new Voronoi computation objects and maximum radius array
    nt=nt_;
    max_r=new double[nt];
    vc=new voro_compute_2d<container_poly_2d>*[nt];
    #pragma omp parallel num_threads(nt)
    {
        vc[t_num()]=new voro_compute_2d<container_poly_2d>(*this,x_prd?2*nx+1:nx,y_prd?2*ny+1:ny);
    }
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle. */
void container_2d::put(int n,double x,double y) {
    int ij;
    if(put_locate_block(ij,x,y)) {
        id[ij][co[ij]]=n;
        double *pp=p[ij]+2*co[ij]++;
        *(pp++)=x;*pp=y;
    }
}

/** Put a particle into the correct region of the container.
 * \param[in] i the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle. */
void container_2d::put_parallel(int i,double x,double y) {
    int ij;

    // Locate the particle that the block is in
    if(put_remap(ij,x,y)) {

        // Find the unique available slot number in the block, using an atomic
        // increment of the counter
        int m;
#pragma omp atomic capture
        m=co[ij]++;

        // If the slot is within the available allocated memory, then add it
        // directly
        if(m<mem[ij]){
            id[ij][m]=i;
            double *pp=p[ij]+2*m;
            *pp=x;pp[1]=y;
        } else {

            // Otherwise, store it into the overflow array to reconcile later.
            // This routine can only be accessed by one thread at a time, in
            // case the overflow buffer needs to be extended.
#pragma omp critical
            {
                if(oflow_co>=oflow_mem) add_overflow_memory();
                int *idp=ij_m_id_oflow+3*oflow_co;
                *idp=ij;idp[1]=m;idp[2]=i;
                double *pp=p_oflow+2*oflow_co++;
                *pp=x;pp[1]=y;
            }
        }
    }
}

/** Increases the size of the overflow buffer for multithreaded insertion of particles
 * into the container. */
void container_base_2d::add_overflow_memory() {

    // Calculate the size of the extended array and check it is within the
    // limit
    oflow_mem<<=1;
    if(oflow_mem>max_overflow_size)
        voro_fatal_error("Maximum overflow memory size exceeded",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=3
    fprintf(stderr,"Overflow memory scaled up to %d\n",oflow_mem);
#endif

    // Copy integer information from old array to new
    int *ij_m_id_new=new int[3*oflow_mem];
    memcpy(ij_m_id_new,ij_m_id_oflow,3*sizeof(int)*oflow_co);
    delete [] ij_m_id_oflow;
    ij_m_id_oflow=ij_m_id_new;

    // Copy floating point information from old array to new
    double *p_new=new double[ps*oflow_mem];
    memcpy(p_new,p_oflow,sizeof(double)*ps*oflow_co);
    delete [] p_oflow;
    p_oflow=p_new;
}

/** Adds an array of particle positions to the container using multithreaded
 * insertion.
 * \param[in] pt_list a pointer to the array of positions, stored as (x,y)
 *                    doublets.
 * \param[in] num the number of particles.
 * \param[in] nt_ the number of threads to use. */
void container_2d::add_parallel(double *pt_list,int num,int nt_) {
#pragma omp parallel for num_threads(nt_)
    for(int i=0;i<num;i++){
        double *pp=pt_list+2*i;
        put_parallel(i,*pp,pp[1]);
    }
}

/** Adds the particles stored in the overflow buffer to the container. */
void container_2d::put_reconcile_overflow() {

    // Consider each of the particles in the buffer
    double *op=p_oflow;
    for(int *idp=ij_m_id_oflow;idp<=ij_m_id_oflow+3*oflow_co;) {

        // Add particle memory if needed to store this particle
        int ij=*(idp++),m=*(idp++);
        if(m>=mem[ij]) add_particle_memory(ij,m);

        // Store the particle information
        id[ij][m]=*(idp++);
        double *pp=p[ij]+2*m;
        *pp=*(op++);pp[1]=*(op++);
    }

    // All particles in the overflow buffer have been considered, so set the
    // overflow counter to zero
    oflow_co=0;
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly_2d::put(int n,double x,double y,double r) {
    int ij;
    if(put_locate_block(ij,x,y)) {
        id[ij][co[ij]]=n;
        double *pp=p[ij]+3*co[ij]++;
        *(pp++)=x;*(pp++)=y;*pp=r;
        if(max_radius<r) max_radius=r;
    }
}

/** Put a particle into the correct region of the container.
 * \param[in] i the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly_2d::put_parallel(int i,double x,double y,double r) {
    int ij;

    // Locate the particle that the block is in
    if(put_remap(ij,x,y)) {

        // Update the per-thread maximum radius
        int tn=t_num();
        if(max_r[tn]<r) {max_r[tn]=r;}

        // Find the unique available slot number in the block, using an atomic
        // increment of the counter
        int m;
#pragma omp atomic capture
        m=co[ij]++;

        // If the slot is within the available allocated memory, then add it
        // directly
        if(m<mem[ij]){
            id[ij][m]=i;
            double *pp=p[ij]+3*m;
            *pp=x;pp[1]=y;pp[2]=r;
        } else {

            // Otherwise, store it into the overflow array to reconcile later.
            // This routine can only be accessed by one thread at a time, in
            // case the overflow buffer needs to be extended.
#pragma omp critical
            {
                if(oflow_co>=oflow_mem) add_overflow_memory();
                int *idp=ij_m_id_oflow+3*oflow_co;
                *idp=ij;idp[1]=m;idp[2]=i;
                double *pp=p_oflow+3*oflow_co++;
                *pp=x;pp[1]=y;pp[2]=r;
            }
        }
    }
}

/** Adds an array of particle positions to the container using multithreaded
 * insertion.
 * \param[in] pt_list a pointer to the array of positions, stored as (x,y,r)
 *                    triplets.
 * \param[in] num the number of particles.
 * \param[in] nt_ the number of threads to use. */
void container_poly_2d::add_parallel(double *pt_list,int num,int nt_){
#pragma omp parallel for num_threads(nt_)
    for(int i=0;i<num;i++) {
        double *pp=pt_list+3*i;
        put_parallel(i,*pp,pp[1],pp[2]);
    }
}

/** Adds the particles stored in the overflow buffer to the container. */
void container_poly_2d::put_reconcile_overflow() {

    // Compute the global maximum radius using the per-thread values
    for(int i=0;i<nt;i++) {
        if(max_radius<max_r[i]) max_radius=max_r[i];
        max_r[i]=0.;
    }

    // Consider each of the particles in the buffer
    double *op=p_oflow;
    for(int *idp=ij_m_id_oflow;idp<=ij_m_id_oflow+3*oflow_co;) {

        // Add particle memory if needed to store this particle
        int ij=*(idp++),m=*(idp++);
        if(m>=mem[ij]) add_particle_memory(ij,m);

        // Store the particle information
        id[ij][m]=*(idp++);
        double *pp=p[ij]+3*m;
        *pp=*(op++);pp[1]=*(op++);pp[2]=*(op++);
    }

    // All particles in the overflow buffer have been considered, so set the
    // overflow counter to zero
    oflow_co=0;
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle. */
void container_2d::put(particle_order &vo,int n,double x,double y) {
    int ij;
    if(put_locate_block(ij,x,y)) {
        id[ij][co[ij]]=n;
        vo.add(ij,co[ij]);
        double *pp=p[ij]+2*co[ij]++;
        *(pp++)=x;*pp=y;
    }
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_poly_2d::put(particle_order &vo,int n,double x,double y,double r) {
    int ij;
    if(put_locate_block(ij,x,y)) {
        id[ij][co[ij]]=n;
        vo.add(ij,co[ij]);
        double *pp=p[ij]+3*co[ij]++;
        *(pp++)=x;*(pp++)=y;*pp=r;
        if(max_radius<r) max_radius=r;
    }
}

/** This routine takes a particle position vector, tries to remap it into the
 * primary domain. If successful, it computes the region into which it can be
 * stored and checks that there is enough memory within this region to store
 * it.
 * \param[out] ij the region index.
 * \param[in,out] (x,y) the particle position, remapped into the primary
 *                      domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base_2d::put_locate_block(int &ij,double &x,double &y) {
    if(put_remap(ij,x,y)) {
        if(co[ij]==mem[ij]) add_particle_memory(ij,co[ij]);
        return true;
    }
#if VOROPP_REPORT_OUT_OF_BOUNDS ==1
    fprintf(stderr,"Out of bounds: (x,y)=(%g,%g)\n",x,y);
#endif
    return false;
}

/** Takes a particle position vector and computes the region index into which
 * it should be stored. If the container is periodic, then the routine also
 * maps the particle position to ensure it is in the primary domain. If the
 * container is not periodic, the routine bails out.
 * \param[out] ij the region index.
 * \param[in,out] (x,y) the particle position, remapped into the primary domain
 *                      if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base_2d::put_remap(int &ij,double &x,double &y) {
    int l;

    ij=step_int((x-ax)*xsp);
    if(x_prd) {l=step_mod(ij,nx);x+=boxx*(l-ij);ij=l;}
    else if(ij<0||ij>=nx) return false;

    int j=step_int((y-ay)*ysp);
    if(y_prd) {l=step_mod(j,ny);y+=boxy*(l-j);j=l;}
    else if(j<0||j>=ny) return false;

    ij+=nx*j;
    return true;
}

/** Takes a position vector and attempts to remap it into the primary domain.
 * \param[out] (ai,aj) the periodic image displacement that the vector is in,
 *                     with (0,0,0) corresponding to the primary domain.
 * \param[out] (ci,cj) the index of the block that the position vector is
 *                     within, once it has been remapped.
 * \param[in,out] (x,y) the position vector to consider, which is remapped into
 *                      the primary domain during the routine.
 * \param[out] ij the block index that the vector is within.
 * \return True if the particle is within the container or can be remapped into
 * it, false if it lies outside of the container bounds. */
inline bool container_base_2d::remap(int &ai,int &aj,int &ci,int &cj,double &x,double &y,int &ij) {
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

    ij=ci+nx*cj;
    return true;
}

/** Takes a vector and finds the particle whose Voronoi cell contains that
 * vector. This is equivalent to finding the particle which is nearest to the
 * vector. Additional wall classes are not considered by this routine.
 * \param[in] (x,y) the vector to test.
 * \param[out] (rx,ry) the position of the particle whose Voronoi cell contains
 *                     the vector. If the container is periodic, this may point
 *                     to a particle in a periodic image of the primary domain.
 * \param[out] pid the ID of the particle.
 * \return True if a particle was found. If the container has no particles,
 * then the search will not find a Voronoi cell and false is returned. */
bool container_2d::find_voronoi_cell(double x,double y,double &rx,double &ry,int &pid) {
    int ai,aj,ci,cj,ij;
    particle_record_2d w;
    double mrs;

    // If the given vector lies outside the domain, but the container is
    // periodic, then remap it back into the domain
    if(!remap(ai,aj,ci,cj,x,y,ij)) return false;
    const int tn=t_num();
    vc[tn]->find_voronoi_cell(x,y,ci,cj,ij,w,mrs);

    if(w.ij!=-1) {

        // Assemble the position vector of the particle to be returned,
        // applying a periodic remapping if necessary
        if(x_prd) {ci+=w.di;if(ci<0||ci>=nx) ai+=step_div(ci,nx);}
        if(y_prd) {cj+=w.dj;if(cj<0||cj>=ny) aj+=step_div(cj,ny);}
        rx=p[w.ij][2*w.l]+ai*(bx-ax);
        ry=p[w.ij][2*w.l+1]+aj*(by-ay);
        pid=id[w.ij][w.l];
        return true;
    }

    // If no particle is found then just return false
    return false;
}

/** Takes a vector and finds the particle whose Voronoi cell contains that
 * vector. Additional wall classes are not considered by this routine.
 * \param[in] (x,y) the vector to test.
 * \param[out] (rx,ry) the position of the particle whose Voronoi cell contains
 *                     the vector. If the container is periodic, this may point
 *                     to a particle in a periodic image of the primary domain.
 * \param[out] pid the ID of the particle.
 * \return True if a particle was found. If the container has no particles,
 * then the search will not find a Voronoi cell and false is returned. */
bool container_poly_2d::find_voronoi_cell(double x,double y,double &rx,double &ry,int &pid) {
    int ai,aj,ci,cj,ij;
    particle_record_2d w;
    double mrs;

    // If the given vector lies outside the domain, but the container is
    // periodic, then remap it back into the domain
    if(!remap(ai,aj,ci,cj,x,y,ij)) return false;
    const int tn=t_num();
    vc[tn]->find_voronoi_cell(x,y,ci,cj,ij,w,mrs);

    if(w.ij!=-1) {

        // Assemble the position vector of the particle to be returned,
        // applying a periodic remapping if necessary
        if(x_prd) {ci+=w.di;if(ci<0||ci>=nx) ai+=step_div(ci,nx);}
        if(y_prd) {cj+=w.dj;if(cj<0||cj>=ny) aj+=step_div(cj,ny);}
        rx=p[w.ij][3*w.l]+ai*(bx-ax);
        ry=p[w.ij][3*w.l+1]+aj*(by-ay);
        pid=id[w.ij][w.l];
        return true;
    }

    // If no particle is found then just return false
    return false;
}

/** Increase memory for a particular region, within the parallel insertion
 * routines.
 * \param[in] i the index of the region to reallocate.
 * \param[in] m a minimum size for the reallocated region. */
void container_base_2d::add_particle_memory(int i,int m) {
    int omem=mem[i];
    do {mem[i]<<=1;} while(m>=mem[i]);

    // Check the memory allocation size and print a status message if requested
    if(mem[i]>max_particle_memory)
        voro_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=3
    fprintf(stderr,"Particle memory in region %d scaled up to %d\n",i,mem[i]);
#endif

    // Allocate new memory and copy in the contents of the old arrays
    int *idp=new int[mem[i]];
    memcpy(idp,id[i],sizeof(int)*omem);
    delete [] id[i];id[i]=idp;
    double *pp=new double[ps*mem[i]];
    memcpy(pp,p[i],ps*sizeof(double)*omem);
    delete [] p[i];p[i]=pp;
}

/** Import a list of particles from an open file stream into the container.
 * Entries of four numbers (Particle ID, x position, y position, z position)
 * are searched for. If the file cannot be successfully read, then the routine
 * causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_2d::import(FILE *fp) {
    int i,j;
    double x,y;
    while((j=fscanf(fp,"%d %lg %lg",&i,&x,&y))==3) put(i,x,y);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position) are searched for. If the file cannot be
 * successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_2d::import(particle_order &vo,FILE *fp) {
    int i,j;
    double x,y;
    while((j=fscanf(fp,"%d %lg %lg",&i,&x,&y))==3) put(vo,i,x,y);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream into the container.
 * Entries of five numbers (Particle ID, x position, y position, z position,
 * radius) are searched for. If the file cannot be successfully read, then the
 * routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_poly_2d::import(FILE *fp) {
    int i,j;
    double x,y,r;
    while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&r))==4) put(i,x,y,r);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position, radius) are searched for. If the file
 * cannot be successfully read, then the routine causes a fatal error.
 * \param[in,out] vo a reference to an ordering class to use.
 * \param[in] fp the file handle to read from. */
void container_poly_2d::import(particle_order &vo,FILE *fp) {
    int i,j;
    double x,y,r;
    while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&r))==4) put(vo,i,x,y,r);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Outputs the a list of all the container regions along with the number of
 * particles stored within each. */
void container_base_2d::region_count() {
    int i,j,*cop=co;
    for(j=0;j<ny;j++) for(i=0;i<nx;i++)
        printf("Region (%d,%d): %d particles\n",i,j,*(cop++));
}

/** Clears a container of particles. */
void container_2d::clear() {
    for(int *cop=co;cop<co+nxy;cop++) *cop=0;
}

/** Clears a container of particles, also clearing resetting the maximum radius
 * to zero. */
void container_poly_2d::clear() {
    for(int *cop=co;cop<co+nxy;cop++) *cop=0;
    max_radius=0;
}

/** This function tests to see if a given vector lies within the container
 * bounds and any walls.
 * \param[in] (x,y) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 * outside. */
bool container_base_2d::point_inside(double x,double y) {
    if(x<ax||x>bx||y<ay||y>by) return false;
    return point_inside_walls(x,y);
}

/** Draws an outline of the domain in Gnuplot format.
 * \param[in] fp the file handle to write to. */
void container_base_2d::draw_domain_gnuplot(FILE *fp) {
    fprintf(fp,"%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n",ax,ay,bx,ay,bx,by,ax,by,ax,ay);
}

/** Draws an outline of the domain in POV-Ray format.
 * \param[in] fp the file handle to write to. */
void container_base_2d::draw_domain_pov(FILE *fp) {
    fprintf(fp,"cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n"
               "cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n"
               "cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n"
               "cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n"
               "sphere{<%g,%g,0>,rr}\nsphere{<%g,%g,0>,rr}\n"
               "sphere{<%g,%g,0>,rr}\nsphere{<%g,%g,0>,rr}\n",
               ax,ay,bx,ay,ax,by,bx,by,ax,ay,ax,by,
               bx,ay,bx,by,ax,ay,bx,ay,ax,by,bx,by);
}

/** Dumps particle IDs and positions to a file.
 * \param[in] fp a file handle to write to. */
void container_2d::draw_particles(FILE *fp) {
    for(iterator cli=begin();cli<end();cli++) {
        int ij=cli->ijk,q=cli->q;
        double *pp=p[ij]+2*q;
        fprintf(fp,"%d %g %g\n",id[ij][q],*pp,pp[1]);
    }
}

/** Dumps particle positions in POV-Ray format.
 * \param[in] fp a file handle to write to. */
void container_2d::draw_particles_pov(FILE *fp) {
    for(iterator cli=begin();cli<end();cli++) {
        int ij=cli->ijk,q=cli->q;
        double *pp=p[ij]+2*q;
        fprintf(fp,"// id %d\nsphere{<%g,%g,0>,s}\n",id[ij][q],*pp,pp[1]);
    }
}

/** Computes Voronoi cells and saves the output in Gnuplot format.
 * \param[in] fp a file handle to write to. */
void container_2d::draw_cells_gnuplot(FILE *fp) {
    voronoicell_2d c;
    for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
        double *pp=p[cli->ijk]+2*cli->q;
        c.draw_gnuplot(*pp,pp[1],fp);
    }
}

/** Computes Voronoi cells and saves the output in POV-Ray format.
 * \param[in] fp a file handle to write to. */
void container_2d::draw_cells_pov(FILE *fp) {
    voronoicell_2d c;
    for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
        int ij=cli->ijk,q=cli->q;
        double *pp=p[ij]+2*q;
        fprintf(fp,"// cell %d\n",id[ij][q]);
        c.draw_pov(*pp,pp[1],fp);
    }
}

/** Computes the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_2d::print_custom(const char *format,FILE *fp) {
    int ij,q;double *pp;
    if(voro_contains_neighbor(format)) {
        voronoicell_neighbor_2d c;
        for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
            ij=cli->ijk;q=cli->q;
            pp=p[ij]+2*q;
            c.output_custom(format,id[ij][q],*pp,pp[1],default_radius,fp);
        }
    } else {
        voronoicell_2d c;
        for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
            ij=cli->ijk;q=cli->q;
            pp=p[ij]+2*q;
            c.output_custom(format,id[ij][q],*pp,pp[1],default_radius,fp);
        }
    }
}

/** Computes all of the Voronoi cells in the container, but does nothing with
 * the output. It is useful for measuring the pure computation time of the
 * Voronoi algorithm, without any additional calculations such as volume
 * evaluation or cell output. */
void container_2d::compute_all_cells() {
    voronoicell_2d c;
    for(iterator cli=begin();cli<end();cli++) compute_cell(c,cli);
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_2d::sum_cell_areas() {
    voronoicell_2d c;
    double area=0;
    for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) area+=c.area();
    return area;
}

/** Dumps particle IDs and positions to a file.
 * \param[in] fp a file handle to write to. */
void container_poly_2d::draw_particles(FILE *fp) {
    for(iterator cli=begin();cli<end();cli++) {
        int ij=cli->ijk,q=cli->q;
        double *pp=p[ij]+3*q;
        fprintf(fp,"%d %g %g %g\n",id[ij][q],*pp,pp[1],pp[2]);
    }
}

/** Dumps particle positions in POV-Ray format.
 * \param[in] fp a file handle to write to. */
void container_poly_2d::draw_particles_pov(FILE *fp) {
    for(iterator cli=begin();cli<end();cli++) {
        int ij=cli->ijk,q=cli->q;
        double *pp=p[ij]+3*q;
        fprintf(fp,"// id %d\nsphere{<%g,%g,0>,%g}\n",id[ij][q],*pp,pp[1],pp[2]);
    }
}

/** Computes Voronoi cells and saves the output in Gnuplot format.
 * \param[in] fp a file handle to write to. */
void container_poly_2d::draw_cells_gnuplot(FILE *fp) {
    voronoicell_2d c;
    for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
        double *pp=p[cli->ijk]+3*cli->q;
        c.draw_gnuplot(*pp,pp[1],fp);;
    }
}

/** Computes all Voronoi cells and saves the output in POV-Ray format.
 * \param[in] fp a file handle to write to. */
void container_poly_2d::draw_cells_pov(FILE *fp) {
    voronoicell_2d c;
    for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
        int ij=cli->ijk,q=cli->q;
        double *pp=p[ij]+3*q;
        fprintf(fp,"// cell %d\n",id[ij][q]);
        c.draw_pov(*pp,pp[1],fp);
    }
}

/** Computes the Voronoi cells and saves customized information about
 * them.
 * \param[in] cli the iterator class to use.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_poly_2d::print_custom(const char *format,FILE *fp) {
    int ij,q;double *pp;
    if(voro_contains_neighbor(format)) {
        voronoicell_neighbor_2d c;
        for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
            ij=cli->ijk;q=cli->q;
            pp=p[ij]+3*q;
            c.output_custom(format,id[ij][q],*pp,pp[1],pp[2],fp);
        }
    } else {
        voronoicell_2d c;
        for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) {
            ij=cli->ijk;q=cli->q;
            pp=p[ij]+3*q;
            c.output_custom(format,id[ij][q],*pp,pp[1],pp[2],fp);
        }
    }
}

/** Computes all of the Voronoi cells in the container, but does nothing with
 * the output. It is useful for measuring the pure computation time of the
 * Voronoi algorithm, without any additional calculations such as volume
 * evaluation or cell output. */
void container_poly_2d::compute_all_cells() {
    voronoicell_2d c;
    for(iterator cli=begin();cli<end();cli++) compute_cell(c,cli);
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_poly_2d::sum_cell_areas() {
    voronoicell_2d c;
    double area=0;
    for(iterator cli=begin();cli<end();cli++) if(compute_cell(c,cli)) area+=c.area();
    return area;
}

}
