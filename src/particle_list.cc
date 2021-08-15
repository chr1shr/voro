// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file particle_list.cc
 * \brief Function implementations for the particle_list_base and related
 * classes. */

#include <cmath>

#include "config.hh"
#include "particle_list.hh"
#include "container.hh"
#include "container_tri.hh"

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction. It allocates an initial
 * chunk into which to store particle information.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (x_prd_,y_prd_,z_prd_) flags setting whether the container is
 *                                   periodic in each coordinate direction.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
particle_list_base::particle_list_base(int ps_) : ps(ps_),
    index_sz(init_chunk_size), pre_id(new int*[index_sz]),
    end_id(pre_id), pre_p(new double*[index_sz]), end_p(pre_p) {
    ch_id=*end_id=new int[particle_list_chunk_size];
    l_id=end_id+index_sz;e_id=ch_id+particle_list_chunk_size;
    ch_p=*end_p=new double[ps*particle_list_chunk_size];
}

/** The destructor frees the dynamically allocated memory. */
particle_list_base::~particle_list_base() {
    delete [] *end_p;
    delete [] *end_id;
    while (end_id!=pre_id) {
        end_p--;
        delete [] *end_p;
        end_id--;
        delete [] *end_id;
    }
    delete [] pre_p;
    delete [] pre_id;
}

/** Guesses the optimal grid of blocks to use for a two-dimensional
 * computation, by assuming that the particles are evenly distributed in space,
 * and aiming for the blocks to be approximately squares.
 * \param[in] (lx,ly) the dimensions of the container.
 * \param[out] (nx,ny) the number of blocks to use. */
void particle_list_base::guess_optimal(double lx,double ly,int &nx,int &ny) {
    double ilscale=sqrt(total_particles()/(optimal_particles_2d*lx*ly));
    nx=int(lx*ilscale+1);
    ny=int(ly*ilscale+1);
}

/** Guesses the optimal grid of blocks to use for a three-dimensional
 * computation, by assuming that the particles are evenly distributed in space,
 * and aiming for the blocks to be approximately cubes.
 * \param[in] (lx,ly,lz) the dimensions of the container.
 * \param[out] (nx,ny,nz) the number of blocks to use. */
void particle_list_base::guess_optimal(double lx,double ly,double lz,int &nx,int &ny,int &nz) {
    double ilscale=pow(total_particles()/(optimal_particles_3d*lx*ly*lz),1/3.0);
    nx=int(lx*ilscale+1);
    ny=int(ly*ilscale+1);
    nz=int(lz*ilscale+1);
}

/** Stores a particle ID and position, allocating a new memory chunk if
 * necessary.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y) the position vector of the inserted particle. */
void particle_list2::put(int n,double x,double y,double c) {
    if(ch_id==e_id) new_chunk();
    *(ch_id++)=n;
    *(ch_p++)=x;*(ch_p++)=y;*(ch_p++)=c;
}

/** Stores a particle ID and position, allocating a new memory chunk if
 * necessary.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,c) the position vector (and possibly radius) of the inserted
 *                    particle. */
void particle_list3::put(int n,double x,double y,double c) {
    if(ch_id==e_id) new_chunk();
    *(ch_id++)=n;
    *(ch_p++)=x;*(ch_p++)=y;*(ch_p++)=c;
}

/** Stores a particle ID and position, allocating a new memory chunk if necessary.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void particle_list4::put(int n,double x,double y,double z,double r) {
    if(ch_id==e_id) new_chunk();
    *(ch_id++)=n;
    *(ch_p++)=x;*(ch_p++)=y;*(ch_p++)=z;*(ch_p++)=r;
}

/** Transfers the particles stored within the class to a container class.
 * \param[in] con the container class to transfer to. */
void particle_list::setup(container &con) {
    int **c_id=pre_id,*idp,*ide,n;
    double **c_p=pre_p,*pp,x,y,z;
    while(c_id<end_id) {
        idp=*(c_id++);ide=idp+particle_list_chunk_size;
        pp=*(c_p++);
        while(idp<ide) {
            n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);
            con.put(n,x,y,z);
        }
    }
    idp=*c_id;
    pp=*c_p;
    while(idp<ch_id) {
        n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);
        con.put(n,x,y,z);
    }
}

/** Transfers the particles stored within the class to a container_poly class.
 * \param[in] con the container_poly class to transfer to. */
void particle_list_poly::setup(container_poly &con) {
    int **c_id=pre_id,*idp,*ide,n;
    double **c_p=pre_p,*pp,x,y,z,r;
    while(c_id<end_id) {
        idp=*(c_id++);ide=idp+particle_list_chunk_size;
        pp=*(c_p++);
        while(idp<ide) {
            n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
            con.put(n,x,y,z,r);
        }
    }
    idp=*c_id;
    pp=*c_p;
    while(idp<ch_id) {
        n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
        con.put(n,x,y,z,r);
    }
}

/** Transfers the particles stored within the class to a container class, also
 * recording the order in which particles were stored.
 * \param[in] vo the ordering class to use.
 * \param[in] con the container class to transfer to. */
void particle_list::setup(particle_order &vo,container &con) {
    int **c_id=pre_id,*idp,*ide,n;
    double **c_p=pre_p,*pp,x,y,z;
    while(c_id<end_id) {
        idp=*(c_id++);ide=idp+particle_list_chunk_size;
        pp=*(c_p++);
        while(idp<ide) {
            n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);
            con.put(vo,n,x,y,z);
        }
    }
    idp=*c_id;
    pp=*c_p;
    while(idp<ch_id) {
        n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);
        con.put(vo,n,x,y,z);
    }
}

/** Transfers the particles stored to a container class, also recording the
 * order in which particles were stored.
 * \param[in] vo the ordering class to use.
 * \param[in] con the container class to transfer to. */
template<class con_class>
void particle_list4::setup(con_class &con) {
    int **c_id=pre_id,*idp,*ide,n;
    double **c_p=pre_p,*pp,x,y,z,r;
    while(c_id<end_id) {
        idp=*(c_id++);ide=idp+particle_list_chunk_size;
        pp=*(c_p++);
        while(idp<ide) {
            n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
            con.put(vo,n,x,y,z,r);
        }
    }
    idp=*c_id;
    pp=*c_p;
    while(idp<ch_id) {
        n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
        con.put(vo,n,x,y,z,r);
    }
}

/** Transfers the particles stored to a container class, also recording the
 * order in which particles were stored.
 * \param[in] vo the ordering class to use.
 * \param[in] con the container class to transfer to. */
template<class con_class>
void particle_list4::setup(particle_order &vo,con_class &con) {
    int **c_id=pre_id,*idp,*ide,n;
    double **c_p=pre_p,*pp,x,y,z,r;
    while(c_id<end_id) {
        idp=*(c_id++);ide=idp+particle_list_chunk_size;
        pp=*(c_p++);
        while(idp<ide) {
            n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
            con.put(vo,n,x,y,z,r);
        }
    }
    idp=*c_id;
    pp=*c_p;
    while(idp<ch_id) {
        n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
        con.put(vo,n,x,y,z,r);
    }
}

/** Imports a list of particles from an open file stream. Entries of three
 * numbers (Particle ID, x position, y position) are searched for. If the file
 * cannot be successfully read, then the routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void particle_list2::import(FILE *fp) {
    int i,j;
    double x,y,z;
    while((j=fscanf(fp,"%d %lg %lg",&i,&x,&y))==3) put(i,x,y);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Imports a list of particles from an open file stream. Entries of four
 * numbers (Particle ID, x position, y position, extra) are searched for. The
 * extra number can either be the z position (for 3D Voronoi computations) or
 * the particle radius (for 2D radical Voronoi computations). If the file
 * cannot be successfully read, then the routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void particle_list3::import(FILE *fp) {
    int i,j;
    double x,y,z;
    while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&c))==4) put(i,x,y,c);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Imports a list of particles from an open file stream. Entries of five
 * numbers (Particle ID, x position, y position, z position, radius) are
 * searched for. If the file cannot be successfully read, then the routine
 * causes a fatal error.
 * \param[in] fp the file handle to read from. */
void particle_list4::import(FILE *fp) {
    int i,j;
    double x,y,z,r;
    while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(i,x,y,z,r);
    if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Allocates a new chunk of memory for storing particles. */
void particle_list_base::new_chunk() {
    end_id++;end_p++;
    if(end_id==l_id) extend_chunk_index();
    ch_id=*end_id=new int[particle_list_chunk_size];
    e_id=ch_id+particle_list_chunk_size;
    ch_p=*end_p=new double[ps*particle_list_chunk_size];
}

/** Extends the index of chunks. */
void particle_list_base::extend_chunk_index() {
    index_sz<<=1;
    if(index_sz>max_chunk_size)
        voro_fatal_error("Absolute memory limit on chunk index reached",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
    fprintf(stderr,"Particle list chunk index scaled up to %d\n",index_sz);
#endif
    int **n_id=new int*[index_sz],**p_id=n_id,**c_id=pre_id;
    double **n_p=new double*[index_sz],**p_p=n_p,**c_p=pre_p;
    while(c_id<end_id) {
        *(p_id++)=*(c_id++);
        *(p_p++)=*(c_p++);
    }
    delete [] pre_id;pre_id=n_id;end_id=p_id;l_id=pre_id+index_sz;
    delete [] pre_p;pre_p=n_p;end_p=p_p;
}

// Explicit instantiation
template void particle_list3::import(container_poly_2d&);
template void particle_list3::import(particle_order&,container_poly_2d&);
template void particle_list3::import(container_3d&);
template void particle_list3::import(particle_order&,container_3d&);
template void particle_list3::import(container_tri&);
template void particle_list3::import(particle_order&,container_tri&);
template void particle_list4::import(container_poly_3d&);
template void particle_list4::import(particle_order&,container_poly_3d&);
template void particle_list4::import(container_poly_tri&);
template void particle_list4::import(particle_order&,container_poly_tri&);

}
