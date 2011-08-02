// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file container.cc
 * \brief Function implementations for the container and related classes. */

#include "container_prd.hh"

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (bx_) exp. 
 * \param[in] (bxy_,by_) exp.
 * \param[in] (bxz_,byz_,bz_) exp.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_periodic_base::container_periodic_base(double bx_,double bxy_,double by_,
		double bxz_,double byz_,double bz_,int nx_,int ny_,int nz_,int init_mem,int ps_)
	: unitcell(bx_,bxy_,by_,bxz_,byz_,bz_), voropp_base(nx_,ny_,nz_,bx_/nx_,by_/ny_,bz_/nz_),
	ey(int(max_uv_y*ysp)), ez(int(max_uv_z*zsp)), wy(ny+ey), wz(nz+ez),
	oy(ny+2*ey), oz(nz+2*ez), oxyz(nx*oy*oz), id(new int*[oxyz]), p(new double*[oxyz]),
	co(new int[oxyz]), mem(new int[oxyz]), img(new char[oxyz]), ps(ps_) {
	int i,j,k,l;

	// Clear the global arrays
	int *pp(co);while(pp<co+oxyz) *(pp++)=0;
	pp=mem;while(pp<mem+oxyz) *(pp++)=0;
	char *cp(img);while(cp<img+oxyz) *(cp++)=0;

	// Set up memory for the blocks in the primary domain
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		mem[l]=init_mem;
		id[l]=new int[init_mem];
		p[l]=new double[ps*init_mem];
	}
}

/** The container destructor frees the dynamically allocated memory. */
container_periodic_base::~container_periodic_base() {
	int l;
	for(l=0;l<oxyz;l++) if(mem[l]>0) delete [] p[l];
	for(l=0;l<oxyz;l++) if(mem[l]>0) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] mem;
	delete [] co;
}

/** The class constructor sets up the geometry of container.
 * \param[in] (bx_) exp. 
 * \param[in] (bxy_,by_) exp.
 * \param[in] (bxz_,byz_,bz_) exp.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *			    coordinate directions.
 * \param[in] init_mem the initial memory allocation for each block. */
container_periodic::container_periodic(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
	int nx_,int ny_,int nz_,int init_mem)
	: container_periodic_base(bx_,bxy_,by_,bxz_,byz_,bz_,nx_,ny_,nz_,init_mem,3),
	vc(*this,2*nx_+1,2*ey+1,2*ez+1) {}

/** The class constructor sets up the geometry of container.
 * \param[in] (bx_) exp. 
 * \param[in] (bxy_,by_) exp.
 * \param[in] (bxz_,byz_,bz_) exp.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *			    coordinate directions.
 * \param[in] init_mem the initial memory allocation for each block. */
container_periodic_poly::container_periodic_poly(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
	int nx_,int ny_,int nz_,int init_mem)
	: container_periodic_base(bx_,bxy_,by_,bxz_,byz_,bz_,nx_,ny_,nz_,init_mem,4),
	max_radius(0), vc(*this,2*nx_+1,2*ey+1,2*ez+1) {}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_periodic::put(int n,double x,double y,double z) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	double *pp(p[ijk]+3*co[ijk]++);
	*(pp++)=x;*(pp++)=y;*pp=z;
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_periodic_poly::put(int n,double x,double y,double z,double r) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	double *pp(p[ijk]+4*co[ijk]++);
	*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
	if(max_radius<r) max_radius=r;
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
void container_periodic::put(voropp_order &vo,int n,double x,double y,double z) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	vo.add(ijk,co[ijk]);
	double *pp(p[ijk]+3*co[ijk]++);
	*(pp++)=x;*(pp++)=y;*pp=z;
}

/** Put a particle into the correct region of the container, also recording
 * into which region it was stored.
 * \param[in] vo the ordering class in which to record the region.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void container_periodic_poly::put(voropp_order &vo,int n,double x,double y,double z,double r) {
	int ijk;
	put_locate_block(ijk,x,y,z);
	id[ijk][co[ijk]]=n;
	vo.add(ijk,co[ijk]);
	double *pp(p[ijk]+4*co[ijk]++);
	*(pp++)=x;*(pp++)=y;*(pp++)=z;*pp=r;
	if(max_radius<r) max_radius=r;
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
void container_periodic_base::put_locate_block(int &ijk,double &x,double &y,double &z) {
	int k=step_int(z*zsp);
	if(k<0||k>=nz) {
		int ak=step_div(k,nz);
		z-=ak*bz;y-=ak*byz;x-=ak*bxz;k-=ak*nz;
	}
	int j=step_int(y*ysp);
	if(j<0||j>=ny) {
		int aj=step_div(j,ny);
		y-=aj*by;x-=aj*bxy;j-=aj*ny;
	}
	ijk=step_int(x*xsp);
	if(ijk<0||ijk>=nx) {
		int ai=step_div(ijk,nx);
		x-=ai*bx;ijk-=ai*nx;
	}
	j+=ey;k+=ez;
	ijk+=nx*(j+oy*k);
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_periodic_base::add_particle_memory(int i) {
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
 * \param[in] fp the file handle to read from. */
void container_periodic::import(FILE *fp) {
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
void container_periodic::import(voropp_order &vo,FILE *fp) {
	int i,j;
	double x,y,z;
	while((j=fscanf(fp,"%d %lg %lg %lg",&i,&x,&y,&z))==4) put(vo,i,x,y,z);
	if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Import a list of particles from an open file stream into the container.
 * Entries of five numbers (Particle ID, x position, y position, z position,
 * radius) are searched for. If the file cannot be successfully read, then the
 * routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void container_periodic_poly::import(FILE *fp) {
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
void container_periodic_poly::import(voropp_order &vo,FILE *fp) {
	int i,j;
	double x,y,z,r;
	while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(vo,i,x,y,z,r);
	if(j!=EOF) voropp_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Outputs the a list of all the container regions along with the number of
 * particles stored within each. */
void container_periodic_base::region_count() {
	int i,j,k,*cop(co);
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
		printf("Region (%d,%d,%d): %d particles\n",i,j,k,*(cop++));
}

/** Clears a container of particles. */
void container_periodic::clear() {
	for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
}

/** Clears a container of particles, also clearing resetting the maximum radius
 * to zero. */
void container_periodic_poly::clear() {
	for(int *cop=co;cop<co+nxyz;cop++) *cop=0;
	max_radius=0;
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_periodic::print_custom(const char *format,FILE *fp) {
	v_loop_all vl(*this);
	print_custom(vl,format,fp);
}

/** Computes all the Voronoi cells and saves customized
 * information about them.
 * \param[in] format the custom output string to use.
 * \param[in] fp a file handle to write to. */
void container_periodic_poly::print_custom(const char *format,FILE *fp) {
	v_loop_all vl(*this);
	print_custom(vl,format,fp);
}

/** Computes all the Voronoi cells and saves customized information about them.
 * \param[in] format the custom output string to use.
 * \param[in] filename the name of the file to write to. */
void container_periodic::print_custom(const char *format,const char *filename) {
	FILE *fp(voropp_safe_fopen(filename,"w"));
	print_custom(format,fp);
	fclose(fp);
}

/** Computes all the Voronoi cells and saves customized
 * information about them
 * \param[in] format the custom output string to use.
 * \param[in] filename the name of the file to write to. */
void container_periodic_poly::print_custom(const char *format,const char *filename) {
	FILE *fp(voropp_safe_fopen(filename,"w"));
	print_custom(format,fp);
	fclose(fp);
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
void container_periodic::compute_all_cells() {
	voronoicell c;
	v_loop_all vl(*this);
	if(vl.start()) do compute_cell(c,vl);
	while(vl.inc());
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
void container_periodic_poly::compute_all_cells() {
	voronoicell c;
	v_loop_all vl(*this);
	if(vl.start()) do compute_cell(c,vl);while(vl.inc());
}

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
double container_periodic::sum_cell_volumes() {
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
double container_periodic_poly::sum_cell_volumes() {
	voronoicell c;
	double vol=0;
	v_loop_all vl(*this);
	if(vl.start()) do if(compute_cell(c,vl)) vol+=c.volume();while(vl.inc());
	return vol;
}

/** This routine creates all periodic images of the particles. It is meant for
 * diagnostic purposes only, since usually periodic images are dynamically
 * created in when they are referenced. */
void container_periodic_base::create_all_images() {
	int i,j,k;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) create_periodic_image(i,j,k);
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] os an output stream to write to. */
void container_periodic_base::check_compartmentalized() {
	int c,l,i,j,k;
	double mix,miy,miz,max,may,maz;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		if(mem[l]==0) continue;
		mix=i*boxx-tolerance;max=mix+boxx+tolerance;
		miy=(j-ey)*boxy-tolerance;may=miy+boxy+tolerance;
		miz=(k-ez)*boxz-tolerance;maz=miz+boxz+tolerance;
		for(c=0;c<co[l];c++) if(p[l][ps*c]<mix||p[l][ps*c]>max
				      ||p[l][ps*c+1]<miy||p[l][ps*c+1]>may
				      ||p[l][ps*c+2]<miz||p[l][ps*c+2]>maz) printf("%d %d %d %f %f %f %f %f %f %f %f %f\n",i,j,k,p[l][ps*c],p[l][ps*c+1],p[l][ps*c+2],mix,max,miy,may,miz,maz);
	}
}

void container_periodic_base::create_side_image(int di,int dj,int dk) {
	double boxx=bx/nx;
	int l,dijk=di+nx*(dj+oy*dk),odijk,ima=step_div(dj-ey,ny);
	int qua=di+step_int(-ima*bxy*xsp),quadiv=step_div(qua,nx);
	int fi=qua-quadiv*nx,fijk=fi+nx*(dj-ima*ny+oy*dk);
	double dis=ima*bxy+quadiv*bx,switchx=di*boxx-ima*bxy-quadiv*bx,adis;

	// Left image computation
	if((img[dijk]&1)==0) {
		if(di>0) {
			odijk=dijk-1;adis=dis;
		} else {
			odijk=dijk+nx-1;adis=dis+bx;
		}
		img[odijk]|=2;
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l]>switchx) quick_put(dijk,fijk,l,dis,by*ima,0);
			else quick_put(odijk,fijk,l,adis,by*ima,0);
		}
	}

	// Right image computation
	if((img[dijk]&2)==0) {
		if(fi==nx-1) {
			fijk+=1-nx;switchx+=(1-nx)*boxx;dis+=bx;
		} else {
			fijk++;switchx+=boxx;
		}
		if(di==nx-1) {
			odijk=dijk-nx+1;adis=dis-bx;
		} else {
			odijk=dijk+1;adis=dis;
		}
		img[odijk]|=1;
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l]<switchx) quick_put(dijk,fijk,l,dis,by*ima,0);
			else quick_put(odijk,fijk,l,adis,by*ima,0);
		}
	}

	// All contributions to the block now added, so set both two bits of
	// the image information
	img[dijk]=3;
}

void container_periodic_base::create_vertical_image(int di,int dj,int dk) {
	int l,dijk=di+nx*(dj+oy*dk),dijkl,dijkr,ima=step_div(dk-ez,nz);
	int qj=dj+step_int(-ima*byz*ysp),qjdiv=step_div(qj-ey,ny);
	int qi=di+step_int((-ima*bxz-qjdiv*bxy)*xsp),qidiv=step_div(qi,nx);
	int fi=qi-qidiv*nx,fj=qj-qjdiv*ny,fijk=fi+nx*(fj+oy*(dk-ima*nz)),fijk2;
	double disy=ima*byz+qjdiv*by,switchy=(dj-ey)*boxy-ima*byz-qjdiv*by;
	double disx=ima*bxz+qjdiv*bxy+qidiv*bx,switchx=di*boxx-ima*bxz-qjdiv*bxy-qidiv*bx;
	double switchx2,disxl,disxr,disx2,disxr2;

	if(di==0) {dijkl=dijk+nx-1;disxl=disx+bx;}
	else {dijkl=dijk-1;disxl=disx;}

	if(di==nx-1) {dijkr=dijk-nx+1;disxr=disx-bx;}
	else {dijkr=dijk+1;disxr=disx;}

	// Down-left image computation
	bool y_exist=dj!=0;
	if((img[dijk]&1)==0) {
		img[dijkl]|=2;
		if(y_exist) {
			img[dijkl-nx]|=8;
			img[dijk-nx]|=4;
		}
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l+1]>switchy) {
				if(p[fijk][ps*l]>switchx) quick_put(dijk,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl,fijk,l,disxl,disy,bz*ima);
			} else {
				if(!y_exist) continue;
				if(p[fijk][ps*l]>switchx) quick_put(dijk-nx,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl-nx,fijk,l,disxl,disy,bz*ima);
			}
		}
	}

	// Down-right image computation
	if((img[dijk]&2)==0) {
		if(fi==nx-1) {
			fijk2=fijk+1-nx;switchx2=switchx+(1-nx)*boxx;disx2=disx+bx;disxr2=disxr+bx;
		} else {
			fijk2=fijk+1;switchx2=switchx+boxx;disx2=disx;disxr2=disxr;
		}
		img[dijkr]|=1;
		if(y_exist) {
			img[dijkr-nx]|=4;
			img[dijk-nx]|=8;
		}
		for(l=0;l<co[fijk2];l++) {
			if(p[fijk2][ps*l+1]>switchy) {
				if(p[fijk2][ps*l]>switchx2) quick_put(dijkr,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk,fijk2,l,disx2,disy,bz*ima);
			} else {
				if(!y_exist) continue;
				if(p[fijk2][ps*l]>switchx2) quick_put(dijkr-nx,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk-nx,fijk2,l,disx2,disy,bz*ima);
			}
		}
	}

	// Recomputation for boundary cases
	if(fj==wy-1) {
		fijk+=nx*(1-ny)-fi;
		switchy+=(1-ny)*boxy;
		disy+=by;
		qi=di+step_int(-(ima*bxz+(qjdiv+1)*bxy)*xsp);
		int dqidiv=step_div(qi,nx)-qidiv;qidiv+=dqidiv;
		fi=qi-qidiv*nx;
		fijk+=fi;
		disx+=bxy+bx*dqidiv;
		disxl+=bxy+bx*dqidiv;
		disxr+=bxy+bx*dqidiv;
		switchx-=bxy+bx*dqidiv;
	} else {
		fijk+=nx;switchy+=boxy;
	}

	// Up-left image computation
	y_exist=dj!=oy-1;
	if((img[dijk]&4)==0) {
		img[dijkl]|=8;
		if(y_exist) {
			img[dijkl+nx]|=2;
			img[dijk+nx]|=1;
		}
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][ps*l+1]>switchy) {
				if(!y_exist) continue;
				if(p[fijk][ps*l]>switchx) quick_put(dijk+nx,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl+nx,fijk,l,disxl,disy,bz*ima);
			} else {
				if(p[fijk][ps*l]>switchx) quick_put(dijk,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl,fijk,l,disxl,disy,bz*ima);
			}
		}
	}

	// Up-right image computation
	if((img[dijk]&8)==0) {
		if(fi==nx-1) {
			fijk2=fijk+1-nx;switchx2=switchx+(1-nx)*boxx;disx2=disx+bx;disxr2=disxr+bx;
		} else {
			fijk2=fijk+1;switchx2=switchx+boxx;disx2=disx;disxr2=disxr;
		}
		img[dijkr]|=4;
		if(y_exist) {
			img[dijkr+nx]|=1;
			img[dijk+nx]|=2;
		}
		for(l=0;l<co[fijk2];l++) {
			if(p[fijk2][ps*l+1]>switchy) {
				if(!y_exist) continue;
				if(p[fijk2][ps*l]>switchx2) quick_put(dijkr+nx,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk+nx,fijk2,l,disx2,disy,bz*ima);
			} else {
				if(p[fijk2][ps*l]>switchx2) quick_put(dijkr,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk,fijk2,l,disx2,disy,bz*ima);
			}
		}
	}
	
	// All contributions to the block now added, so set all four bits of
	// the image information
	img[dijk]=15;
}

inline void container_periodic_base::quick_put(int reg,int fijk,int l,double dx,double dy,double dz) {
	if(co[reg]==mem[reg]) add_particle_memory(reg);
	double *p1(p[reg]+ps*co[reg]),*p2(p[fijk]+ps*l);
	*(p1++)=*(p2++)+dx;
	*(p1++)=*(p2++)+dy;
	*p1=*p2+dz;
	if(ps==4) *(++p1)=*(++p2);
	id[reg][co[reg]++]=id[fijk][l];
}
