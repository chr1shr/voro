// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file container.cc
 * \brief Function implementations for the container_periodic_base template and related
 * classes. */

#include "cell.hh"
#include "container.hh"

/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (xa,xb) the minimum and maximum x coordinates.
 * \param[in] (ya,yb) the minimum and maximum y coordinates.
 * \param[in] (za,zb) the minimum and maximum z coordinates.
 * \param[in] (xn,yn,zn) the number of grid blocks in each of the three
 *                       coordinate directions.
 * \param[in] memi the initial memory allocation for each block. */
template<class r_option>
container_periodic_base<r_option>::container_periodic_base(fpoint xb,fpoint xyb,fpoint yb,
		fpoint xzb,fpoint yzb,fpoint zb,int xn,int yn,int zn,int memi)
	: bx(xb),bxy(xyb),by(yb),bxz(xzb),byz(yzb),bz(zb),
	xsp(xn/xb),ysp(yn/yb),zsp(zn/zb),nx(xn),ny(yn),nz(zn),
	nxyz(xn*yn*zn),init_mem(memi),
	mv(0),wall_number(0),current_wall_size(init_wall_size),radius(this),
	sz(radius.mem_size),mrad(new fpoint[hgridsq*seq_length]),
	walls(new wall*[init_wall_size]) {
	int i,j,k,l;

	// Additional memory allocation for network
	pts=new double*[nxyz];
	idmem=new int*[nxyz];
	ptsc=new int[nxyz];
	ptsmem=new int[nxyz];
	for(l=0;l<nxyz;l++) {
		pts[l]=new double[4*memi];
		idmem[l]=new int[memi];
		ptsc[l]=0;ptsmem[l]=memi;
	}

	edc=0;edmem=memi*nxyz;
	ed=new int*[edmem];
	pered=new unsigned int*[edmem];
	raded=new double*[edmem];
	nu=new int[edmem];
	numem=new int[edmem];

	reg=new int[edmem];
	regp=new int[edmem]; 
	for(l=0;l<edmem;l++) {
		ed[l]=new int[4];
		nu[l]=0;numem[l]=4;
	}
	for(l=0;l<edmem;l++) raded[l]=new double[4];
	for(l=0;l<edmem;l++) pered[l]=new unsigned int[4];

	nett=new int[init_vertices];
	perio=new unsigned int[init_vertices];
	netmem=init_vertices;

	compute_unit_cell();
	if(unitcell.p==0) voropp_fatal_error("Periodic cell vanished",VOROPP_INTERNAL_ERROR);
	fpoint *upts=unitcell.pts;
	fpoint mx=upts[0],my=upts[1],mz=upts[2];
	for(l=3;l<3*unitcell.p;l+=3) {
		if(upts[l]>mx) mx=upts[l];
		if(upts[l+1]>my) my=upts[l+1];
		if(upts[l+2]>mz) mz=upts[l+2];
	}

	ey=1+int(my*ysp);
	ez=1+int(mz*zsp);

	wy=ny+ey;
	wz=nz+ez;

	oy=ny+2*ey;
	oz=nz+2*ez;
	int oxyz=nx*oy*oz;

	hx=1+2*nx;
	hy=1+2*ey;
	hz=1+2*ez;
	hxy=hx*hy;hxyz=hxy*hz;

	id=new int*[oxyz];
	p=new fpoint*[oxyz];
	co=new int[oxyz];
	mem=new int[oxyz];
	img=new char[oxyz];
	
	for(l=0;l<oxyz;l++) {
		co[l]=mem[l]=0;
		img[l]=0;
	}
	//printf("%d %d %d %d %d\n",ez,wz,ey,wy,nz);
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		mem[l]=memi;
		id[l]=new int[memi];
		p[l]=new fpoint[sz*memi];
	}

	s_size=3*(3+hxy+hz*(hx+hy));
	sl=new int[s_size];
	mask=new unsigned int[hxyz];
	for(l=0;l<hxyz;l++) mask[l]=0;

	// Precompute the radius table used in the cell construction
	initialize_radii();
}

template<class r_option>
void container_periodic_base<r_option>::compute_unit_cell() {
	unitcell.init(-10*bx,10*bx,-10*by,10*by,-10*bz,10*bz);
	int i,j,l=1;
	while(l<20) {
		//cout << l << endl;
		if(unit_cell_intersect(l)) {
			unit_cell_apply(l,0,0);
			for(i=1;i<l;i++) {
				unit_cell_apply(l,i,0);
				unit_cell_apply(-l,i,0);
			}
			for(i=-l;i<=l;i++) unit_cell_apply(i,l,0);
			for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
				unit_cell_apply(l,j,i);
				unit_cell_apply(-j,l,i);
				unit_cell_apply(-l,-j,i);
				unit_cell_apply(j,-l,i);
			}
			for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) unit_cell_apply(i,j,l);
		} else return;
		l++;
	}
	voropp_fatal_error("Periodic cell computation failed",VOROPP_INTERNAL_ERROR);
}

template<class r_option>
inline void container_periodic_base<r_option>::unit_cell_apply(int i,int j,int k) {
	double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
	unitcell.plane(x,y,z);
	unitcell.plane(-x,-y,-z);
}

template<class r_option>
bool container_periodic_base<r_option>::unit_cell_intersect(int l) {
	int i,j;
	if(unit_cell_test(l,0,0)) return true;
	for(i=1;i<l;i++) {
		if(unit_cell_test(l,i,0)) return true;
		if(unit_cell_test(-l,i,0)) return true;
	}
	for(i=-l;i<=l;i++) if(unit_cell_test(i,l,0)) return true;
	for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
		if(unit_cell_test(l,j,i)) return true;
		if(unit_cell_test(-j,l,i)) return true;
		if(unit_cell_test(-l,-j,i)) return true;
		if(unit_cell_test(j,-l,i)) return true;
	}
	for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) if(unit_cell_test(i,j,l)) return true;
	return false;
}

template<class r_option>
inline bool container_periodic_base<r_option>::unit_cell_test(int i,int j,int k) {
	double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
	double rsq=x*x+y*y+z*z;
	return unitcell.plane_intersects(x,y,z,rsq);
}

/** Increase network memory for a particular region. */
template<class r_option>
void container_periodic_base<r_option>::add_network_memory(int l) {
	ptsmem[l]<<=1;
	if(ptsmem[l]>max_container_vertex_memory)
		voropp_fatal_error("Container vertex maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
	double *npts(new double[4*ptsmem[l]]);
	int *nidmem(new int[ptsmem[l]]);
	for(int i=0;i<4*ptsc[l];i++) npts[i]=pts[l][i];
	for(int i=0;i<ptsc[l];i++) nidmem[i]=idmem[l][i];
	delete [] pts[l];
	delete [] idmem[l];
	pts[l]=npts;
	idmem[l]=nidmem;
}

/** Increase edge network memory. */
template<class r_option>
void container_periodic_base<r_option>::add_edge_network_memory() {
	int i;
	edmem<<=1;
	int **ned(new int*[edmem]);
	double **nraded(new double*[edmem]);
	unsigned int **npered(new unsigned int*[edmem]);
	int *nnu(new int[edmem]);
	int *nnumem(new int[edmem]);
	int *nreg(new int[edmem]);
	int *nregp(new int[edmem]);
	for(i=0;i<edc;i++) {
		ned[i]=ed[i];
		nraded[i]=raded[i];
		npered[i]=pered[i];
		nnu[i]=nu[i];
		nnumem[i]=numem[i];
		nreg[i]=reg[i];
		nregp[i]=regp[i];
	}
	while(i<edmem) {
		ned[i]=new int[4];
		nnu[i]=0;nnumem[i]=4;
		nraded[i]=new double[4];
		npered[i++]=new unsigned int[4];
	}
	delete [] ed;ed=ned;
	delete [] raded;raded=nraded;
	delete [] pered;pered=npered;
	delete [] nu;nu=nnu;
	delete [] numem;numem=nnumem;
	delete [] reg;reg=nreg;
	delete [] regp;regp=nregp;
}

/** Increase a particular vertex memory. */
template<class r_option>
void container_periodic_base<r_option>::add_particular_vertex_memory(int l) {
	numem[l]<<=1;
	if(numem[l]>1024)
		voropp_fatal_error("Particular vertex maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
	int *ned(new int[numem[l]]);
	double *nraded(new double[numem[l]]);
	unsigned int *npered(new unsigned int[numem[l]]);
	for(int i=0;i<nu[l];i++) {
		ned[i]=ed[l][i];
		nraded[i]=raded[l][i];
		npered[i]=pered[l][i];
	}
	delete [] ed[l];
	delete [] raded[l];
	delete [] pered[l];
	ed[l]=ned;
	raded[l]=nraded;
	pered[l]=npered;
}


/** The container destructor frees the dynamically allocated memory. */
template<class r_option>
container_periodic_base<r_option>::~container_periodic_base() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] p;
	delete [] id;
	delete [] walls;
	delete [] mrad;
	delete [] sl;
	delete [] mask;
	delete [] mem;
	delete [] co;
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] os an output stream to write to. */
template<class r_option>
void container_periodic_base<r_option>::draw_particles(ostream &os) {
	int c,l,i,j,k,ll;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		if(mem[l]==0) continue;
		for(c=0;c<co[l];c++) {
			os << id[l][c];
			for(ll=sz*c;ll<sz*(c+1);ll++) os << " " << p[l][ll];
			os << "\n";
		}
	}
}

template<class r_option>
void container_periodic_base<r_option>::draw_domain(const char* filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	os << "0 0 0\n";
	os << bx << " 0 0\n";
	os << bx+bxy << " " << by << " 0\n";
	os << bxy << " " << by << " 0\n";
	os << bxy+bxz << " " << by+byz << " " << bz << "\n";
	os << bx+bxy+bxz << " " << by+byz << " " << bz << "\n";
	os << bx+bxz << " " << byz << " " << bz << "\n";
	os << bxz << " " << byz << " " << bz << "\n";
	os << "0 0 0\n";
	os << bxy << " " << by << " 0\n\n";
	
	os << bxz << " " << byz << " " << bz << "\n";
	os << bxy+bxz << " " << by+byz << " " << bz << "\n\n";
	
	os << bx << " 0 0\n";
	os << bx+bxz << " " << byz << " " << bz << "\n\n";
	
	os << bx+bxy << " " << by << " 0\n";
	os << bx+bxy+bxz << " " << by+byz << " " << bz << "\n\n";
	
	os.close();
}

template<class r_option>
void container_periodic_base<r_option>::create_all_images() {
	int i,j,k;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) create_periodic_image(i,j,k);
}

/** Dumps all the particle positions and identifies to a file.
 * \param[in] os an output stream to write to. */
template<class r_option>
void container_periodic_base<r_option>::check_compartmentalized() {
	const fpoint boxx=bx/nx,boxy=by/ny,boxz=bz/nz;	
	int c,l,i,j,k;
	fpoint mix,miy,miz,max,may,maz;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		if(mem[l]==0) continue;
		mix=i*boxx-tolerance;max=mix+boxx+tolerance;
		miy=(j-ey)*boxy-tolerance;may=miy+boxy+tolerance;
		miz=(k-ez)*boxz-tolerance;maz=miz+boxz+tolerance;
		for(c=0;c<co[l];c++) if(p[l][sz*c]<mix||p[l][sz*c]>max
				      ||p[l][sz*c+1]<miy||p[l][sz*c+1]>may
				      ||p[l][sz*c+2]<miz||p[l][sz*c+2]>maz) printf("%d %d %d %f %f %f %f %f %f %f %f %f\n",i,j,k,p[l][sz*c],p[l][sz*c+1],p[l][sz*c+2],mix,max,miy,may,miz,maz);
	}
}

/** An overloaded version of the draw_particles() routine, that just prints
 * to standard output. */
template<class r_option>
void container_periodic_base<r_option>::draw_particles() {
	draw_particles(cout);
}

/** An overloaded version of the draw_particles() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
template<class r_option>
void container_periodic_base<r_option>::draw_particles(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles(os);
	os.close();
}

/** Dumps all the particle positions in POV-Ray format.
 * \param[in] os an output stream to write to. */
template<class r_option>
void container_periodic_base<r_option>::draw_particles_pov(ostream &os) {
	int c,i,j,k,l;
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		l=i+nx*(j+oy*k);
		for(c=0;c<co[l];c++) {
			os << "// id " << id[l][c] << "\n";
			os << "sphere{<" << p[l][sz*c] << "," << p[l][sz*c+1] << ","
			   << p[l][sz*c+2] << ">,";
			radius.rad(os,l,c);
		 	os << "}\n";
		}
	}
}

/** An overloaded version of the draw_particles_pov() routine, that just prints
 * to standard output. */
template<class r_option>
void container_periodic_base<r_option>::draw_particles_pov() {
	draw_particles_pov(cout);
}

/** An overloaded version of the draw_particles_pov() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
template<class r_option>
void container_periodic_base<r_option>::draw_particles_pov(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_particles_pov(os);
	os.close();
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle. */
template<class r_option>
void container_periodic_base<r_option>::put(int n,fpoint x,fpoint y,fpoint z) {
	put(n,x,y,z,0.5);
}

/** Put a particle into the correct region of the container.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle.*/
template<class r_option>
void container_periodic_base<r_option>::put(int n,fpoint x,fpoint y,fpoint z,fpoint r) {
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
	int i=step_int(x*xsp);
	if(i<0||i>=nx) {
		int ai=step_div(i,nx);
		x-=ai*bx;i-=ai*nx;
	}
	j+=ey;k+=ez;
	i+=nx*(j+oy*k);
	if(co[i]==mem[i]) add_particle_memory(i);
	p[i][sz*co[i]]=x;p[i][sz*co[i]+1]=y;p[i][sz*co[i]+2]=z;
	radius.store_radius(i,co[i],r);
	cout << p[i][sz*co[i]] << " " << p[i][sz*co[i]+1] << " " << p[i][sz*co[i]+2] << " " << p[i][sz*co[i]+3] << endl;
	id[i][co[i]++]=n;
}

/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
template<class r_option>
void container_periodic_base<r_option>::add_particle_memory(int i) {
	if(mem[i]==0) {
		mem[i]=init_mem;
		id[i]=new int[init_mem];
		p[i]=new fpoint[sz*init_mem];
		return;
	}
	int *idp;fpoint *pp;
	int l,nmem=2*mem[i];
#if VOROPP_VERBOSE >=3
	cerr << "Particle memory in region " << i << " scaled up to " << nmem << endl;
#endif
	if(nmem>max_particle_memory)
		voropp_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new fpoint[sz*nmem];
	for(l=0;l<sz*co[i];l++) pp[l]=p[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}

/** Add list memory. */
template<class r_option>
inline void container_periodic_base<r_option>::add_list_memory() {
	int i,j=0,*ps;
	ps=new int[s_size*2];
#if VOROPP_VERBOSE >=2
	cerr << "List memory scaled up to " << s_size*2 << endl;
#endif
	if(s_start<=s_end) {
		for(i=s_start;i<s_end;i++) ps[j++]=sl[i];
	} else {
		for(i=s_start;i<s_size;i++) ps[j++]=sl[i];
		for(i=0;i<s_end;i++) ps[j++]=sl[i];
	}
	s_size*=2;
	s_start=0;s_end=j;
	delete [] sl;sl=ps;
}

/** Import a list of particles from standard input.
 * \param[in] is a standard input stream to read from. */
template<class r_option>
void container_periodic_base<r_option>::import(istream &is) {
	radius.import(is);
}

/** An overloaded version of the import routine, that reads the standard input.
 */
template<class r_option>
inline void container_periodic_base<r_option>::import() {
	import(cin);
}

/** An overloaded version of the import routine, that reads in particles from
 * a particular file.
 * \param[in] filename the name of the file to read from. */
template<class r_option>
inline void container_periodic_base<r_option>::import(const char *filename) {
	ifstream is;
	is.open(filename,ifstream::in);
	if(is.fail()) voropp_fatal_error("Unable to open file for import",VOROPP_FILE_ERROR);
	import(is);
	is.close();
}

/** Outputs the number of particles within each region. */
template<class r_option>
void container_periodic_base<r_option>::region_count() {
	int i,j,k,ijk;
	for(k=0;k<oz;k++) for(j=0;j<oy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		printf("Region (%d,%d,%d): ijk=%d, p=%d, mem=%d, img=%d\n",i,j-ey,k-ez,ijk,co[ijk],mem[ijk],img[ijk]);
	}
}

/** Clears a container of particles. */
template<class r_option>
void container_periodic_base<r_option>::clear() {
	voropp_fatal_error("Not supported at present",VOROPP_INTERNAL_ERROR);
/*	for(int ijk=0;ijk<oxyz;ijk++) co[ijk]=0;
	radius.clear_max();*/
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in gnuplot format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box. */
template<class r_option>
void container_periodic_base<r_option>::draw_cells_gnuplot(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px;
	voropp_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;y=p[s][sz*q+1];z=p[s][sz*q+2];
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				if(compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) c.draw_gnuplot(os,x,y,z);
			}
		}
	} while((s=l1.inc(px))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_gnuplot() that computes the Voronoi
 * cells for the entire simulation region and saves the output in gnuplot
 * format.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_periodic_base<r_option>::draw_cells_gnuplot(const char *filename) {
	voronoicell c;
	int i,j,k,s=0,q;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		s=i+nx*(j+oy*k);
		for(q=0;q<co[s];q++) if(compute_cell(c,i,j,k,s,q)) {
			c.draw_gnuplot(os,p[s][sz*q],p[s][sz*q+1],p[s][sz*q+2]);
		}
	}
	os.close();
}

/** Computes the Voronoi cells for all particles within a rectangular box,
 * and saves the output in POV-Ray format.
 * \param[in] filename the name of the file to write to.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box. */
template<class r_option>
void container_periodic_base<r_option>::draw_cells_pov(const char *filename,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	fpoint x,y,z,px;
	voropp_loop l1(this);
	int q,s;
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px);
	do {
		for(q=0;q<co[s];q++) {
			os << "// cell " << id[s][q] << "\n";
			x=p[s][sz*q]+px;
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				if(compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) c.draw_pov(os,x,y,z);
			}
		}
	} while((s=l1.inc(px))!=-1);
	os.close();
}

/** An overloaded version of draw_cells_pov() that computes the Voronoi
 * cells for the entire simulation region and saves the output in POV-Ray
 * format.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_periodic_base<r_option>::draw_cells_pov(const char *filename) {
	draw_cells_pov(filename,0,bx,0,by,0,bz);
}

/** Computes all of the Voronoi cells in the container, but does nothing
 * with the output. It is useful for measuring the pure computation time
 * of the Voronoi algorithm, without any additional calculations such as
 * volume evaluation or cell output. */
template<class r_option>
void container_periodic_base<r_option>::compute_all_cells() {
	voronoicell c;
	int i,j,k,ijk,q;
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		for(q=0;q<co[ijk];q++) compute_cell(c,i,j,k,ijk,q);
	}
}

/** Computes the Voronoi volumes for all the particles, and stores the results
 * according to the particle ID numbers in a floating point array that the user
 * has supplied. No bounds checking on the array is performed, so it is up to
 * the user to ensure that the array is large enough to store the computed
 * numbers.
 * \param[in] bb a pointer to an array to store the volumes. The volume of the
 *               particle with ID number n will be stored at bb[n]. */
template<class r_option>
void container_periodic_base<r_option>::store_cell_volumes(fpoint *bb) {
	voronoicell c;
	int i,j,k,ijk,q;
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		for(q=0;q<co[ijk];q++) bb[id[ijk][q]]=compute_cell(c,i,j,k,ijk,q)?c.volume():0;
	}
}

/** Computes the local packing fraction at a point, by summing the volumes
 * of all particles within a test sphere, and dividing by the sum of their
 * Voronoi volumes that were previously computed using the store_cell_volumes()
 * function.
 * \param[in] bb a pointer to an array holding the Voronoi volumes of the
 *               particles.
 * \param[in] (cx,cy,cz) the center of the test sphere.
 * \param[in] r the radius of the test sphere. */
template<class r_option>
fpoint container_periodic_base<r_option>::packing_fraction(fpoint *bb,fpoint cx,fpoint cy,fpoint cz,fpoint r) {
	voropp_loop l1(this);
	fpoint px,x,y,z,rsq=r*r,pvol=0,vvol=0;
	int q,s;
	s=l1.init(cx,cy,cz,r,px);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px-cx;
			y=p[s][sz*q+1]-cy;
			z=p[s][sz*q+2]-cz;
			if(x*x+y*y+z*z<rsq) {
				pvol+=radius.volume(s,q);
				vvol+=bb[id[s][q]];
			}
		}
	} while((s=l1.inc(px))!=-1);
	return vvol>tolerance?pvol/vvol*4.1887902047863909846168578443726:0;
}

/** Computes the local packing fraction at a point, by summing the volumes of
 * all particles within test box, and dividing by the sum of their Voronoi
 * volumes that were previous computed using the store_cell_volumes() function.
 * \param[in] bb an array holding the Voronoi volumes of the particles.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box. */
template<class r_option>
fpoint container_periodic_base<r_option>::packing_fraction(fpoint *bb,fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
	voropp_loop l1(this);
	fpoint x,y,z,px,pvol=0,vvol=0;
	int q,s;
	s=l1.init(xmin,xmax,ymin,ymax,zmin,zmax,px);
	do {
		for(q=0;q<co[s];q++) {
			x=p[s][sz*q]+px;
			y=p[s][sz*q+1];
			z=p[s][sz*q+2];
			if(x>xmin&&x<xmax&&y>ymin&&y<ymax&&z>zmin&&z<zmax) {
				pvol+=radius.volume(s,q);
				vvol+=bb[id[s][q]];
			}
		}
	} while((s=l1.inc(px))!=-1);
	return vvol>tolerance?pvol/vvol*4.1887902047863909846168578443726:0;
}

template<class r_option>
void container_periodic_base<r_option>::clear_network() {
	int l;
	edc=0;
	for(l=0;l<nxyz;l++) ptsc[l]=0;
	for(l=0;l<edmem;l++) nu[l]=0;
}

template<class r_option>
void container_periodic_base<r_option>::draw_network(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_network(os);
	os.close();
}

template<class r_option>
void container_periodic_base<r_option>::draw_network() {
	draw_network(cout);
}

template<class r_option>
void container_periodic_base<r_option>::draw_network(ostream &os) {
	voronoicell c;
	int i,j,k,l,ijk=0,q,ai,aj,ak;
	double x,y,z;
	clear_network();
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		for(q=0;q<co[ijk];q++) {
			x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
			if(compute_cell(c,i,j,k,ijk,q,x,y,z)) add_to_network(c,x,y,z);
		}
	}
	for(l=0;l<edc;l++) {
		for(q=0;q<nu[l];q++) {
			unpack_periodicity(pered[l][q],ai,aj,ak);
			if(ed[l][q]<l&&ai==0&&aj==0&&ak==0) continue;
			i=reg[l];j=4*regp[l];
			os << pts[i][j] << " " << pts[i][j+1] << " " << pts[i][j+2] << "\n";
			i=reg[ed[l][q]];j=3*regp[ed[l][q]];
			os << pts[i][j]+bx*ai+bxy*aj+bxz*ak << " " << pts[i][j+1]+by*aj+byz*ak << " " << pts[i][j+2]+bz*ak << "\n\n\n";
		}
	}
}

template<class r_option>
void container_periodic_base<r_option>::print_network(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_network(os);
	os.close();
}

template<class r_option>
void container_periodic_base<r_option>::print_network() {
	print_network(cout);
}

template<class r_option>
void container_periodic_base<r_option>::print_network(ostream &os) {
	voronoicell c;
	int i,j,k,l,ijk=0,q,ai,aj,ak;
	double *ptsp;
	double x,y,z;
	clear_network();
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		for(q=0;q<co[ijk];q++) {
			x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
			if(compute_cell(c,i,j,k,ijk,q,x,y,z)) add_to_network(c,x,y,z);
		}
	}
	os << "Vertex table:\n";
	for(l=0;l<edc;l++) {
		ptsp=pts[reg[l]];j=4*regp[l];
		os << l << " " << ptsp[j] << " " << ptsp[j+1] << " " << ptsp[j+2] << " " << ptsp[j+3] << "\n";
	}
	
	os << "\nEdge table:\n";
	for(l=0;l<edc;l++) {
		ptsp=pts[reg[l]];j=4*regp[l];
		x=ptsp[j];y=ptsp[j+1];z=ptsp[j+2];

		for(q=0;q<nu[l];q++) {
			unpack_periodicity(pered[l][q],ai,aj,ak);

			// Uncomment the next line to remove reverse edges
			// if(ed[l][q]<l&&ai==0&&aj==0&&ak==0) continue;

			os << l << " -> " << ed[l][q] << " " << raded[l][q] << " " << ai << " " << aj << " " << ak;

			// This section computes the length of the edge
			{
				ptsp=pts[reg[ed[l][q]]];j=4*regp[ed[l][q]];
				double x2=ptsp[j]+ai*bx+aj*bxy+ak*bxz-x,
				       y2=ptsp[j+1]+aj*by+ak*byz-y,
				       z2=ptsp[j+2]+ak*bz-z;
				os << " " << (sqrt(x2*x2+y2*y2+z2*z2));
			}
			os << "\n"; 
		}
	}
}

template<class r_option>
unsigned int container_periodic_base<r_option>::pack_periodicity(int i,int j,int k) {
	unsigned int pa=((unsigned int) (127+i));
	pa<<=8;pa+=((unsigned int) (127+j));
	pa<<=8;pa+=((unsigned int) (127+k));
	return pa;
}

template<class r_option>
inline void container_periodic_base<r_option>::unpack_periodicity(unsigned int pa,int &i,int &j,int &k) {
	i=((signed int) (pa>>16))-127;
	j=((signed int) ((pa>>8)&255))-127;
	k=((signed int) (pa&255))-127;
}

template<class r_option>
template<class n_option>
void container_periodic_base<r_option>::add_to_network(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
	int i,j,k,ijk,l,q,ai,aj,ak,bi,bj,bk;unsigned int cper;
	fpoint vx,vy,vz,wx,wy,wz,dx,dy,dz,dis;fpoint *pp;
	if(c.p>netmem) {
		do {
			netmem<<=1;
		} while(netmem<c.p);
		delete [] nett;
		delete [] perio;
		nett=new int[netmem];
		perio=new unsigned int[netmem];
	}
	for(l=0;l<c.p;l++) {
		vx=x+c.pts[3*l]*0.5;vy=y+c.pts[3*l+1]*0.5;vz=z+c.pts[3*l+2]*0.5;
		if(safe_search_previous(vx,vy,vz,ijk,q,cper)) {
			nett[l]=idmem[ijk][q];
			perio[l]=cper;
		} else {
			k=step_int(vz*zsp);
			if(k<0||k>=nz) {
				ak=step_div(k,nz);
				vz-=ak*bz;vy-=ak*byz;vx-=ak*bxz;k-=ak*nz;
			} else ak=0;
			j=step_int(vy*ysp);
			if(j<0||j>=ny) {
				aj=step_div(j,ny);
				vy-=aj*by;vx-=aj*bxy;j-=aj*ny;
			} else aj=0;
			i=step_int(vx*xsp);
			if(i<0||i>=nx) {
				ai=step_div(i,nx);
				vx-=ai*bx;i-=ai*nx;
			} else ai=0;
			perio[l]=pack_periodicity(ai,aj,ak);
			ijk=i+nx*(j+ny*k);
			if(edc==edmem) add_edge_network_memory();
			if(ptsc[ijk]==ptsmem[ijk]) add_network_memory(ijk);
			reg[edc]=ijk;regp[edc]=ptsc[ijk];
			pts[ijk][4*ptsc[ijk]]=vx;
			pts[ijk][4*ptsc[ijk]+1]=vy;
			pts[ijk][4*ptsc[ijk]+2]=vz;
			pts[ijk][4*ptsc[ijk]+3]=0.5*sqrt(c.pts[3*l]*c.pts[3*l]+c.pts[3*l+1]*c.pts[3*l+1]+c.pts[3*l+2]*c.pts[3*l+2]);
			idmem[ijk][ptsc[ijk]++]=edc;
			nett[l]=edc++;
		}
	}
	for(l=0;l<c.p;l++) {
		k=nett[l];
		unpack_periodicity(perio[l],ai,aj,ak);
		pp=pts[reg[k]]+(4*regp[k]);
		vx=pp[0]+ai*bx+aj*bxy+ak*bxz;vy=pp[1]+aj*by+ak*byz;vz=pp[2]+ak*bz;
		for(q=0;q<c.nu[l];q++) {
			j=nett[c.ed[l][q]];
			unpack_periodicity(perio[c.ed[l][q]],bi,bj,bk);
			cper=pack_periodicity(bi-ai,bj-aj,bk-ak);
			if(not_already_there(k,j,cper)) {
				pp=pts[reg[j]]+(4*regp[j]);
				wx=pp[0]+bi*bx+bj*bxy+bk*bxz;wy=pp[1]+bj*by+bk*byz;wz=pp[2]+bk*bz;
				dx=wx-vx;dy=wy-vy;dz=wz-vz;
				dis=(x-vx)*dx+(y-vy)*dy+(z-vz)*dz;
				dis/=dx*dx+dy*dy+dz*dz;
				if(dis<0) dis=0;if(dis>1) dis=1;
				wx=vx-x+dis*dx;wy=vy-y+dis*dy;wz=vz-z+dis*dz;
				if(nu[k]==numem[k]) add_particular_vertex_memory(k);
				ed[k][nu[k]]=j;
				raded[k][nu[k]]=sqrt(wx*wx+wy*wy+wz*wz);
				pered[k][nu[k]++]=cper;
			}
		}
	}
}

template<class r_option> 
bool container_periodic_base<r_option>::not_already_there(int k,int j,unsigned int cper) {
	for(int i=0;i<nu[k];i++) if(ed[k][i]==j&&pered[k][i]==cper) return false;
	return true;
}

template<class r_option>
bool container_periodic_base<r_option>::safe_search_previous(fpoint x,fpoint y,fpoint z,int &ijk,int &q,unsigned int &cper) {
	const fpoint tol=0.5*tolerance;
	if(search_previous(x+tol,y+tol,z+tol,ijk,q,cper)) return true;
	if(search_previous(x-tol,y+tol,z+tol,ijk,q,cper)) return true;
	if(search_previous(x+tol,y-tol,z+tol,ijk,q,cper)) return true;
	if(search_previous(x-tol,y-tol,z+tol,ijk,q,cper)) return true;
	if(search_previous(x+tol,y+tol,z-tol,ijk,q,cper)) return true;
	if(search_previous(x-tol,y+tol,z-tol,ijk,q,cper)) return true;
	if(search_previous(x+tol,y-tol,z-tol,ijk,q,cper)) return true;
	return search_previous(x-tolerance,y-tolerance,z-tolerance,ijk,q,cper);
}


template<class r_option>
bool container_periodic_base<r_option>::search_previous(fpoint x,fpoint y,fpoint z,int &ijk,int &q,unsigned int &cper) {
		
	int k=step_int(z*zsp);
	int ai,aj,ak;

	if(k<0||k>=nz) {
		ak=step_div(k,nz);
		z-=ak*bz;y-=ak*byz;x-=ak*bxz;k-=ak*nz;
	} else ak=0;

	int j=step_int(y*ysp);
	if(j<0||j>=ny) {
		aj=step_div(j,ny);
		y-=aj*by;x-=aj*bxy;j-=aj*ny;
	} else aj=0;

	int i=step_int(x*xsp);
	if(i<0||i>=nx) {
		ai=step_div(i,nx);
		x-=ai*bx;i-=ai*nx;
	} else ai=0;

	ijk=i+nx*(j+ny*k);
	
	for(q=0;q<ptsc[ijk];q++) if(abs(pts[ijk][4*q]-x)<tolerance&&abs(pts[ijk][4*q+1]-y)<tolerance&&abs(pts[ijk][4*q+2]-z)<tolerance) {
		cper=pack_periodicity(ai,aj,ak);
		return true;
	}
	return false;
}


/** Calculates all of the Voronoi cells and sums their volumes. In most cases
 * without walls, the sum of the Voronoi cell volumes should equal the volume
 * of the container to numerical precision.
 * \return The sum of all of the computed Voronoi volumes. */
template<class r_option>
fpoint container_periodic_base<r_option>::sum_cell_volumes() {
	voronoicell c;
	int i,j,k,ijk,q;fpoint vol=0;
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		for(q=0;q<co[ijk];q++) if (compute_cell(c,i,j,k,ijk,q)) vol+=c.volume();
	}
	return vol;
}

/** Computes the Voronoi cells for all particles in the container, and for each
 * cell, outputs a line containing custom information about the cell structure.
 * The output format is specified using an input string with control sequences
 * similar to the standard C printf() routine. Full information about the
 * control sequences is available at http://math.lbl.gov/voro++/doc/custom.html
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] os an open output stream to write to. */
template<class r_option>
void container_periodic_base<r_option>::print_all_custom(const char *format,ostream &os) {
	int fp=0;

	// Check to see if the sequence "%n" appears in the format sequence
	while(format[fp]!=0) {
		if(format[fp]=='%') {
			fp++;
			if(format[fp]=='n') {

				// If a "%n" is detected, then we're going to
				// need neighbor information during the custom
				// output, so use the voronoicell_neighbor
				// class
				voronoicell_neighbor c;
				print_all_custom_internal(c,format,os);
				return;
			} else if(format[fp]==0) break;
		}
		fp++;
	}

	// No "%n" was detected, so we can just use the regular voronoicell
	// class without computing neighbor information.
	voronoicell c;
	print_all_custom_internal(c,format,os);
}

/** An overloaded version of print_all_custom() that prints to standard output.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics. */
template<class r_option>
void container_periodic_base<r_option>::print_all_custom(const char *format) {
	print_all_custom(format,cout);
}

/** An overloaded version of print_all_custom(), which outputs the result to a
 * particular file.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
void container_periodic_base<r_option>::print_all_custom(const char *format,const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_custom(format,os);
	os.close();
}

/** The internal part of the print_all_custom() routine, that can be called
 * with either a voronoicell class (if no neighbor computations are needed) or
 * with a voronoicell_neighbor class (if neighbor computations are needed).
 * \param[in,out] c a Voronoi cell object to use for the computation.
 * \param[in] format the format of the output lines, using control sequences to
 *                   denote the different cell statistics.
 * \param[in] os an open output stream to write to. */
template<class r_option>
template<class n_option>
void container_periodic_base<r_option>::print_all_custom_internal(voronoicell_base<n_option> &c,const char *format,ostream &os) {
	fpoint x,y,z;
	int i,j,k,ijk=0,q,fp;
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) for(ijk=i+nx*(j+oy*k),q=0;q<co[ijk];q++) {
		x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
		if(!compute_cell(c,i,j,k,ijk,q,x,y,z)) continue;
		fp=0;
		while(format[fp]!=0) {
			if(format[fp]=='%') {
				fp++;
				switch(format[fp]) {

					// Particle-related output
					case 'i': os << id[ijk][q];break;
					case 'x': os << x;break;
					case 'y': os << y;break;
					case 'z': os << z;break;
					case 'q': os << x << " " << y << " " << z;break;
					case 'r': radius.print(os,ijk,q,false);break;

					// Vertex-related output
					case 'w': os << c.p;break;
					case 'p': c.output_vertices(os);break;
					case 'P': c.output_vertices(os,x,y,z);break;
					case 'o': c.output_vertex_orders(os);
					case 'm': os << 0.25*c.max_radius_squared();break;

					// Edge-related output
					case 'g': os << c.number_of_edges();break;
					case 'E': os << c.total_edge_distance();break;
					case 'e': c.output_face_perimeters(os);break;

					// Face-related output
					case 's': os << c.number_of_faces();break;
					case 'F': os << c.surface_area();break;
					case 'A': c.output_face_freq_table(os);break;
					case 'a': c.output_face_orders(os);break;
					case 'f': c.output_face_areas(os);break;
					case 't': c.output_face_vertices(os);break;
					case 'l': c.output_normals(os);break;
					case 'n': c.output_neighbors(os);break;

					// Volume-related output
					case 'v': os << c.volume();break;
					case 'c': {
							  fpoint cx,cy,cz;
							  c.centroid(cx,cy,cz);
							  os << cx << " " << cy << " " << cz;
						  } break;
					case 'C': {
							  fpoint cx,cy,cz;
							  c.centroid(cx,cy,cz);
							  os << x+cx << " " << y+cy << " " << z+cz;
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

/** The internal part of the print_all() and print_all_neighbor() routines, that
 * computes all of the Voronoi cells, and then outputs simple information about
 * them. The routine can be called with either a voronoicell class (if no
 * neighbor computations are needed) or with a voronoicell_neighbor class (if
 * neighbor computations are needed).
 * \param[in,out] c a Voronoi cell object to use for the computation.
 * \param[in] os an open output stream to write to. */
template<class r_option>
template<class n_option>
inline void container_periodic_base<r_option>::print_all_internal(voronoicell_base<n_option> &c,ostream &os) {
	fpoint x,y,z;
	int i,j,k,ijk,q;
	for(k=ez;k<wz;k++) for(j=ey;j<wy;j++) for(i=0;i<nx;i++) {
		ijk=i+nx*(j+oy*k);
		for(q=0;q<co[ijk];q++) {
			x=p[ijk][sz*q];y=p[ijk][sz*q+1];z=p[ijk][sz*q+2];
			os << id[ijk][q] << " " << x << " " << y << " " << z;
			radius.print(os,ijk,q);
			if(compute_cell(c,i,j,k,ijk,q,x,y,z)) {
				os << " " << c.volume();
				c.output_neighbors(os,true);
				os << "\n";
			} else os << " 0\n";
		}
	}
}

/** Prints a list of all particle labels, positions, and Voronoi volumes to the
 * standard output.
 * \param[in] os the output stream to print to. */
template<class r_option>
void container_periodic_base<r_option>::print_all(ostream &os) {
	voronoicell c;
	print_all_internal(c,os);
}

/** An overloaded version of print_all(), which just prints to standard output. */
template<class r_option>
void container_periodic_base<r_option>::print_all() {
	voronoicell c;
	print_all_internal(c,cout);
}

/** An overloaded version of print_all(), which outputs the result to a particular
 * file.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
inline void container_periodic_base<r_option>::print_all(const char* filename) {
	voronoicell c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_internal(c,os);
	os.close();
}

/** Prints a list of all particle labels, positions, Voronoi volumes, and a list
 * of neighboring particles to an output stream.
 * \param[in] os the output stream to print to.*/
template<class r_option>
void container_periodic_base<r_option>::print_all_neighbor(ostream &os) {
	voronoicell_neighbor c;
	print_all_internal(c,os);
}

/** An overloaded version of print_all_neighbor(), which just prints to
 * standard output. */
template<class r_option>
void container_periodic_base<r_option>::print_all_neighbor() {
	voronoicell_neighbor c;
	print_all_internal(c,cout);
}

/** An overloaded version of print_all_neighbor(), which outputs the results to
 * a particular file.
 * \param[in] filename the name of the file to write to. */
template<class r_option>
inline void container_periodic_base<r_option>::print_all_neighbor(const char* filename) {
	voronoicell_neighbor c;
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	print_all_internal(c,os);
	os.close();
}

/** Initialize the Voronoi cell to be the entire container. For non-periodic
 * coordinates, this is set by the position of the walls. For periodic
 * coordinates, the space is equally divided in either direction from the
 * particle's initial position. That makes sense since those boundaries would
 * be made by the neighboring periodic images of this particle. It also applies
 * plane cuts made by any walls that have been added to the container.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (x,y,z) the position of the particle.
 * \return False if the plane cuts applied by walls completely removed the
 *         cell, true otherwise. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::initialize_voronoicell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
/*	fpoint x1,x2,y1,y2,z1,z2;
	x1=-(x2=0.5*bx);
	y1=-(y2=0.5*by);
	z1=-(z2=0.5*bz);
	c.init(x1,x2,y1,y2,z1,z2);*/
	c.init(unitcell);
	for(int j=0;j<wall_number;j++) if(!(walls[j]->cut_cell(c,x,y,z))) return false;
	return true;
}

/** This function tests to see if a given vector lies within the container
 * bounds and any walls.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
template<class r_option>
bool container_periodic_base<r_option>::point_inside(fpoint x,fpoint y,fpoint z) {
	if(x<0||x>bx||y<0||y>by||z<0||z>bz) return false;
	return point_inside_walls(x,y,z);
}

/** This function tests to see if a give vector lies within the walls that have
 * been added to the container, but does not specifically check whether the
 * vector lies within the container bounds.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
template<class r_option>
bool container_periodic_base<r_option>::point_inside_walls(fpoint x,fpoint y,fpoint z) {
	for(int j=0;j<wall_number;j++) if(!walls[j]->point_inside(x,y,z)) return false;
	return true;
}

/** This routine is a simpler alternative to compute_cell(), that constructs
 * the cell by testing over successively larger spherical shells of particles.
 * For a container that is homogeneously filled with particles, this routine
 * runs as fast as compute_cell(). However, it rapidly becomes inefficient
 * for cases when the particles are not homogeneously distributed, or where
 * parts of the container might not be filled. In that case, the spheres may
 * grow very large before being cut off, leading to poor efficiency.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \param[in] (x,y,z) the coordinates of the particle.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
bool container_periodic_base<r_option>::compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z) {

	// This length scale determines how large the spherical shells should
	// be, and it should be set to approximately the particle diameter
	const fpoint length_scale=1;
	fpoint x1,y1,z1,qx,lr=0,lrs=0,ur,urs,rs;
	int q,t;
	voropp_loop l(this);
	if(!initialize_voronoicell(c,x,y,z)) return false;

	// Now the cell is cut by testing neighboring particles in concentric
	// shells. Once the test shell becomes twice as large as the Voronoi
	// cell we can stop testing.
	radius.init(ijk,s);
	while(radius.cutoff(lrs)<c.max_radius_squared()) {
		ur=lr+0.5*length_scale;urs=ur*ur;
		t=l.init(x,y,z,ur,qx);
		do {
			for(q=0;q<co[t];q++) {
				x1=p[t][sz*q]+qx-x;y1=p[t][sz*q+1]-y;z1=p[t][sz*q+2]-z;
				rs=x1*x1+y1*y1+z1*z1;
				if(lrs-tolerance<rs&&rs<urs&&(q!=s||ijk!=t)) {
					if(!c.nplane(x1,y1,z1,radius.scale(rs,t,q),id[t][q])) return false;
				}
			}
		} while((t=l.inc(qx))!=-1);
		lr=ur;lrs=urs;
	}
	return true;
}

/** A overloaded version of compute_cell_sphere(), that sets up the x, y, and z
 * variables.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::compute_cell_sphere(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s) {
	fpoint x=p[ijk][sz*s],y=p[ijk][sz*s+1],z=p[ijk][sz*s+2];
	return compute_cell_sphere(c,i,j,k,ijk,s,x,y,z);
}

/** A overloaded version of compute_cell, that sets up the x, y, and z variables.
 * It can be run by the user, and it is also called multiple times by the
 * functions print_all(), store_cell_volumes(), and the output routines.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s) {
	fpoint x=p[ijk][sz*s],y=p[ijk][sz*s+1],z=p[ijk][sz*s+2];
	return compute_cell(c,i,j,k,ijk,s,x,y,z);
}

/** This routine computes a Voronoi cell for a single particle in the
 * container. It can be called by the user, but is also forms the core part of
 * several of the main functions, such as store_cell_volumes(), print_all(),
 * and the drawing routines. The algorithm constructs the cell by testing over
 * the neighbors of the particle, working outwards until it reaches those
 * particles which could not possibly intersect the cell. For maximum
 * efficiency, this algorithm is divided into three parts. In the first
 * section, the algorithm tests over the blocks which are in the immediate
 * vicinity of the particle, by making use of one of the precomputed worklists.
 * The code then continues to test blocks on the worklist, but also begins to
 * construct a list of neighboring blocks outside the worklist which may need
 * to be test. In the third section, the routine starts testing these
 * neighboring blocks, evaluating whether or not a particle in them could
 * possibly intersect the cell. For blocks that intersect the cell, it tests
 * the particles in that block, and then adds the block neighbors to the list
 * of potential places to consider.
 * \param[in,out] c a reference to a voronoicell object.
 * \param[in] (i,j,k) the coordinates of the block that the test particle is
 *                    in.
 * \param[in] ijk the index of the block that the test particle is in, set to
 *                i+nx*(j+ny*k).
 * \param[in] s the index of the particle within the test block.
 * \param[in] (x,y,z) the coordinates of the particle.
 * \return False if the Voronoi cell was completely removed during the
 *         computation and has zero volume, true otherwise. */
template<class r_option>
template<class n_option>
bool container_periodic_base<r_option>::compute_cell(voronoicell_base<n_option> &c,int i,int j,int k,int ijk,int s,fpoint x,fpoint y,fpoint z) {
	const fpoint boxx=bx/nx,boxy=by/ny,boxz=bz/nz;
	fpoint x1,y1,z1,qx=0;
	fpoint xlo,ylo,zlo,xhi,yhi,zhi,rs;
	int ci,cj,ck,di,dj,dk,dijk,ei,ej,ek,eijk,si,sj,sk,sijk,iv;
	fpoint gxs,gys,gzs,*radp;
	int f,g,l;unsigned int q,*e;
	const unsigned int b1=1<<21,b2=1<<22,b3=1<<24,b4=1<<25,b5=1<<27,b6=1<<28;

	radius.init(ijk,s);

	// Initialize the Voronoi cell to fill the entire container
	if(!initialize_voronoicell(c,x,y,z)) return false;
	fpoint crs,mrs;

	int next_count=3,list_index=0,list_size=8;
	int count_list[]={7,11,15,19,26,35,45,59};

	// Test all particles in the particle's local region first
	for(l=0;l<s;l++) {
		x1=p[ijk][sz*l]-x;
		y1=p[ijk][sz*l+1]-y;
		z1=p[ijk][sz*l+2]-z;
		rs=radius.scale(x1*x1+y1*y1+z1*z1,ijk,l);
		if(!c.nplane(x1,y1,z1,rs,id[ijk][l])) return false;
	}
	l++;
	while(l<co[ijk]) {
		x1=p[ijk][sz*l]-x;
		y1=p[ijk][sz*l+1]-y;
		z1=p[ijk][sz*l+2]-z;
		rs=radius.scale(x1*x1+y1*y1+z1*z1,ijk,l);
		if(!c.nplane(x1,y1,z1,rs,id[ijk][l])) return false;
		l++;
	}

	// Now compute the maximum distance squared from the cell center to a
	// vertex. This is used to cut off the calculation since we only need
	// to test out to twice this range.
	mrs=c.max_radius_squared();

	// Now compute the fractional position of the particle within its
	// region and store it in (fx,fy,fz). We use this to compute an index
	// (si,sj,sk) of which subregion the particle is within.
	unsigned int m1,m2;
	fpoint fx=x-boxx*i,fy=y-boxy*(j-ey),fz=z-boxz*(k-ez);
	si=int(fx*xsp*fgrid);sj=int(fy*ysp*fgrid);sk=int(fz*zsp*fgrid);

	// The indices (si,sj,sk) tell us which worklist to use, to test the
	// blocks in the optimal order. But we only store worklists for the
	// eighth of the region where si, sj, and sk are all less than half the
	// full grid. The rest of the cases are handled by symmetry. In this
	// section, we detect for these cases, by reflecting high values of si,
	// sj, and sk. For these cases, a mask is constructed in m1 and m2
	// which is used to flip the worklist information when it is loaded.
	if(si>=hgrid) {
		gxs=fx;
		m1=127+(3<<21);si=fgrid-1-si;m2=1+(1<<21);
	} else {m1=m2=0;gxs=boxx-fx;}
	if(sj>=hgrid) {
		gys=fy;
		m1|=(127<<7)+(3<<24);sj=fgrid-1-sj;m2|=(1<<7)+(1<<24);
	} else gys=boxy-fy;
	if(sk>=hgrid) {
		gzs=fz;
		m1|=(127<<14)+(3<<27);sk=fgrid-1-sk;m2|=(1<<14)+(1<<27);
	} else gzs=boxz-fz;
	gxs*=gxs;gys*=gys;gzs*=gzs;

	// It's possible that a problem with the int() function could lead to
	// spurious values with particles lying on the boundaries of the regions.
	// In this section we correct for that.
	if(si<0) si=0;if(sj<0) sj=0;if(sk<0) sk=0;

	// Now compute which worklist we are going to use, and set radp and e to
	// point at the right offsets
	sijk=si+hgrid*(sj+hgrid*sk);
	radp=mrad+sijk*seq_length;
	e=(const_cast<unsigned int*> (wl))+sijk*seq_length;

	// Read in how many items in the worklist can be tested without having to
	// worry about writing to the mask
	f=e[0];g=0;
	do {

		// At the intervals specified by count_list, we recompute the
		// maximum radius squared
		if(g==next_count) {
			mrs=c.max_radius_squared();
			if(list_index!=list_size) next_count=count_list[list_index++];
		}

		// If mrs is less than the minimum distance to any untested
		// block, then we are done
		if(mrs<radius.cutoff(radp[g])) return true;
		g++;

		// Load in a block off the worklist, permute it with the
		// symmetry mask, and decode its position. These are all
		// integer bit operations so they should run very fast.
		q=e[g];q^=m1;q+=m2;
		di=q&127;di-=64;
		dj=(q>>7)&127;dj-=64;
		dk=(q>>14)&127;dk-=64;

		// Check that the worklist position is in range
		if(di<-nx) continue;else if(di>nx) continue;
		if(dj<-ey) continue;else if(dj>ey) continue;
		if(dk<-ez) continue;else if(dk>ez) continue;

		// Call the compute_min_max_radius() function. This returns
		// true if the minimum distance to the block is bigger than the
		// current mrs, in which case we skip this block and move on.
		// Otherwise, it computes the maximum distance to the block and
		// returns it in crs.
		if(compute_min_max_radius(di,dj,dk,fx,fy,fz,gxs,gys,gzs,crs,mrs)) continue;

		// Now compute which region we are going to loop over, adding a
		// displacement for the periodic cases
		di+=i;dj+=j;dk+=k;
		iv=step_div(di,nx);if(iv!=0) {qx=iv*bx;di-=nx*iv;} else qx=0;
		dijk=di+nx*(dj+oy*dk);
		create_periodic_image(di,dj,dk);

		// If mrs is bigger than the maximum distance to the block,
		// then we have to test all particles in the block for
		// intersections. Otherwise, we do additional checks and skip
		// those particles which can't possibly intersect the block.
		if(mrs>radius.cutoff(crs)) {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]-y;
				z1=p[dijk][sz*l+2]-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
			}
		} else {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]-y;
				z1=p[dijk][sz*l+2]-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(rs<mrs) {
					if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
				}
			}
		}
	} while(g<f);

	// If we reach here, we were unable to compute the entire cell using
	// the first part of the worklist. This section of the algorithm
	// continues the worklist, but it now starts preparing the mask that we
	// need if we end up going block by block. We do the same as before,
	// but we put a mark down on the mask for every block that's tested.
	// The worklist also contains information about which neighbors of each
	// block are not also on the worklist, and we start storing those
	// points in a list in case we have to go block by block.
	ci=nx;
	cj=ey;
	ck=ez;

	// Update the mask counter, and if it wraps around then reset the
	// whole mask; that will only happen once every 2^32 tries
	mv++;
	if(mv==0) {
		for(l=0;l<hxyz;l++) mask[l]=0;
		mv=1;
	}

	// Reset the block by block counters
	s_start=s_end=0;

	while(g<seq_length-1) {

		// At the intervals specified by count_list, we recompute the
		// maximum radius squared
		if(g==next_count) {
			mrs=c.max_radius_squared();
			if(list_index!=list_size) next_count=count_list[list_index++];
		}

		// If mrs is less than the minimum distance to any untested
		// block, then we are done
		if(mrs<radius.cutoff(radp[g])) return true;
		g++;

		// Load in a block off the worklist, permute it with the
		// symmetry mask, and decode its position. These are all
		// integer bit operations so they should run very fast.
		q=e[g];q^=m1;q+=m2;
		di=q&127;di-=64;
		dj=(q>>7)&127;dj-=64;
		dk=(q>>14)&127;dk-=64;

		// Compute the position in the mask of the current block. If
		// this lies outside the mask, then skip it. Otherwise, mark
		// it.
		ei=ci+di;
		ej=cj+dj;
		ek=ck+dk;
		if(ei<0) continue;else if(ei>=hx) continue;
		if(ej<0) continue;else if(ej>=hy) continue;
		if(ek<0) continue;else if(ek>=hz) continue;
		eijk=ei+hx*(ej+hy*ek);
		mask[eijk]=mv;

		// Call the compute_min_max_radius() function. This returns
		// true if the minimum distance to the block is bigger than the
		// current mrs, in which case we skip this block and move on.
		// Otherwise, it computes the maximum distance to the block and
		// returns it in crs.
		if(compute_min_max_radius(di,dj,dk,fx,fy,fz,gxs,gys,gzs,crs,mrs)) continue;

		// Now compute which region we are going to loop over, adding a
		// displacement for the periodic cases
		di+=i;dj+=j;dk+=k;
		iv=step_div(di,nx);if(iv!=0) {qx=iv*bx;di-=nx*iv;} else qx=0;
		dijk=di+nx*(dj+oy*dk);
		create_periodic_image(di,dj,dk);

		// If mrs is bigger than the maximum distance to the block,
		// then we have to test all particles in the block for
		// intersections. Otherwise, we do additional checks and skip
		// those particles which can't possibly intersect the block.
		if(mrs>radius.cutoff(crs)) {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]-y;
				z1=p[dijk][sz*l+2]-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
			}
		} else {
			for(l=0;l<co[dijk];l++) {
				x1=p[dijk][sz*l]+qx-x;
				y1=p[dijk][sz*l+1]-y;
				z1=p[dijk][sz*l+2]-z;
				rs=radius.scale(x1*x1+y1*y1+z1*z1,dijk,l);
				if(rs<mrs) {
					if(!c.nplane(x1,y1,z1,rs,id[dijk][l])) return false;
				}
			}
		}

		// If there might not be enough memory on the list for these
		// additions, then add more
		if(s_end+18>s_size) add_list_memory();

		// Test the parts of the worklist element which tell us what
		// neighbors of this block are not on the worklist. Store them
		// on the block list, and mark the mask.
		if((q&b2)==b2) {
			if(ei>0) if(mask[eijk-1]!=mv) {mask[eijk-1]=mv;sl[s_end++]=ei-1;sl[s_end++]=ej;sl[s_end++]=ek;}
			if((q&b1)==0) if(ei<hx-1) if(mask[eijk+1]!=mv) {mask[eijk+1]=mv;sl[s_end++]=ei+1;sl[s_end++]=ej;sl[s_end++]=ek;}
		} else if((q&b1)==b1) {if(ei<hx-1) if(mask[eijk+1]!=mv) {mask[eijk+1]=mv;sl[s_end++]=ei+1;sl[s_end++]=ej;sl[s_end++]=ek;}}
		if((q&b4)==b4) {if(ej>0) if(mask[eijk-hx]!=mv) {mask[eijk-hx]=mv;sl[s_end++]=ei;sl[s_end++]=ej-1;sl[s_end++]=ek;}
			if((q&b3)==0) if(ej<hy-1) if(mask[eijk+hx]!=mv) {mask[eijk+hx]=mv;sl[s_end++]=ei;sl[s_end++]=ej+1;sl[s_end++]=ek;}
		} else if((q&b3)==b3) {if(ej<hy-1) if(mask[eijk+hx]!=mv) {mask[eijk+hx]=mv;sl[s_end++]=ei;sl[s_end++]=ej+1;sl[s_end++]=ek;}}
		if((q&b6)==b6) {if(ek>0) if(mask[eijk-hxy]!=mv) {mask[eijk-hxy]=mv;sl[s_end++]=ei;sl[s_end++]=ej;sl[s_end++]=ek-1;}
			if((q&b5)==0) if(ek<hz-1) if(mask[eijk+hxy]!=mv) {mask[eijk+hxy]=mv;sl[s_end++]=ei;sl[s_end++]=ej;sl[s_end++]=ek+1;}
		} else if((q&b5)==b5) if(ek<hz-1) if(mask[eijk+hxy]!=mv) {mask[eijk+hxy]=mv;sl[s_end++]=ei;sl[s_end++]=ej;sl[s_end++]=ek+1;}
	}

	// Do a check to see if we've reached the radius cutoff
	if(mrs<radius.cutoff(radp[g])) return true;

	// Update the mask counter, and if it has wrapped around, then
	// reset the mask
	fx+=boxx*ci;fy+=boxy*cj;fz+=boxz*ck;

	// We were unable to completely compute the cell based on the blocks in
	// the worklist, so now we have to go block by block, reading in items
	// off the list
	while(s_start!=s_end) {

		// If we reached the end of the list memory loop back to the
		// start
		if(s_start==s_size) s_start=0;

		// Read in a block off the list, and compute the upper and lower
		// coordinates in each of the three dimensions
		di=sl[s_start++];dj=sl[s_start++];dk=sl[s_start++];
		xlo=di*boxx-fx;xhi=xlo+boxx;
		ylo=dj*boxy-fy;yhi=ylo+boxy;
		zlo=dk*boxz-fz;zhi=zlo+boxz;

		// Carry out plane tests to see if any particle in this block
		// could possibly intersect the cell
		if(di>ci) {
			if(dj>cj) {
				if(dk>ck) {
					if(corner_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if(edge_z_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if(corner_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xlo,yhi,zhi,xhi,ylo,zlo)) continue;
				} else {
					if(edge_z_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if(edge_y_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(edge_y_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if(face_x_test(c,xlo,ylo,zlo,yhi,zhi)) continue;
				}
			}
		} else if(di<ci) {
			if(dj>cj) {
				if(dk>ck) {
					if(corner_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xhi,ylo,zhi,xlo,yhi,zlo)) continue;
				} else {
					if(edge_z_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if(corner_test(c,xhi,yhi,zlo,xlo,ylo,zhi)) continue;
				} else if(dk<ck) {
					if(corner_test(c,xhi,yhi,zhi,xlo,ylo,zlo)) continue;
				} else {
					if(edge_z_test(c,xhi,yhi,zlo,xlo,ylo,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if(edge_y_test(c,xhi,ylo,zlo,xlo,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(edge_y_test(c,xhi,ylo,zhi,xlo,yhi,zlo)) continue;
				} else {
					if(face_x_test(c,xhi,ylo,zlo,yhi,zhi)) continue;
				}
			}
		} else {
			if(dj>cj) {
				if(dk>ck) {
					if(edge_x_test(c,xlo,ylo,zlo,xhi,yhi,zhi)) continue;
				} else if(dk<ck) {
					if(edge_x_test(c,xlo,ylo,zhi,xhi,yhi,zlo)) continue;
				} else {
					if(face_y_test(c,xlo,ylo,zlo,xhi,zhi)) continue;
				}
			} else if(dj<cj) {
				if(dk>ck) {
					if(edge_x_test(c,xlo,yhi,zlo,xhi,ylo,zhi)) continue;
				} else if(dk<ck) {
					if(edge_x_test(c,xlo,yhi,zhi,xhi,ylo,zlo)) continue;
				} else {
					if(face_y_test(c,xlo,yhi,zlo,xhi,zhi)) continue;
				}
			} else {
				if(dk>ck) {
					if(face_z_test(c,xlo,ylo,zlo,xhi,yhi)) continue;
				} else if(dk<ck) {
					if(face_z_test(c,xlo,ylo,zhi,xhi,yhi)) continue;
				} else voropp_fatal_error("Compute cell routine revisiting central block, which should never\nhappen.",VOROPP_INTERNAL_ERROR);
			}
		}

		// Now compute the region that we are going to test over, and
		// set a displacement vector for the periodic cases
		ei=i+di-nx;
		iv=step_div(ei,nx);if(iv!=0) {qx=iv*bx;ei-=nx*iv;} else qx=0;
		ej=j+dj-ey;
		ek=k+dk-ez;
		eijk=ei+nx*(ej+oy*ek);
		create_periodic_image(ei,ej,ek);

		// Loop over all the elements in the block to test for cuts. It
		// would be possible to exclude some of these cases by testing
		// against mrs, but this will probably not save time.
		for(l=0;l<co[eijk];l++) {
			x1=p[eijk][sz*l]+qx-x;
			y1=p[eijk][sz*l+1]-y;
			z1=p[eijk][sz*l+2]-z;
			rs=radius.scale(x1*x1+y1*y1+z1*z1,eijk,l);
			if(!c.nplane(x1,y1,z1,rs,id[eijk][l])) return false;
		}

		// If there's not much memory on the block list then add more
		if((s_start<=s_end?s_size-s_end+s_start:s_end-s_start)<18) add_list_memory();

		// Test the neighbors of the current block, and add them to the
		// block list if they haven't already been tested
		dijk=di+hx*(dj+hy*dk);
		if(di>0) if(mask[dijk-1]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-1]=mv;sl[s_end++]=di-1;sl[s_end++]=dj;sl[s_end++]=dk;}
		if(dj>0) if(mask[dijk-hx]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-hx]=mv;sl[s_end++]=di;sl[s_end++]=dj-1;sl[s_end++]=dk;}
		if(dk>0) if(mask[dijk-hxy]!=mv) {if(s_end==s_size) s_end=0;mask[dijk-hxy]=mv;sl[s_end++]=di;sl[s_end++]=dj;sl[s_end++]=dk-1;}
		if(di<hx-1) if(mask[dijk+1]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+1]=mv;sl[s_end++]=di+1;sl[s_end++]=dj;sl[s_end++]=dk;}
		if(dj<hy-1) if(mask[dijk+hx]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+hx]=mv;sl[s_end++]=di;sl[s_end++]=dj+1;sl[s_end++]=dk;}
		if(dk<hz-1) if(mask[dijk+hxy]!=mv) {if(s_end==s_size) s_end=0;mask[dijk+hxy]=mv;sl[s_end++]=di;sl[s_end++]=dj;sl[s_end++]=dk+1;}
	}

	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is at a corner.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (xl,yl,zl) the relative coordinates of the corner of the block
 *                       closest to the cell center.
 * \param[in] (xh,yh,zh) the relative coordinates of the corner of the block
 *                       furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::corner_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint zl,fpoint xh,fpoint yh,fpoint zh) {
	if(c.plane_intersects_guess(xh,yl,zl,radius.cutoff(xl*xh+yl*yl+zl*zl))) return false;
	if(c.plane_intersects(xh,yh,zl,radius.cutoff(xl*xh+yl*yh+zl*zl))) return false;
	if(c.plane_intersects(xl,yh,zl,radius.cutoff(xl*xl+yl*yh+zl*zl))) return false;
	if(c.plane_intersects(xl,yh,zh,radius.cutoff(xl*xl+yl*yh+zl*zh))) return false;
	if(c.plane_intersects(xl,yl,zh,radius.cutoff(xl*xl+yl*yl+zl*zh))) return false;
	if(c.plane_intersects(xh,yl,zh,radius.cutoff(xl*xh+yl*yl+zl*zh))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on an edge which points along the x
 * direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
 *                    block.
 * \param[in] (yl,zl) the relative y and z coordinates of the corner of the
 *                    block closest to the cell center.
 * \param[in] (yh,zh) the relative y and z coordinates of the corner of the
 *                    block furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::edge_x_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint zl,fpoint x1,fpoint yh,fpoint zh) {
	if(c.plane_intersects_guess(x0,yl,zh,radius.cutoff(yl*yl+zl*zh))) return false;
	if(c.plane_intersects(x1,yl,zh,radius.cutoff(yl*yl+zl*zh))) return false;
	if(c.plane_intersects(x1,yl,zl,radius.cutoff(yl*yl+zl*zl))) return false;
	if(c.plane_intersects(x0,yl,zl,radius.cutoff(yl*yl+zl*zl))) return false;
	if(c.plane_intersects(x0,yh,zl,radius.cutoff(yl*yh+zl*zl))) return false;
	if(c.plane_intersects(x1,yh,zl,radius.cutoff(yl*yh+zl*zl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on an edge which points along the y
 * direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
 *                    block.
 * \param[in] (xl,zl) the relative x and z coordinates of the corner of the
 *                    block closest to the cell center.
 * \param[in] (xh,zh) the relative x and z coordinates of the corner of the
 *                    block furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::edge_y_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint zl,fpoint xh,fpoint y1,fpoint zh) {
	if(c.plane_intersects_guess(xl,y0,zh,radius.cutoff(xl*xl+zl*zh))) return false;
	if(c.plane_intersects(xl,y1,zh,radius.cutoff(xl*xl+zl*zh))) return false;
	if(c.plane_intersects(xl,y1,zl,radius.cutoff(xl*xl+zl*zl))) return false;
	if(c.plane_intersects(xl,y0,zl,radius.cutoff(xl*xl+zl*zl))) return false;
	if(c.plane_intersects(xh,y0,zl,radius.cutoff(xl*xh+zl*zl))) return false;
	if(c.plane_intersects(xh,y1,zl,radius.cutoff(xl*xh+zl*zl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on an edge which points along the z
 * direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the block.
 * \param[in] (xl,yl) the relative x and y coordinates of the corner of the
 *                    block closest to the cell center.
 * \param[in] (xh,yh) the relative x and y coordinates of the corner of the
 *                    block furthest away from the cell center.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::edge_z_test(voronoicell_base<n_option> &c,fpoint xl,fpoint yl,fpoint z0,fpoint xh,fpoint yh,fpoint z1) {
	if(c.plane_intersects_guess(xl,yh,z0,radius.cutoff(xl*xl+yl*yh))) return false;
	if(c.plane_intersects(xl,yh,z1,radius.cutoff(xl*xl+yl*yh))) return false;
	if(c.plane_intersects(xl,yl,z1,radius.cutoff(xl*xl+yl*yl))) return false;
	if(c.plane_intersects(xl,yl,z0,radius.cutoff(xl*xl+yl*yl))) return false;
	if(c.plane_intersects(xh,yl,z0,radius.cutoff(xl*xh+yl*yl))) return false;
	if(c.plane_intersects(xh,yl,z1,radius.cutoff(xl*xh+yl*yl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on a face aligned with the x direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] xl the minimum distance from the cell center to the face.
 * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
 *                    block.
 * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
 *                    block.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::face_x_test(voronoicell_base<n_option> &c,fpoint xl,fpoint y0,fpoint z0,fpoint y1,fpoint z1) {
	if(c.plane_intersects_guess(xl,y0,z0,radius.cutoff(xl*xl))) return false;
	if(c.plane_intersects(xl,y0,z1,radius.cutoff(xl*xl))) return false;
	if(c.plane_intersects(xl,y1,z1,radius.cutoff(xl*xl))) return false;
	if(c.plane_intersects(xl,y1,z0,radius.cutoff(xl*xl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on a face aligned with the y direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] yl the minimum distance from the cell center to the face.
 * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
 *                    block.
 * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
 *                    block.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::face_y_test(voronoicell_base<n_option> &c,fpoint x0,fpoint yl,fpoint z0,fpoint x1,fpoint z1) {
	if(c.plane_intersects_guess(x0,yl,z0,radius.cutoff(yl*yl))) return false;
	if(c.plane_intersects(x0,yl,z1,radius.cutoff(yl*yl))) return false;
	if(c.plane_intersects(x1,yl,z1,radius.cutoff(yl*yl))) return false;
	if(c.plane_intersects(x1,yl,z0,radius.cutoff(yl*yl))) return false;
	return true;
}

/** This function checks to see whether a particular block can possibly have
 * any intersection with a Voronoi cell, for the case when the closest point
 * from the cell center to the block is on a face aligned with the z direction.
 * \param[in,out] c a reference to a Voronoi cell.
 * \param[in] zl the minimum distance from the cell center to the face.
 * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
 *                    block.
 * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
 *                    block.
 * \return False if the block may intersect, true if does not. */
template<class r_option>
template<class n_option>
inline bool container_periodic_base<r_option>::face_z_test(voronoicell_base<n_option> &c,fpoint x0,fpoint y0,fpoint zl,fpoint x1,fpoint y1) {
	if(c.plane_intersects_guess(x0,y0,zl,radius.cutoff(zl*zl))) return false;
	if(c.plane_intersects(x0,y1,zl,radius.cutoff(zl*zl))) return false;
	if(c.plane_intersects(x1,y1,zl,radius.cutoff(zl*zl))) return false;
	if(c.plane_intersects(x1,y0,zl,radius.cutoff(zl*zl))) return false;
	return true;
}

template<class r_option>
inline void container_periodic_base<r_option>::create_periodic_image(int di,int dj,int dk) {
	if(di<0||di>=nx||dj<0||dj>=oy||dk<0||dk>=oz) 
		voropp_fatal_error("Constructing periodic image for nonexistent point",VOROPP_INTERNAL_ERROR);
	if(dk>=ez&&dk<wz) {
		if(dj<ey||dj>=wy) create_side_image(di,dj,dk); 
	} else create_vertical_image(di,dj,dk);
}

template<class r_option>
void container_periodic_base<r_option>::create_side_image(int di,int dj,int dk) {
	fpoint boxx=bx/nx;
	int l,dijk=di+nx*(dj+oy*dk),odijk,ima=step_div(dj-ey,ny);
	int qua=di+step_int(-ima*bxy*xsp),quadiv=step_div(qua,nx);
	int fi=qua-quadiv*nx,fijk=fi+nx*(dj-ima*ny+oy*dk);
	fpoint dis=ima*bxy+quadiv*bx,switchx=di*boxx-ima*bxy-quadiv*bx,adis;

	if((img[dijk]&1)==0) {
		if(di>0) {
			odijk=dijk-1;adis=dis;
		} else {
			odijk=dijk+nx-1;adis=dis+bx;
		}
		img[odijk]|=2;
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][sz*l]>switchx) quick_put(dijk,fijk,l,dis,by*ima,0);
			else quick_put(odijk,fijk,l,adis,by*ima,0);
		}
	}
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
			if(p[fijk][sz*l]<switchx) quick_put(dijk,fijk,l,dis,by*ima,0);
			else quick_put(odijk,fijk,l,adis,by*ima,0);
		}
	}
	img[dijk]=3;
}

template<class r_option>
inline void container_periodic_base<r_option>::quick_put(int reg,int fijk,int l,fpoint dx,fpoint dy,fpoint dz) {
	if(co[reg]==mem[reg]) add_particle_memory(reg);
	p[reg][sz*co[reg]]=p[fijk][sz*l]+dx;
	p[reg][sz*co[reg]+1]=p[fijk][sz*l+1]+dy;
	p[reg][sz*co[reg]+2]=p[fijk][sz*l+2]+dz;
	if(sz==4) p[reg][sz*co[reg]+3]=p[fijk][sz*l+3];
	id[reg][co[reg]++]=id[fijk][l];
}

template<class r_option>
void container_periodic_base<r_option>::create_vertical_image(int di,int dj,int dk) {
	fpoint boxx=bx/nx,boxy=by/ny;
	int l,dijk=di+nx*(dj+oy*dk),dijkl,dijkr,ima=step_div(dk-ez,nz);
	int qj=dj+step_int(-ima*byz*ysp),qjdiv=step_div(qj-ey,ny);
	int qi=di+step_int((-ima*bxz-qjdiv*bxy)*xsp),qidiv=step_div(qi,nx);
	int fi=qi-qidiv*nx,fj=qj-qjdiv*ny,fijk=fi+nx*(fj+oy*(dk-ima*nz)),fijk2;
	fpoint disy=ima*byz+qjdiv*by,switchy=(dj-ey)*boxy-ima*byz-qjdiv*by;
	fpoint disx=ima*bxz+qjdiv*bxy+qidiv*bx,switchx=di*boxx-ima*bxz-qjdiv*bxy-qidiv*bx,switchx2,disxl,disxr,disx2,disxr2;

//	printf("%d %d %d %d ima=%d qi=%d qidiv=%d qj=%d qjdiv=%d fi=%d fj=%d fijk=%d\n",di,dj,dk,dijk,ima,qi,qidiv,qj,qjdiv,fi,fj,fijk);

	if(di==0) {
		dijkl=dijk+nx-1;disxl=disx+bx;
	} else {
		dijkl=dijk-1;disxl=disx;
	}

	if(di==nx-1) {
		dijkr=dijk-nx+1;disxr=disx-bx;
	} else {
		dijkr=dijk+1;disxr=disx;
	}

	bool y_exist=dj!=0;

//	printf("%d %d %d\n",dijkl,dijk,dijkr);
//	printf("%f %f %f %f %f\n",disxl,disx,disy,switchx,switchy);
	
	if((img[dijk]&1)==0) {
		img[dijkl]|=2;
		if(y_exist) {
			img[dijkl-nx]|=8;
			img[dijk-nx]|=4;
		}
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][sz*l+1]>switchy) {
				if(p[fijk][sz*l]>switchx) quick_put(dijk,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl,fijk,l,disxl,disy,bz*ima);
			} else {
				if(!y_exist) continue;
				if(p[fijk][sz*l]>switchx) quick_put(dijk-nx,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl-nx,fijk,l,disxl,disy,bz*ima);
			}
		}
	}
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
			if(p[fijk2][sz*l+1]>switchy) {
				if(p[fijk2][sz*l]>switchx2) quick_put(dijkr,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk,fijk2,l,disx2,disy,bz*ima);
			} else {
				if(!y_exist) continue;
				if(p[fijk2][sz*l]>switchx2) quick_put(dijkr-nx,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk-nx,fijk2,l,disx2,disy,bz*ima);
			}
		}
	}



//Recomputation of some stuff
	if(fj==wy-1) {
		fijk+=nx*(1-ny)-fi;
		switchy+=(1-ny)*boxy;
		disy+=by;
	//	int qj=dj+step_int(-ima*byz*ysp),qjdiv=step_div(qj,ny);
		qi=di+step_int(-(ima*bxz+(qjdiv+1)*bxy)*xsp);
		int dqidiv=step_div(qi,nx)-qidiv;qidiv+=dqidiv;
		fi=qi-qidiv*nx;
		fijk+=fi;
		disx+=bxy+bx*dqidiv;
		disxl+=bxy+bx*dqidiv;
		disxr+=bxy+bx*dqidiv;
		switchx-=bxy+bx*dqidiv;
	//	switchx=di*boxx-ima*bxz+qjdiv*bxy-qidiv*bx;
	} else {
		fijk+=nx;switchy+=boxy;
	}

	y_exist=dj!=oy-1;
	
	if((img[dijk]&4)==0) {
		img[dijkl]|=8;
		if(y_exist) {
			img[dijkl+nx]|=2;
			img[dijk+nx]|=1;
		}
		for(l=0;l<co[fijk];l++) {
			if(p[fijk][sz*l+1]>switchy) {
				if(!y_exist) continue;
				if(p[fijk][sz*l]>switchx) quick_put(dijk+nx,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl+nx,fijk,l,disxl,disy,bz*ima);
			} else {
				if(p[fijk][sz*l]>switchx) quick_put(dijk,fijk,l,disx,disy,bz*ima);
				else quick_put(dijkl,fijk,l,disxl,disy,bz*ima);
			}
		}
	}
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
			if(p[fijk2][sz*l+1]>switchy) {
				if(!y_exist) continue;
				if(p[fijk2][sz*l]>switchx2) quick_put(dijkr+nx,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk+nx,fijk2,l,disx2,disy,bz*ima);
			} else {
				if(p[fijk2][sz*l]>switchx2) quick_put(dijkr,fijk2,l,disxr2,disy,bz*ima);
				else quick_put(dijk,fijk2,l,disx2,disy,bz*ima);
			}
		}
	}
	
	img[dijk]=15;
}

/** Creates a voropp_loop object, by setting the necessary constants about the
 * container geometry from a pointer to the current container class.
 * \param[in] q a pointer to the current container class. */
template<class r_option>
voropp_loop::voropp_loop(container_periodic_base<r_option> *q) : sx(q->bx), sy(q->by), sz(q->bz),
	xsp(q->xsp),ysp(q->ysp),zsp(q->zsp),
	nx(q->nx),ny(q->ny),nz(q->nz),nxyz(q->nxyz),
	ey(q->ey),ez(q->ez),oy(q->oy),oz(q->oz) {}

/** Initializes a voropp_loop object, by finding all blocks which are within a
 * given sphere. It calculates the index of the first block that needs to be
 * tested and sets the periodic displacement vector accordingly.
 * \param[in] (vx,vy,vz) the position vector of the center of the sphere.
 * \param[in] r the radius of the sphere.
 * \param[out] (px,py,pz) the periodic displacement vector for the first block
 *                        to be tested.
 * \return The index of the first block to be tested. */
inline int voropp_loop::init(fpoint vx,fpoint vy,fpoint vz,fpoint r,fpoint &px) {
	ai=step_int((vx-r)*xsp);
	bi=step_int((vx+r)*xsp);
	aj=step_int((vy-r)*ysp)+ey;if(aj<0) aj=0;
	bj=step_int((vy+r)*ysp)+ey;if(bj>=oy) bj=oy-1;
	ak=step_int((vz-r)*zsp)+ez;if(ak<0) ak=0;
	bk=step_int((vz+r)*zsp)+ez;if(bk>=oz) bk=oz-1;
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(oy+aj-bj)+inc1;
	inc1+=nx;
	s=aip+nx*(aj+oy*ak);
	return s;
}

/** Initializes a voropp_loop object, by finding all blocks which overlap a given
 * rectangular box. It calculates the index of the first block that needs to be
 * tested and sets the periodic displacement vector (px,py,pz) accordingly.
 * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
 * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
 * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box.
 * \param[out] (px,py,pz) the periodic displacement vector for the first block
 *                        to be tested.
 * \return The index of the first block to be tested. */
inline int voropp_loop::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax,fpoint &px) {
	ai=step_int(xmin*xsp);
	bi=step_int(xmax*xsp);
	aj=step_int(ymin*ysp);if(aj<0) aj=0;
	bj=step_int(ymax*ysp);if(bj>=oy) bj=oy-1;
	ak=step_int(zmin*zsp);if(ak<0) ak=0;
	bk=step_int(zmax*zsp);if(bk>=oz) bk=oz-1;
	i=ai;j=aj;k=ak;
	aip=ip=step_mod(i,nx);apx=px=step_div(i,nx)*sx;
	inc1=aip-step_mod(bi,nx);
	inc2=nx*(oy+aj-bj)+inc1;
	inc1+=nx;
	s=aip+nx*(aj+oy*ak);
	return s;
}

/** Returns the next block to be tested in a loop, and updates the periodicity
 * vector if necessary.
 * \param[in,out] (px,py,pz) the current block on entering the function, which
 *                           is updated to the next block on exiting the
 *                           function. */
inline int voropp_loop::inc(fpoint &px) {
	if(i<bi) {
		i++;
		if(ip<nx-1) {ip++;s++;} else {ip=0;s+=1-nx;px+=sx;}
		return s;
	} else if(j<bj) {
		i=ai;ip=aip;px=apx;j++;
		s+=inc1;
		return s;
	} else if(k<bk) {
		i=ai;ip=aip;j=aj;px=apx;k++;
		s+=inc2;
		return s;
	} else return -1;
}

/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1). */
inline int voropp_loop::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
}

/** Custom modulo function, that gives consistent stepping for negative
 * numbers. */
inline int voropp_loop::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

/** Custom integer division function, that gives consistent stepping for
 * negative numbers. */
inline int voropp_loop::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}

/** Adds a wall to the container.
 * \param[in] w a wall object to be added.*/
template<class r_option>
void container_periodic_base<r_option>::add_wall(wall& w) {
	if(wall_number==current_wall_size) {
		current_wall_size*=2;
		if(current_wall_size>max_wall_size)
			voropp_fatal_error("Wall memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
		wall **pwall;
		pwall=new wall*[current_wall_size];
		for(int i=0;i<wall_number;i++) pwall[i]=walls[i];
		delete [] walls;walls=pwall;
	}
	walls[wall_number++]=&w;
}

/** Sets the radius of the jth particle in region i to r, and updates the
 * maximum particle radius.
 * \param[in] i the region of the particle to consider.
 * \param[in] j the number of the particle within the region.
 * \param[in] r the radius to set. */
inline void radius_poly::store_radius(int i,int j,fpoint r) {
	cc->p[i][4*j+3]=r;
	if(r>max_radius) max_radius=r;
}

/** Clears the stored maximum radius. */
inline void radius_poly::clear_max() {
	max_radius=0;
}

/** Imports a list of particles from an input stream for the monodisperse case
 * where no radius information is expected.
 * \param[in] is an input stream to read from. */
inline void radius_mono::import(istream &is) {
	int n;fpoint x,y,z;
	is >> n >> x >> y >> z;
	while(!is.eof()) {
		cc->put(n,x,y,z);
		is >> n >> x >> y >> z;
	}
}

/** Imports a list of particles from an input stream for the polydisperse case,
 * where both positions and particle radii are both stored.
 * \param[in] is an input stream to read from. */
inline void radius_poly::import(istream &is) {
	int n;fpoint x,y,z,r;
	is >> n >> x >> y >> z >> r;
	while(!is.eof()) {
		cc->put(n,x,y,z,r);
		is >> n >> x >> y >> z >> r;
	}
}

/** Initializes the radius_poly class for a new Voronoi cell calculation, by
 * computing the radial cut-off value, based on the current particle's radius
 * and the maximum radius of any particle in the packing.
 * \param[in] ijk the region to consider.
 * \param[in] s the number of the particle within the region. */
inline void radius_poly::init(int ijk,int s) {
	crad=cc->p[ijk][4*s+3];
	mul=1+(crad*crad-max_radius*max_radius)/((max_radius+crad)*(max_radius+crad));
	crad*=crad;
}

/** This routine is called when deciding when to terminate the computation of a
 * Voronoi cell. For the Voronoi radical tessellation for a polydisperse case,
 * this routine multiplies the cutoff value by the scaling factor that was
 * precomputed in the init() routine.
 * \param[in] lrs a cutoff radius for the cell computation.
 * \return The value scaled by the factor mul. */
inline fpoint radius_poly::cutoff(fpoint lrs) {
	return mul*lrs;
}

/** This routine is called when deciding when to terminate the computation of a
 * Voronoi cell. For the monodisperse case, this routine just returns the same
 * value that is passed to it.
 * \param[in] lrs a cutoff radius for the cell computation.
 * \return The same value passed to it. */
inline fpoint radius_mono::cutoff(fpoint lrs) {
	return lrs;
}

/** Prints the radius of particle, by just supplying a generic value of "s".
 * \param[in] os the output stream to write to.
 * \param[in] l the region to consider.
 * \param[in] c the number of the particle within the region. */
inline void radius_mono::rad(ostream &os,int l,int c) {
	os << "s";
}

/** Prints the radius of a particle to an open output stream.
 * \param[in] os the output stream to write to.
 * \param[in] l the region to consider.
 * \param[in] c the number of the particle within the region. */
inline void radius_poly::rad(ostream &os,int l,int c) {
	os << cc->p[l][4*c+3];
}

/** Returns the scaled volume of a particle, which is always set to 0.125 for
 * the monodisperse case where particles are taken to have unit diameter.
 * \param[in] ijk the region to consider.
 * \param[in] s the number of the particle within the region.
 * \return The cube of the radius of the particle, which is 0.125 in this case.
 */
inline fpoint radius_mono::volume(int ijk,int s) {
	return 0.125;
}

/** Returns the scaled volume of a particle.
 * \param[in] ijk the region to consider.
 * \param[in] s the number of the particle within the region.
 * \return The cube of the radius of the particle. */
inline fpoint radius_poly::volume(int ijk,int s) {
	fpoint a=cc->p[ijk][4*s+3];
	return a*a*a;
}

/** Scales the position of a plane according to the relative sizes
 * of the particle radii.
 * \param[in] rs the distance between the Voronoi cell and the cutting plane.
 * \param[in] t the region to consider
 * \param[in] q the number of the particle within the region.
 * \return The scaled position. */
inline fpoint radius_poly::scale(fpoint rs,int t,int q) {
	return rs+crad-cc->p[t][4*q+3]*cc->p[t][4*q+3];
}

/** Applies a blank scaling to the position of a cutting plane.
 * \param[in] rs the distance between the Voronoi cell and the cutting plane.
 * \param[in] t the region to consider
 * \param[in] q the number of the particle within the region.
 * \return The scaled position, which for this case, is equal to rs. */
inline fpoint radius_mono::scale(fpoint rs,int t,int q) {
	return rs;
}

/** Prints the radius of a particle to an open file stream.
 * \param[in] os an open file stream.
 * \param[in] ijk the region to consider.
 * \param[in] q the number of the particle within the region.
 * \param[in] later A boolean value to determine whether or not to write a
 *                  space character before the first entry. */
inline void radius_poly::print(ostream &os,int ijk,int q,bool later) {
	if(later) os << " ";
	os << cc->p[ijk][4*q+3];
}

/** This function is called during container construction. The routine scans
 * all of the worklists in the wl[] array. For a given worklist of blocks
 * labeled \f$w_1\f$ to \f$w_n\f$, it computes a sequence \f$r_0\f$ to
 * \f$r_n\f$ so that $r_i$ is the minimum distance to all the blocks
 * \f$w_{j}\f$ where \f$j>i\f$ and all blocks outside the worklist. The values
 * of \f$r_n\f$ is calculated first, as the minimum distance to any block in
 * the shell surrounding the worklist. The \f$r_i\f$ are then computed in
 * reverse order by considering the distance to \f$w_{i+1}\f$. */
template<class r_option>
void container_periodic_base<r_option>::initialize_radii() {
	const unsigned int b1=1<<21,b2=1<<22,b3=1<<24,b4=1<<25,b5=1<<27,b6=1<<28;
	const fpoint xstep=bx/nx/fgrid;
	const fpoint ystep=by/ny/fgrid;
	const fpoint zstep=bz/nz/fgrid;
	int i,j,k,lx,ly,lz,l=0,q;
	unsigned int *e=const_cast<unsigned int*> (wl);
	fpoint xlo,ylo,zlo,xhi,yhi,zhi,minr;fpoint *radp=mrad;
	unsigned int f;
	for(zlo=0,zhi=zstep,lz=0;lz<hgrid;zlo=zhi,zhi+=zstep,lz++) {
		for(ylo=0,yhi=ystep,ly=0;ly<hgrid;ylo=yhi,yhi+=ystep,ly++) {
			for(xlo=0,xhi=xstep,lx=0;lx<hgrid;xlo=xhi,xhi+=xstep,l++,lx++) {
				minr=large_number;
				for(q=e[0]+1;q<seq_length;q++) {
					f=e[q];
					i=(f&127)-64;
					j=(f>>7&127)-64;
					k=(f>>14&127)-64;
					if((f&b2)==b2) {
						compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i-1,j,k);
						if((f&b1)==0) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i+1,j,k);
					} else if((f&b1)==b1) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i+1,j,k);
					if((f&b4)==b4) {
						compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j-1,k);
						if((f&b3)==0) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j+1,k);
					} else if((f&b3)==b3) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j+1,k);
					if((f&b6)==b6) {
						compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k-1);
						if((f&b5)==0) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k+1);
					} else if((f&b5)==b5) compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k+1);
				}
				q--;
				while(q>0) {
					radp[q]=minr;
					f=e[q];
					i=(f&127)-64;
					j=(f>>7&127)-64;
					k=(f>>14&127)-64;
					compute_minimum(minr,xlo,xhi,ylo,yhi,zlo,zhi,i,j,k);
					q--;
				}
				radp[0]=minr;
				e+=seq_length;
				radp+=seq_length;
			}
		}
	}
}

/** Computes the minimum distance from a subregion to a given block. If this distance
 * is smaller than the value of minr, then it passes
 * \param[in,out] minr a pointer to the current minimum distance. If the distance
 *                     computed in this function is smaller, then this distance is
 *                     set to the new one.
 * \param[out] (xlo,ylo,zlo) the lower coordinates of the subregion being
 *                           considered.
 * \param[out] (xhi,yhi,zhi) the upper coordinates of the subregion being
 *                           considered.
 * \param[in] (ti,tj,tk) the coordinates of the block. */
template<class r_option>
inline void container_periodic_base<r_option>::compute_minimum(fpoint &minr,fpoint &xlo,fpoint &xhi,fpoint &ylo,fpoint &yhi,fpoint &zlo,fpoint &zhi,int ti,int tj,int tk) {
	const fpoint boxx=bx/nx,boxy=by/ny,boxz=bz/nz;
	fpoint radsq,temp;
	if(ti>0) {temp=boxx*ti-xhi;radsq=temp*temp;}
	else if(ti<0) {temp=xlo-boxx*(1+ti);radsq=temp*temp;}
	else radsq=0;

	if(tj>0) {temp=boxy*tj-yhi;radsq+=temp*temp;}
	else if(tj<0) {temp=ylo-boxy*(1+tj);radsq+=temp*temp;}

	if(tk>0) {temp=boxz*tk-zhi;radsq+=temp*temp;}
	else if(tk<0) {temp=zlo-boxz*(1+tk);radsq+=temp*temp;}

	if(radsq<minr) minr=radsq;
}

/** This routine checks to see whether a point is within a particular distance
 * of a nearby region. If the point is within the distance of the region, then
 * the routine returns true, and computes the maximum distance from the point
 * to the region. Otherwise, the routine returns false.
 * \param[in] (di,dj,dk) the position of the nearby region to be tested,
 *                       relative to the region that the point is in.
 * \param[in] (fx,fy,fz) the displacement of the point within its region.
 * \param[in] (gxs,gys,gzs) the minimum squared distances from the point to the
 *                          sides of its region.
 * \param[out] crs a reference in which to return the maximum distance to the
 *                 region (only computed if the routine returns positive).
 * \param[in] mrs the distance to be tested.
 * \return False if the region is further away than mrs, true if the region in
 *         within mrs.*/
template<class r_option>
inline bool container_periodic_base<r_option>::compute_min_max_radius(int di,int dj,int dk,fpoint fx,fpoint fy,fpoint fz,fpoint gxs,fpoint gys,fpoint gzs,fpoint &crs,fpoint mrs) {
	fpoint xlo,ylo,zlo;
	const fpoint boxx=bx/nx,boxy=by/ny,boxz=bz/nz;
	const fpoint bxsq=boxx*boxx+boxy*boxy+boxz*boxz;
	if(di>0) {
		xlo=di*boxx-fx;
		crs=xlo*xlo;
		if(dj>0) {
			ylo=dj*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo+boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo+boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(2*xlo+boxx)+boxy*(2*ylo+boxy)+gzs;
			}
		} else if(dj<0) {
			ylo=(dj+1)*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo-boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(boxx*xlo-boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(2*xlo+boxx)+boxy*(-2*ylo+boxy)+gzs;
			}
		} else {
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=gys+boxx*(2*xlo+boxx);
		}
	} else if(di<0) {
		xlo=(di+1)*boxx-fx;
		crs=xlo*xlo;
		if(dj>0) {
			ylo=dj*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo+boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo+boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(-2*xlo+boxx)+boxy*(2*ylo+boxy)+gzs;
			}
		} else if(dj<0) {
			ylo=(dj+1)*boxy-fy;
			crs+=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo-boxy*ylo+boxz*zlo);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=bxsq+2*(-boxx*xlo-boxy*ylo-boxz*zlo);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=boxx*(-2*xlo+boxx)+boxy*(-2*ylo+boxy)+gzs;
			}
		} else {
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=gys+boxx*(-2*xlo+boxx);
		}
	} else {
		if(dj>0) {
			ylo=dj*boxy-fy;
			crs=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=boxy*(2*ylo+boxy);
		} else if(dj<0) {
			ylo=(dj+1)*boxy-fy;
			crs=ylo*ylo;
			if(dk>0) {
				zlo=dk*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;
				crs+=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				if(radius.cutoff(crs)>mrs) return true;
				crs+=gzs;
			}
			crs+=boxy*(-2*ylo+boxy);
		} else {
			if(dk>0) {
				zlo=dk*boxz-fz;crs=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(2*zlo+boxz);
			} else if(dk<0) {
				zlo=(dk+1)*boxz-fz;crs=zlo*zlo;if(radius.cutoff(crs)>mrs) return true;
				crs+=boxz*(-2*zlo+boxz);
			} else {
				crs=0;
				voropp_fatal_error("Min/max radius function called for central block, which should never\nhappen.",VOROPP_INTERNAL_ERROR);
			}
			crs+=gys;
		}
		crs+=gxs;
	}
	return false;
}


/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).*/
template<class r_option>
inline int container_periodic_base<r_option>::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
}

/** Custom mod function, that gives consistent stepping for negative numbers. */
template<class r_option>
inline int container_periodic_base<r_option>::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

/** Custom div function, that gives consistent stepping for negative numbers. */
template<class r_option>
inline int container_periodic_base<r_option>::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}

#include "worklist.cc"
