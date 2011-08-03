#include "unitcell.hh"
#include "cell.hh"

unitcell::unitcell(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_)
	: bx(bx_), bxy(bxy_), by(by_), bxz(bxz_), byz(byz_), bz(bz_) {
	int i,j,l=1;

	// Initialize the Voronoi cell to be a very large rectangular box
	const double ucx=max_unit_voro_shells*bx,ucy=max_unit_voro_shells*by,ucz=max_unit_voro_shells*bz;
	unit_voro.init(-ucx,ucx,-ucy,ucy,-ucz,ucz);

	// Repeatedly cut the cell by shells of periodic image particles
	while(l<2*max_unit_voro_shells) {

		// Check to see if any of the planes from the current shell
		// will cut the cell
		if(unit_voro_intersect(l)) {

			// If they do, apply the plane cuts from the current
			// shell
			unit_voro_apply(l,0,0);
			for(i=1;i<l;i++) {
				unit_voro_apply(l,i,0);
				unit_voro_apply(-l,i,0);
			}
			for(i=-l;i<=l;i++) unit_voro_apply(i,l,0);
			for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
				unit_voro_apply(l,j,i);
				unit_voro_apply(-j,l,i);
				unit_voro_apply(-l,-j,i);
				unit_voro_apply(j,-l,i);
			}
			for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) unit_voro_apply(i,j,l);
		} else {
			max_uv_z=max_uv_y=sqrt(unit_voro.max_radius_squared());
			return;	
		}
		l++;
	}

	// If the routine makes it here, then the unit cell still hasn't been
	// completely bounded by the plane cuts. Give the memory error code,
	// because this is mainly a case of hitting a safe limit, than any
	// inherent problem.
	voropp_fatal_error("Periodic cell computation failed",VOROPP_MEMORY_ERROR);
}

/** Applies a pair of opposing plane cuts from a periodic image point
 * to the unit Voronoi cell.
 * \param[in] (i,j,k) the index of the periodic image to consider. */
inline void unitcell::unit_voro_apply(int i,int j,int k) {
	double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
	unit_voro.plane(x,y,z);
	unit_voro.plane(-x,-y,-z);
}

bool unitcell::intersects_image(double dx,double dy,double dz,double &vol) {
	const double bxinv=1/bx,byinv=1/by,bzinv=1/bz,ivol=bxinv*byinv*bzinv;
	voronoicell c;
	c=unit_voro;
	if(!c.plane(0,0,bzinv,dx+1)) return false;
	if(!c.plane(0,0,-bzinv,-dx+1)) return false;
	if(!c.plane(0,byinv,-byz*byinv*bzinv,dy+1)) return false;
	if(!c.plane(0,-byinv,byz*byinv*bzinv,-dy+1)) return false;
	if(!c.plane(bxinv,-bxy*bxinv*byinv,(bxy*byz-by*bxz)*ivol,dz+1)) return false;
	if(!c.plane(-bxinv,bxy*bxinv*byinv,(-bxy*byz+by*bxz)*ivol,-dz+1)) return false;
	vol=c.volume()*ivol;
	return true;
}

/*bool unitcell::images(vector &vi,vector &vd) {
	for(i=0;i<ms2*ms2*ms2;i++) a[i]=true;

	*ac=false;
	queue<int> q;
	q.push(0);q.push(0);q.push(0);

	while(!q.empty()) {
		i=q.front();q.pop();
		j=q.front();q.pop();
		k=q.front();q.pop();
		if(intersects(i,j,k,vol)) {
			vi.push_back(i);
			vi.push_back(j);
			vi.push_back(k);
			vd.push_back(vol);
			printf("%d %d %d %g\n",i,j,k,vol);
			acp=ac+i+ms2*(j+ms2*k);
			if(k>-max_shells&&*(acp-ms2*ms2)) {q.push(i);q.push(j);q.push(k-1);*(acp-ms2*ms2)=false;}
			if(j>-max_shells&&*(acp-ms2)) {q.push(i);q.push(j-1);q.push(k);*(acp-ms2)=false;}
			if(i>-max_shells&&*(acp-1)) {q.push(i-1);q.push(j);q.push(k);*(acp-1)=false;}
			if(i<max_shells&&*(acp+1)) {q.push(i+1);q.push(j);q.push(k);*(acp+1)=false;}
			if(j<max_shells&&*(acp+ms2)) {q.push(i);q.push(j+1);q.push(k);*(acp+ms2)=false;}
			if(k<max_shells&&*(acp+ms2*ms2)) {q.push(i);q.push(j);q.push(k+1);*(acp+ms2*ms2)=false;}
		}
	}
}*/

/** Tests to see if a shell of periodic images could possibly cut the periodic
 * unit cell.
 * \param[in] l the index of the shell to consider. 
 * \return True if a point in the shell cuts the cell, false otherwise. */
bool unitcell::unit_voro_intersect(int l) {
	int i,j;
	if(unit_voro_test(l,0,0)) return true;
	for(i=1;i<l;i++) {
		if(unit_voro_test(l,i,0)) return true;
		if(unit_voro_test(-l,i,0)) return true;
	}
	for(i=-l;i<=l;i++) if(unit_voro_test(i,l,0)) return true;
	for(i=1;i<l;i++) for(j=-l+1;j<=l;j++) {
		if(unit_voro_test(l,j,i)) return true;
		if(unit_voro_test(-j,l,i)) return true;
		if(unit_voro_test(-l,-j,i)) return true;
		if(unit_voro_test(j,-l,i)) return true;
	}
	for(i=-l;i<=l;i++) for(j=-l;j<=l;j++) if(unit_voro_test(i,j,l)) return true;
	return false;
}

/** Tests to see if a plane cut from a particular periodic image will cut th
 * unit Voronoi cell.
 * \param[in] (i,j,k) the index of the periodic image to consider.
 * \return True if the image cuts the cell, false otherwise. */
inline bool unitcell::unit_voro_test(int i,int j,int k) {
	double x=i*bx+j*bxy+k*bxz,y=j*by+k*byz,z=k*bz;
	double rsq=x*x+y*y+z*z;
	return unit_voro.plane_intersects(x,y,z,rsq);
}

/** Draws the periodic domain in gnuplot format.
 * \param[in] filename the filename to write to. */
void unitcell::draw_domain_gnuplot(FILE *fp) {
	fprintf(fp,"0 0 0\n%g 0 0\n%g %g 0\n%g %g 0\n",bx,bx+bxy,by,bxy,by);
	fprintf(fp,"%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n",bxy+bxz,by+byz,bz,bx+bxy+bxz,by+byz,bz,bx+bxz,byz,bz,bxz,byz,bz);
	fprintf(fp,"0 0 0\n%g %g 0\n\n%g %g %g\n%g %g %g\n\n",bxy,by,bxz,byz,bz,bxy+bxz,by+byz,bz);
	fprintf(fp,"%g 0 0\n%g %g %g\n\n%g %g 0\n%g %g %g\n\n",bx,bx+bxz,byz,bz,bx+bxy,by,bx+bxy+bxz,by+byz,bz);
}

/** Draws the periodic domain in POV-Ray format.
 * \param[in] filename the filename to write to. */
void unitcell::draw_domain_pov(FILE *fp) {
	fprintf(fp,"cylinder{0,0,0>,<%g,0,0>,rr}\n"
		   "cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n",bx,bxy,by,bx+bxy,by);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",bxz,byz,bz,bx+bxz,byz,bz,bxy+bxz,by+byz,bz,bx+bxy+bxz,by+byz,bz);
	fprintf(fp,"cylinder{<0,0,0>,<%g,%g,0>,rr}\n"
		   "cylinder{<%g,0,0>,<%g,%g,0>,rr}\n",bxy,by,bx,bx+bxy,by);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",bxz,byz,bz,bxy+bxz,by+byz,bz,bx+bxz,byz,bz,bx+bxy+bxz,by+byz,bz);
	fprintf(fp,"cylinder{<0,0,0>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,0,0>,<%g,%g,%g>,rr}\n",bxz,byz,bz,bx,bx+bxz,byz,bz);
	fprintf(fp,"cylinder{<%g,%g,0>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,0>,<%g,%g,%g>,rr}\n",bxy,by,bxy+bxz,by+byz,bz,bx+bxy,by,bx+bxy+bxz,by+byz,bz);
	fprintf(fp,"sphere{<0,0,0>,rr}\nsphere{<%g,0,0>,rr}\n"
		   "sphere{<%g,%g,0>,rr}\nsphere{<%g,%g,0>,rr}\n",bx,bxy,by,bx+bxy,by);
	fprintf(fp,"sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
		   "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n",bxz,byz,bz,bx+bxz,byz,bz,bxy+bxz,by+byz,bz,bx+bxy+bxz,by+byz,bz);
}
