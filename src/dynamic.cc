// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file dynamic.hh
 * \brief Header file for the dynamic extension classes, which add
 * functionality for a variety of dynamic particle motions. */

#include "dynamic.hh"

template<class r_option>
container_dynamic_base<r_option>::container_dynamic_base(fpoint xa,fpoint xb,fpoint ya,
		fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,const bool xper,const bool yper,
		const bool zper,int memi) : container_base<r_option> (xa,xb,ya,yb,za,zb,xn,yn,zn,xper,yper,zper,memi), v_inter(ve) {
	gh=new int[nxyz];
	ve=new fpoint*[nxyz];
	for(int i=0;i<nxyz;i++) ve[i]=new fpoint[3*memi];
}

template<class r_option>
container_dynamic_base<r_option>::~container_dynamic_base() {
	for(int i=0;i<nxyz;i++) delete [] ve[i];
	delete [] ve;
	delete [] gh;
}

template<class r_option>
void container_dynamic_base<r_option>::wall_diagnostic() {
	fpoint x,y,z,dx,dy,dz;
	int i,j,k,ijk=0,l,w;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(l=0;l<co[ijk];l++) {
			x=p[ijk][sz*l];y=p[ijk][sz*l+1];z=p[ijk][sz*l+2];
			cout << id[ijk][l] << " " << x << " " << y << " " << z;
			radius.print(cout,ijk,l);
			for(w=0;w<wall_number;w++) {
				walls[w]->min_distance(x,y,z,dx,dy,dz);
				cout << " " << dx << " " << dy << " " << dz << "\n";
			}
		}
	}	
}

template<class r_option>
int container_dynamic_base<r_option>::count(fpoint x,fpoint y,fpoint z,fpoint r) {
	fpoint xx,yy,zz,px,py,pz;
	int cou=0,s,q;
	voropp_loop l1(this);
	s=l1.init(x,y,z,r,px,py,pz);r*=r;
        do {
                for(q=0;q<co[s];q++) {
                        xx=p[s][sz*q]+px-x;yy=p[s][sz*q+1]+py-y;zz=p[s][sz*q+2]+pz-z;
			if (xx*xx+yy*yy+zz*zz<r) cou++;
		}
        } while((s=l1.inc(px,py,pz))!=-1);
	return cou;
}

template<class r_option>
void container_dynamic_base<r_option>::spot(fpoint cx,fpoint cy,fpoint cz,fpoint dx,fpoint dy,fpoint dz,fpoint rad) {
	velocity_constant vcl(dx,dy,dz);
	local_move(vcl,cx,cy,cz,rad);
}

template<class r_option>
void container_dynamic_base<r_option>::gauss_spot(fpoint cx,fpoint cy,fpoint cz,fpoint dx,fpoint dy,fpoint dz,fpoint dec,fpoint rad) {
	velocity_gaussian vcl(cx,cy,cz,dx,dy,dz,dec);
	local_move(vcl,cx,cy,cz,rad);
}

template<class r_option>
template<class v_class>
void container_dynamic_base<r_option>::local_move(v_class &vcl,fpoint cx,fpoint cy,fpoint cz,fpoint rad) {
	int l,ll,s,ni,nj,nk;
	fpoint x,y,z,ex,ey,ez,px,py,pz;
	voropp_loop l1(this);
	s=l1.init(cx,cy,cz,rad,px,py,pz);rad*=rad;
	do {gh[s]=co[s];} while ((s=l1.inc(px,py,pz))!=-1);

	s=l1.reset(px,py,pz);
	do {
		l=0;
		while(l<gh[s]) {
			x=p[s][sz*l]+px;ex=x-cx;
			y=p[s][sz*l+1]+py;ey=y-cy;
			z=p[s][sz*l+2]+pz;ez=z-cz;
			if(ex*ex+ey*ey+ez*ez>=rad) {l++;continue;}
			vcl.vel(s,l,x,y,z);

			ni=step_int((x-ax)*xsp);
			nj=step_int((y-ay)*ysp);
			nk=step_int((z-az)*zsp);

			if(ni==l1.ip&&nj==l1.jp&&nk==l1.kp) {
				p[s][sz*l]=x;
				p[s][sz*l+1]=y;
				p[s][sz*l+2]=z;
				l++;
			} else {
				if((xperiodic||(ni>=0&&ni<nx))
				 &&(yperiodic||(nj>=0&&nj<ny))
				 &&(zperiodic||(nk>=0&&nk<nz))) {
					if(xperiodic) {x-=step_div(ni,nx)*(bx-ax);ni=step_mod(ni,nx);}
					if(yperiodic) {y-=step_div(nj,ny)*(by-ay);nj=step_mod(nj,ny);}
					if(zperiodic) {z-=step_div(nk,nz)*(bz-az);nk=step_mod(nk,nz);}
					ni+=nx*(nj+ny*nk);
					if(co[ni]==mem[ni]) add_particle_memory(ni);
					id[ni][co[ni]]=id[s][l];
					p[ni][co[ni]*sz]=x;
					p[ni][co[ni]*sz+1]=y;
					p[ni][co[ni]*sz+2]=z;
					if(sz==4) p[ni][co[ni]*sz+3]=p[s][l*sz+3];
					co[ni]++;
				}
				co[s]--;
				id[s][l]=id[s][co[s]];
				for(ll=0;ll<sz;ll++) p[s][sz*l+ll]=p[s][sz*co[s]+ll];
				if(co[s]+1==gh[s]) {
					if (vcl.track_ve) {
						ve[s][3*l]=ve[s][3*co[s]];
						ve[s][3*l+1]=ve[s][3*co[s]+1];
						ve[s][3*l+2]=ve[s][3*co[s]+2];
					}
					gh[s]--;
				} else l++;
			}
		}
	} while ((s=l1.inc(px,py,pz))!=-1);
}

template<class r_option>
void container_dynamic_base<r_option>::relax(fpoint cx,fpoint cy,fpoint cz,fpoint rad,fpoint alpha) {
	int l,q,s1,s2;
	fpoint x,y,z,ox,oy,oz,px,py,pz,ex,ey,ez;
	fpoint dx,dy,dz,dd,rr,radsq=rad*rad;
	voropp_loop l1(this),l2(this);
	s1=l1.init(cx,cy,cz,rad,ox,oy,oz);
	do {
		for(l=0;l<3*co[s1];l++) ve[s1][l]=0;
	} while ((s1=l1.inc(ox,oy,oz))!=-1);

	s1=l1.reset(ox,oy,oz);
	do {
		for(l=0;l<co[s1];l++) {
			ex=p[s1][sz*l];if(xperiodic) ex+=ox;dx=ex-cx;
			ey=p[s1][sz*l+1];if(yperiodic) ey+=oy;dy=ey-cy;
			ez=p[s1][sz*l+2];if(zperiodic) ez+=oz;dz=ez-cz;
			dd=dx*dx+dy*dy+dz*dz;
			if(dd>=radsq) {l++;continue;}

			s2=l2.init(ex,ey,ez,1,px,py,pz);
			while(!l2.reached(l1.i,l1.j,l1.k)) {
				for(q=0;q<co[s2];q++) {
					x=p[s2][sz*q]-ex;if(xperiodic) x+=px;
					y=p[s2][sz*q+1]-ey;if(yperiodic) y+=py;
					z=p[s2][sz*q+2]-ez;if(zperiodic) z+=pz;
					rr=x*x+y*y+z*z;
					if (rr<1) {
						if(dd+rr+2*(dx*x+dy*y+dz*z)<radsq) {
							rr=alpha*(0.5-0.5/sqrt(rr));
							ve[s2][3*q]-=x*rr;
							ve[s2][3*q+1]-=y*rr;
							ve[s2][3*q+2]-=z*rr;
						} else rr=alpha*(1-1/sqrt(rr));
						ve[s1][3*l]+=x*rr;
						ve[s1][3*l+1]+=y*rr;
						ve[s1][3*l+2]+=z*rr;
					}
				}
				s2=l2.inc(px,py,pz);
			}
			for(q=0;q<l;q++) {
				x=p[s2][sz*q]-ex;if(xperiodic) x+=px;
				y=p[s2][sz*q+1]-ey;if(yperiodic) y+=py;
				z=p[s2][sz*q+2]-ez;if(zperiodic) z+=pz;
				rr=x*x+y*y+z*z;
				if (rr<1) {
					if(dd+rr+2*(dx*x+dy*y+dz*z)<radsq) {
						rr=alpha*(0.5-0.5/sqrt(rr));
						ve[s2][3*q]-=x*rr;
						ve[s2][3*q+1]-=y*rr;
						ve[s2][3*q+2]-=z*rr;
					} else rr=alpha*(1-1/sqrt(rr));
					ve[s1][3*l]+=x*rr;
					ve[s1][3*l+1]+=y*rr;
					ve[s1][3*l+2]+=z*rr;
				}
			}
			q++;
			while(q<co[s2]) {
				x=p[s2][sz*q]-ex;if(xperiodic) x+=px;
				y=p[s2][sz*q+1]-ey;if(yperiodic) y+=py;
				z=p[s2][sz*q+2]-ez;if(zperiodic) z+=pz;
				rr=x*x+y*y+z*z;
				if (rr<1) {
					if(dd+rr+2*(dx*x+dy*y+dz*z)>=radsq) {
						rr=alpha*(1-1/sqrt(rr));
						ve[s1][3*l]+=x*rr;
						ve[s1][3*l+1]+=y*rr;
						ve[s1][3*l+2]+=z*rr;
					}
				}
				q++;
			}
			while((s2=l2.inc(px,py,pz))!=-1) {
				for(q=0;q<co[s2];q++) {
					x=p[s2][sz*q]-ex;if(xperiodic) x+=px;
					y=p[s2][sz*q+1]-ey;if(yperiodic) y+=py;
					z=p[s2][sz*q+2]-ez;if(zperiodic) z+=pz;
					rr=x*x+y*y+z*z;
					if (rr<1) {
						if(dd+rr+2*(dx*x+dy*y+dz*z)>=radsq) {
							rr=alpha*(1-1/sqrt(rr));
							ve[s1][3*l]+=x*rr;
							ve[s1][3*l+1]+=y*rr;
							ve[s1][3*l+2]+=z*rr;
						}
					}
				}
			}
			wall_contribution(s1,l,ex,ey,ez,alpha);
		}
	} while ((s1=l1.inc(ox,oy,oz))!=-1);

	local_move(v_inter,cx,cy,cz,rad);
}

template<class r_option>
inline void container_dynamic_base<r_option>::wall_contribution(int s,int l,fpoint cx,fpoint cy,fpoint cz,fpoint alpha) {
	fpoint x,y,z,rr;
	for(int q=0;q<wall_number;q++) {
		walls[q]->min_distance(cx,cy,cz,x,y,z);
		rr=x*x+y*y+z*z;
		if (rr<0.25) {
			rr=0.5*alpha*(1-0.5/sqrt(rr));
			ve[s][3*l]+=x*rr;
			ve[s][3*l+1]+=y*rr;
			ve[s][3*l+2]+=z*rr;
		}
	}
#ifdef VOROPP_AUTO_X_WALL
	if(!xperiodic) {
		if(cx-ax<0.5) ve[s][3*l]+=alpha*(0.5+ax-cx);
		if(bx-cx<0.5) ve[s][3*l]+=alpha*(bx-cx-0.5);
	}
#endif
#ifdef VOROPP_AUTO_Y_WALL
	if(!yperiodic) {
		if(cy-ay<0.5) ve[s][3*l+1]+=alpha*(0.5+ay-cy);
		if(by-cy<0.5) ve[s][3*l+1]+=alpha*(by-cy-0.5);
	}
#endif
#ifdef VOROPP_AUTO_Z_WALL
	if(!zperiodic) {
		if(cz-az<0.5) ve[s][3*l+2]+=alpha*(0.5+az-cz);
		if(bz-cz<0.5) ve[s][3*l+2]+=alpha*(bz-cz-0.5);
	}
#endif	
}

template<class r_option>
inline int container_dynamic_base<r_option>::full_count() {
	int ijk=0,count=0;
	while(ijk<nxyz) count+=co[ijk++];
	return count;
}

template<class r_option>
template<class cond_class>
void container_dynamic_base<r_option>::neighbor_distribution(int *bb,fpoint dr,int max) {
	int i,j,k,ijk,l,s,q,ll;
	cond_class c1;
	voropp_loop l1(this);
	fpoint idr=1/dr,maxr=dr*max,maxrsq=maxr*maxr;
	fpoint x,y,z,px,py,pz,cx,cy,cz,rr;
	for(ijk=k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(l=0;l<co[ijk];l++) {
			cx=p[ijk][sz*l];
			cy=p[ijk][sz*l+1];
			cz=p[ijk][sz*l+2];
			if(c1.test(cx,cy,cz)) {
				s=l1.init(cx,cy,cz,maxr,px,py,pz);
				do {
					for(q=0;q<co[s];q++) {
						x=p[s][sz*q]+px-cx;
						y=p[s][sz*q+1]+py-cy;
						z=p[s][sz*q+2]+pz-cz;
						rr=x*x+y*y+z*z;
						if(rr<maxrsq) {ll=int(idr*sqrt(rr));if(ll<max) bb[ll]++;}
					}
				} while ((s=l1.inc(px,py,pz))!=-1);
			}
		}
	}
}


template<class r_option>
template<class cond_class>
fpoint container_dynamic_base<r_option>::packing_badness() {
	int i,j,k,ijk,l,s,q;
	cond_class c1;
	voropp_loop l1(this);
	int pcount=0;
	fpoint x,y,z,px,py,pz,cx,cy,cz,rr,badcount=0;
	for(ijk=k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(l=0;l<co[ijk];l++) {
			cx=p[ijk][sz*l];
			cy=p[ijk][sz*l+1];
			cz=p[ijk][sz*l+2];
			if(c1.test(cx,cy,cz)) {
				pcount++;
				s=l1.init(cx,cy,cz,1,px,py,pz);
				do {
					for(q=0;q<co[s];q++) {
						x=p[s][sz*q]+px-cx;
						y=p[s][sz*q+1]+py-cy;
						z=p[s][sz*q+2]+pz-cz;
						rr=x*x+y*y+z*z;
						if(rr<1&&rr>tolerance) badcount+=1-2*sqrt(rr)+rr;
					}
				} while ((s=l1.inc(px,py,pz))!=-1);
				wall_badness(cx,cy,cz,badcount);
			}
		}
	}
	return pcount==0?0:sqrt(badcount/pcount);
}

template<class r_option>
inline void container_dynamic_base<r_option>::wall_badness(fpoint cx,fpoint cy,fpoint cz,fpoint &badcount) {
	fpoint x,y,z,rr;
	for(int q=0;q<wall_number;q++) {
		walls[q]->min_distance(cx,cy,cz,x,y,z);
		rr=x*x+y*y+z*z;
		if (rr<0.25) badcount+=0.25-sqrt(rr)+rr;
	}
#ifdef VOROPP_AUTO_X_WALL
	if(!xperiodic) {
		if(cx-ax<0.5) badcount+=(0.5+ax-cx)*(0.5+ax-cx);
		if(bx-cx<0.5) badcount+=(bx-cx-0.5)*(bx-cx-0.5);
	}
#endif
#ifdef VOROPP_AUTO_Y_WALL
	if(!yperiodic) {
		if(cy-ay<0.5) badcount+=(0.5+ay-cy)*(0.5+ay-cy);
		if(by-cy<0.5) badcount+=(by-cy-0.5)*(by-cy-0.5);
	}
#endif
#ifdef VOROPP_AUTO_Z_WALL
	if(!zperiodic) {
		if(cz-az<0.5) badcount+=(0.5+az-cz)*(0.5+az-cz);
		if(bz-cz<0.5) badcount+=(bz-cz-0.5)*(bz-cz-0.5);
	}
#endif	
}

#ifdef YEAST_ROUTINES
/** An overloaded version of the draw_yeast_pov() routine, that outputs
 * the particle positions to a file.
 * \param[in] filename the file to write to. */
template<class r_option>
void container_dynamic_base<r_option>::draw_yeast_pov(const char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	draw_yeast_pov(os);
	os.close();
}

/** Dumps all the particle positions in the POV-Ray format. */
template<class r_option>
void container_dynamic_base<r_option>::draw_yeast_pov(ostream &os) {
	int l,c;
	for(l=0;l<nxyz;l++) {
		for(c=0;c<co[l];c++) {
			os << "// id " << id[l][c] << "\n";
			os << "sphere{<" << p[l][sz*c] << "," << p[l][sz*c+1] << ","
			   << p[l][sz*c+2] << ">,";
			radius.rad(os,l,c);
			if(id[l][c]<stickycut) {
				os << " pigment{rgb <0.9,0.3,0.6>} finish{reflection 0.15 ambient 0.4 specular 0.3}}\n";
			} else {
				os << " pigment{rgb <0.9,0.85,0.35>} finish{reflection 0.15 specular 0.3 ambient 0.42}}\n";
			}
		}
	}
}

template<class r_option>
void container_dynamic_base<r_option>::stick(fpoint alpha) {
	int i,j,k,ijk,l,q,s;
	fpoint bigrad,smallrad,currad;
	fpoint bigfac,smallfac,curfac;
	fpoint x,y,z,px,py,pz,cx,cy,cz,rr;
	voropp_loop l1(this);

	for(ijk=0;ijk<nxyz;ijk++) for(l=0;l<3*co[ijk];l++) ve[ijk][l]=0;
	
	for(ijk=k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(l=0;l<co[ijk];l++) {
			cx=p[ijk][sz*l];
			cy=p[ijk][sz*l+1];
			cz=p[ijk][sz*l+2];
			if(id[ijk][l]<stickycut) {
				bigrad=1.6;smallrad=1.2;
				bigfac=1/0.8;smallfac=1/0.4;
			} else {
				bigrad=1.2;smallrad=0.799999999;
				bigfac=1/0.4;smallfac=0;
			}
			s=l1.init(cx,cy,cz,bigrad,px,py,pz);
			while(!l1.reached(i,j,k)) {
				for(q=0;q<co[s];q++) {
					x=p[s][sz*q]+px-cx;
					y=p[s][sz*q+1]+py-cy;
					z=p[s][sz*q+2]+pz-cz;
					rr=x*x+y*y+z*z;
					if(id[s][q]<stickycut) {
						currad=bigrad;curfac=bigfac;
					} else {
						currad=smallrad;curfac=smallfac;
					}
					if (rr<currad*currad) {
						if(rr<0.8*0.8) rr=0.5*alpha*(1-0.8/sqrt(rr));
						else {rr=sqrt(rr);rr=0.5*alpha*(rr-0.8)*(1-(rr-0.8)*curfac*(2-(rr-0.8)*curfac))/rr;}
						ve[ijk][3*l]+=x*rr;
						ve[ijk][3*l+1]+=y*rr;
						ve[ijk][3*l+2]+=z*rr;
						ve[s][3*q]-=x*rr;
						ve[s][3*q+1]-=y*rr;
						ve[s][3*q+2]-=z*rr;
					}
				}
				s=l1.inc(px,py,pz);
			}
			if (s!=ijk) throw fatal_error("s and ijk should be aligned and they're not");
			for(q=0;q<l;q++) {
				x=p[s][sz*q]+px-cx;
				y=p[s][sz*q+1]+py-cy;
				z=p[s][sz*q+2]+pz-cz;
				rr=x*x+y*y+z*z;
				if(id[s][q]<stickycut) {
					currad=bigrad;curfac=bigfac;
				} else {
					currad=smallrad;curfac=smallfac;
				}
				if (rr<currad*currad) {
					if(rr<0.8*0.8) rr=0.5*alpha*(1-0.8/sqrt(rr));
					else {rr=sqrt(rr);rr=0.5*alpha*(rr-0.8)*(1-(rr-0.8)*curfac*(2-(rr-0.8)*curfac))/rr;}
					ve[ijk][3*l]+=x*rr;
					ve[ijk][3*l+1]+=y*rr;
					ve[ijk][3*l+2]+=z*rr;
					ve[s][3*q]-=x*rr;
					ve[s][3*q+1]-=y*rr;
					ve[s][3*q+2]-=z*rr;
				}
			}
			wall_contribution(ijk,l,cx,cy,cz,alpha);
		}
	}

	move(v_inter);
}
#endif

template<class r_option>
inline void container_dynamic_base<r_option>::clear_velocities() {
	for(int l,ijk=0;ijk<nxyz;ijk++) for(l=0;l<3*co[ijk];l++) ve[ijk][l]=0;
}

template<class r_option>
inline void container_dynamic_base<r_option>::damp_velocities(fpoint damp) {
	for(int l,ijk=0;ijk<nxyz;ijk++) for(l=0;l<3*co[ijk];l++) ve[ijk][l]*=damp;
}

template<class r_option>
void container_dynamic_base<r_option>::full_relax(fpoint alpha) {
	int i,j,k,ijk,l,q,s;
	fpoint x,y,z,px,py,pz,cx,cy,cz,rr;
	voropp_loop l1(this);

	clear_velocities();

	for(ijk=k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		for(l=0;l<co[ijk];l++) {
			cx=p[ijk][sz*l];
			cy=p[ijk][sz*l+1];
			cz=p[ijk][sz*l+2];
			s=l1.init(cx,cy,cz,1,px,py,pz);
			while(!l1.reached(i,j,k)) {
				for(q=0;q<co[s];q++) {
					x=p[s][sz*q]+px-cx;
					y=p[s][sz*q+1]+py-cy;
					z=p[s][sz*q+2]+pz-cz;
					rr=x*x+y*y+z*z;
					if (rr<1) {
						rr=alpha*(0.5-0.5/sqrt(rr));
						ve[ijk][3*l]+=x*rr;
						ve[ijk][3*l+1]+=y*rr;
						ve[ijk][3*l+2]+=z*rr;
						ve[s][3*q]-=x*rr;
						ve[s][3*q+1]-=y*rr;
						ve[s][3*q+2]-=z*rr;
					}
				}
				s=l1.inc(px,py,pz);
			}
			if (s!=ijk) throw fatal_error("s and ijk should be aligned and they're not");
			for(q=0;q<l;q++) {
				x=p[s][sz*q]+px-cx;
				y=p[s][sz*q+1]+py-cy;
				z=p[s][sz*q+2]+pz-cz;
				rr=x*x+y*y+z*z;
				if (rr<1) {
					rr=alpha*(0.5-0.5/sqrt(rr));
					ve[ijk][3*l]+=x*rr;
					ve[ijk][3*l+1]+=y*rr;
					ve[ijk][3*l+2]+=z*rr;
					ve[s][3*q]-=x*rr;
					ve[s][3*q+1]-=y*rr;
					ve[s][3*q+2]-=z*rr;
				}
			}
			wall_contribution(ijk,l,cx,cy,cz,alpha);
		}
	}

	move(v_inter);
}

/** Increase memory for a particular region. */
template<class r_option>
void container_dynamic_base<r_option>::add_particle_memory(int i) {
	int *idp;fpoint *pp,*vep;
	int l,nmem=2*mem[i];
#if VOROPP_VERBOSE >=3
	cerr << "Particle memory in region " << i << " scaled up to " << nmem << endl;
#endif
	if(nmem>max_particle_memory) throw fatal_error("Absolute maximum memory allocation exceeded");
	idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	pp=new fpoint[sz*nmem];
	for(l=0;l<sz*co[i];l++) pp[l]=p[i][l];
	vep=new fpoint[3*nmem];
	for(l=0;l<3*co[i];l++) vep[l]=ve[i][l];
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
	delete [] ve[i];ve[i]=vep;
}

/** Custom int function, that gives consistent stepping for negative numbers.
 * With normal int, we have (-1.5,-0.5,0.5,1.5) -> (-1,0,0,1).
 * With this routine, we have (-1.5,-0.5,0.5,1.5) -> (-2,-1,0,1).*/
template<class r_option>
inline int container_dynamic_base<r_option>::step_int(fpoint a) {
	return a<0?int(a)-1:int(a);
}

/** Custom mod function, that gives consistent stepping for negative numbers. */
template<class r_option>
inline int container_dynamic_base<r_option>::step_mod(int a,int b) {
	return a>=0?a%b:b-1-(b-1-a)%b;
}

/** Custom div function, that gives consistent stepping for negative numbers. */
template<class r_option>
inline int container_dynamic_base<r_option>::step_div(int a,int b) {
	return a>=0?a/b:-1+(a+1)/b;
}

template<class r_option>
template<class v_class>
void container_dynamic_base<r_option>::move(v_class &vcl) {
	int i,j,k,ijk,l,ll,ni,nj,nk;
	fpoint x,y,z;

	// For each block, introduce a second counter that counts the number of
	// particles initially in this block, without including any that have
	// moved into the block
	for(ijk=0;ijk<nxyz;ijk++) gh[ijk]=co[ijk];

	// Loop over the blocks within the counter
	for(ijk=k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		l=0;

		// Scan over the particles that started within this block
		while(l<gh[ijk]) {

			// Get the particle's initial position
			x=p[ijk][sz*l];y=p[ijk][sz*l+1];z=p[ijk][sz*l+2];
			
			// Call an external routine to compute where this
			// particle should move to
			vcl.vel(ijk,l,x,y,z);
			
			// Calculate which block the particle's new position is
			// within
			ni=step_int((x-ax)*xsp);
			nj=step_int((y-ay)*ysp);
			nk=step_int((z-az)*zsp);
			
			// See if the particle is within the same block
			if(ni==i&&nj==j&&nk==k) {

				// If so, just update its position
				p[ijk][sz*l]=x;
				p[ijk][sz*l+1]=y;
				p[ijk][sz*l+2]=z;l++;
			} else {

				// Check whether the particle's new position is
				// within the container
				if((xperiodic||(ni>=0&&ni<nx))&&(yperiodic||(nj>=0&&nj<ny))&&(zperiodic||(nk>=0&&nk<nz))) {

					// If periodic boundary conditions are
					// being used, then remap the particle
					// back into the container
					if(xperiodic) {x-=step_div(ni,nx)*(bx-ax);ni=step_mod(ni,nx);}
					if(yperiodic) {y-=step_div(nj,ny)*(by-ay);nj=step_mod(nj,ny);}
					if(zperiodic) {z-=step_div(nk,nz)*(bz-az);nk=step_mod(nk,nz);}
					
					// Calculate the index the new block
					// where the particle is to be stored,
					// and allocate memory if necessary
					ni+=nx*(nj+ny*nk);
					if(co[ni]==mem[ni]) add_particle_memory(ni);

					// Add the particle to the new block
					id[ni][co[ni]]=id[ijk][l];
					p[ni][co[ni]*sz]=x;
					p[ni][co[ni]*sz+1]=y;
					p[ni][co[ni]*sz+2]=z;

					// This class can support a velocity
					// array. Move it too, if it is being
					// used.
					if(vcl.track_ve) {
						ve[ni][3*co[ni]]=ve[ijk][3*l];
						ve[ni][3*co[ni]+1]=ve[ijk][3*l+1];
						ve[ni][3*co[ni]+2]=ve[ijk][3*l+2];
					}

					// Copy particle radius information if
					// it being used, and add one to the
					// number of particles in the block
					if(sz==4) p[ni][co[ni]*sz+3]=p[ijk][l*sz+3];
					co[ni]++;
				}

				// Delete the particle from the original block,
				// by copying the final particle in the block's
				// memory on top of it
				co[ijk]--;
				id[ijk][l]=id[ijk][co[ijk]];
				for(ll=0;ll<sz;ll++) p[ijk][sz*l+ll]=p[ijk][sz*co[ijk]+ll];
				if (vcl.track_ve) {
					ve[ijk][3*l]=ve[ijk][3*co[ijk]];
					ve[ijk][3*l+1]=ve[ijk][3*co[ijk]+1];
					ve[ijk][3*l+2]=ve[ijk][3*co[ijk]+2];
				}

				// If the copied particle was one that had
				// initially started off in this block, then
				// prepare to check that one next. If it was a
				// particle that had moved into this block,
				// then skip it.
				if(co[ijk]+1==gh[ijk]) gh[ijk]--;
				else l++;
			}
		}
	}
}
