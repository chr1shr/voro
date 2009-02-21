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
		fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,
		bool zper,int memi) : container_base<r_option> (xa,xb,ya,yb,za,zb,xn,yn,zn,xper,yper,zper,memi), v_inter(ve) {
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
void container_dynamic_base<r_option>::full_relax(fpoint alpha) {
	int i,j,k,ijk,l,q,s;
	fpoint x,y,z,px,py,pz,cx,cy,cz,rr;
	voropp_loop l1(this);

	for(ijk=0;ijk<nxyz;ijk++) for(l=0;l<3*co[ijk];l++) ve[ijk][l]=0;
	
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
					if (x*x+y*y+z*z<1) {
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
				if (x*x+y*y+z*z<1) {
					rr=alpha*(0.5-0.5/sqrt(rr));
					ve[ijk][3*l]+=x*rr;
					ve[ijk][3*l+1]+=y*rr;
					ve[ijk][3*l+2]+=z*rr;
					ve[s][3*q]-=x*rr;
					ve[s][3*q+1]-=y*rr;
					ve[s][3*q+2]-=z*rr;
				}
			}
			for(q=0;q<wall_number;q++) {
				walls[q]->min_distance(cx,cy,cz,x,y,z);
				rr=x*x+y*y+z*z;
				if (x*x+y*y+z*z<0.25) {
					rr=0.5*alpha*(1-0.5/sqrt(rr));
					ve[ijk][3*l]+=x*rr;
					ve[ijk][3*l+1]+=y*rr;
					ve[ijk][3*l+2]+=z*rr;
				}
			}
		}
	}

	move(&v_inter);
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

template<class r_option>
template<class v_class>
void container_dynamic_base<r_option>::move(v_class *vcl) {
	int i,j,k,ijk,l,ll,ni,nj,nk;
	fpoint x,y,z;

	for(ijk=0;ijk<nxyz;ijk++) {
		gh[ijk]=co[ijk];
	}

	for(ijk=k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++,ijk++) {
		l=0;
		while(l<gh[ijk]) {
			x=p[ijk][sz*l];y=p[ijk][sz*l+1];z=p[ijk][sz*l+2];
			vcl->vel(ijk,l,x,y,z);
			ni=step_int((x-ax)*xsp);
			nj=step_int((y-ay)*ysp);
			nk=step_int((z-az)*zsp);
			if(ni==i&&nj==j&&nk==k) {
				p[ijk][sz*l]=x;
				p[ijk][sz*l+1]=y;
				p[ijk][sz*l+2]=z;l++;
			} else {
				if(xperiodic) ni=step_mod(ni,nx);
				if(yperiodic) nj=step_mod(nj,ny);
				if(zperiodic) nk=step_mod(nk,nz);
				if((xperiodic||(ni>=0&&ni<nx))&&(yperiodic||(nj>=0&&nj<ny))&&(zperiodic||(nk>=0&&nk<nz))) {
					ni+=nx*(nj+ny*nk);
					if(co[ni]==mem[ni]) add_particle_memory(ni);
					id[ni][co[ni]]=id[ijk][l];
					p[ni][co[ni]*sz]=x;
					p[ni][co[ni]*sz+1]=y;
					p[ni][co[ni]*sz+2]=z;
					if(sz==4) p[ni][co[ni]*sz+3]=p[ijk][l*sz+3];
					co[ni]++;
				}
				co[ijk]--;
				id[ijk][l]=id[ijk][co[ijk]];
				for(ll=0;ll<sz;ll++) p[ijk][sz*l+ll]=p[ijk][sz*co[ijk]+ll];
				if(co[ijk]+1==gh[ijk]) {
					gh[ijk]--;
					if (vcl->track_ve) {
						ve[ijk][3*l]=ve[ijk][3*co[ijk]];
						ve[ijk][3*l+1]=ve[ijk][3*co[ijk]+1];
						ve[ijk][3*l+2]=ve[ijk][3*co[ijk]+2];
					}
				} else l++;
			}
		}
	}
}
