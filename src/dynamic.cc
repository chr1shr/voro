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
		bool zper,int memi) : container_base<r_option> (xa,xb,ya,yb,za,zb,xn,yn,zn,xper,yper,zper,memi) {
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
	fpoint px,py,pz;
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
void full_relax() {

}
