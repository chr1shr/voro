// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"
#include "container.hh"

facets_search::facets_search(container *q) : v(0),
	sx(q->bx-q->ax), sy(q->by-q->ay), sz(q->bz-q->az),
	xsp(q->xsp),ysp(q->ysp),zsp(q->zsp),
	ax(q->ax),ay(q->ay),az(q->az),
	nx(q->nx),ny(q->ny),nz(q->nz),nxy(q->nxy),nxyz(q->nxyz),
	xperiodic(q->xperiodic),yperiodic(q->yperiodic),zperiodic(q->zperiodic),
	hx(xperiodic?2*nx+1:nx),
	hy(yperiodic?2*ny+1:ny),
	hz(zperiodic?2*nz+1:nz),hxy(hx*hy),hxyz(hxy*hz) {
		m=new unsigned int[hxyz];
		for(int i=0;i<hxyz;i++) m[i]=0;
		s=new int[s_size];
	}
}

inline void facets_search::add_list_memory() {
	int i,j=0,*ps;
	s_size*=2;
	ps=new int[s_size];


	for(i=s_start;i<s_size;i++) ps[j++]=s[i];
	for(i=0;i<s_end;i++) ps[j++]=
	s_start=0;s_end=j
	delete [] s;s=ps;
}

inline int facets_search::init() {
	v++;if (v==0) {
		for(int i=0;i<hxyz;i++) m[i]=0;
		v++;
	}
}

inline int facets_search::inc(fpoint &px,fpoint &py,fpoint &pz) {
	



}
