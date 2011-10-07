// Voro++, a cell-based Voronoi library
//
// Authors  : Chris H. Rycroft (LBL / UC Berkeley)
//            Cody Robert Dance (UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file cell_nc_2d.cc
 * \brief Function implementations for the non-convex 2D Voronoi classes. */

#include "cell_nc_2d.hh"

namespace voro {

void voronoicell_nonconvex_neighbor_2d::init(double xmin,double xmax,double ymin,double ymax) {
	nonconvex=false;
	init_base(xmin,xmax,ymin,ymax);
	*ne=-3;ne[1]=-2;ne[2]=-4;ne[3]=-1;
}

void voronoicell_nonconvex_base_2d::init_nonconvex_base(double xmin,double xmax,double ymin,double ymax,double wx0,double wy0,double wx1,double wy1) {
	nonconvex=true;
	xmin*=2;xmax*=2;ymin*=2;ymax*=2;
	int f0=face(xmin,xmax,ymin,ymax,wx0,wy0),
	    f1=face(xmin,xmax,ymin,ymax,wx1,wy1);

	*pts=0;pts[1]=0;
	pts[2]=wx0;pts[3]=wy0;p=4;
	if(f0!=f1||wx0*wy1<wx1*wy0) {
		do {
			if(f0>1) {
				if(f0==2) {pts[p++]=xmin;pts[p++]=ymin;}
				else {pts[p++]=xmax;pts[p++]=ymin;}
			} else {
				if(f0==0) {pts[p++]=xmax;pts[p++]=ymax;}
				else {pts[p++]=xmin;pts[p++]=ymax;}
			}
			f0++;f0&=3;
		} while(f0!=f1);
	}
	pts[p++]=wx1;pts[p++]=wy1;
	p>>=1;

	int i,*q=ed;
	*(q++)=1;*(q++)=p-1;
	for(i=1;i<p-1;i++) {*(q++)=i+1;*(q++)=i-1;}
	*(q++)=0;*q=p-2;
}

void voronoicell_nonconvex_neighbor_2d::init_nonconvex(double xmin,double xmax,double ymin,double ymax,double wx0,double wy0,double wx1,double wy1) {
	init_nonconvex_base(xmin,xmax,ymin,ymax,wx0,wy0,wx1,wy1);
	*ne=-5;
	for(int i=1;i<p-1;i++) ne[i]=-99;
	ne[p-1]=-5;
}

inline int voronoicell_nonconvex_base_2d::face(double xmin,double xmax,double ymin,double ymax,double &wx,double &wy) {
	if(wy>0) {
		if(xmin*wy>ymax*wx) {wy*=xmin/wx;wx=xmin;return 2;}
		if(xmax*wy>ymax*wx) {wx*=ymax/wy;wy=ymax;return 1;}
		wy*=xmax/wx;wx=xmax;return 0;
	}
	if(xmax*wy>ymin*wx) {wy*=xmax/wx;wx=xmax;return 0;}
	if(xmin*wy>ymin*wx) {wx*=ymin/wy;wy=ymin;return 3;}
	wy*=xmin/wx;wx=xmin;return 2;
}

}
