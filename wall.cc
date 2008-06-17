// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#include "wall.hh"

wall_sphere::wall_sphere(fpoint ixc,fpoint iyc,fpoint izc,fpoint irc) : xc(ixc),yc(iyc),zc(izc),rc(irc) {}

bool wall_sphere::point_inside(fpoint x,fpoint y,fpoint z) {
	return (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)<rc*rc;
}

template<class n_option>
inline void wall_sphere::cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
	fpoint xd=x-xc,yd=y-yc,zd=z-zc,dq;
	dq=xd*xd+yd*yd+zd*zd;
	if (dq>1e-5) {
		dq=2*(sqrt(dq)*rc-dq);
		c.nplane(xd,yd,zd,dq,10);
	}
}

wall_plane::wall_plane(fpoint ixc,fpoint iyc,fpoint izc,fpoint iac) : xc(ixc),yc(iyc),zc(izc),ac(iac) {}

bool wall_plane::point_inside(fpoint x,fpoint y,fpoint z) {
	return x*xc+y*yc+z*zc<ac;
}
