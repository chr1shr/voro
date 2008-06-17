// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#ifndef FACETS_WALL_HH
#define FACETS_WALL_HH

struct wall_sphere : public wall {
	const fpoint xc,yc,zc,rc;
	wall_sphere(fpoint ixc,fpoint iyc,fpoint izc,fpoint irc);
	bool point_inside(fpoint x,fpoint y,fpoint z);
	template<class n_option>
	inline void cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
	void cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
	void cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
};

struct wall_plane : public wall {
	const fpoint xc,yc,zc,ac;
	wall_plane(fpoint ixc,fpoint iyc,fpoint izc,fpoint iac);
	virtual bool point_inside(fpoint x,fpoint y,fpoint z);
	template <class n_option>
	void cut_cell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
};

struct wall_open_cylinder : public wall {
	const fpoint xc,yz,zc,rc;
	wall_open_cylinder(fpoint ixc,fpoint iyc,fpoint izc,fpoint irc);
	virtual bool point_inside(fpoint x,fpoint y,fpoint z);
	template <class n_option>
	void cut_cell(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
};

#endif
