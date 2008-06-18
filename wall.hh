// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : October 19th 2007

#ifndef FACETS_WALL_HH
#define FACETS_WALL_HH

struct wall_sphere : public wall {
	const int w_id;
	const fpoint xc,yc,zc,rc;
	wall_sphere(fpoint ixc,fpoint iyc,fpoint izc,fpoint irc)
		: w_id(-99), xc(ixc), yc(iyc), zc(izc), rc(irc) {};
	wall_sphere(int iw_id,fpoint ixc,fpoint iyc,fpoint izc,fpoint irc)
		: w_id(iw_id), xc(ixc), yc(iyc), zc(izc), rc(irc) {};
	bool point_inside(fpoint x,fpoint y,fpoint z);
	template<class n_option>
	inline void cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
	void cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
	void cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
};

struct wall_plane : public wall {
	const int w_id; 
	const fpoint xc,yc,zc,ac;
	wall_plane(fpoint ixc,fpoint iyc,fpoint izc,fpoint iac)
		: w_id(-99), xc(ixc), yc(iyc), zc(izc), ac(iac) {};
	wall_plane(int iw_id,fpoint ixc,fpoint iyc,fpoint izc,fpoint iac)
		: w_id(iw_id), xc(ixc), yc(iyc), zc(izc), ac(iac) {};
	bool point_inside(fpoint x,fpoint y,fpoint z);
	template<class n_option>
	inline void cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
	void cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
	void cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
};

struct wall_cylinder : public wall {
	const int w_id;
	const fpoint xc,yc,zc,xa,ya,za,asi,rc;
	wall_cylinder(fpoint ixc,fpoint iyc,fpoint izc,fpoint ixa,fpoint iya,fpoint iza,fpoint irc)
		: w_id(-99), xc(ixc), yc(iyc), zc(izc), xa(ixa), ya(iya), za(iza),
		asi(1/(ixa*ixa+iya*iya+iza*iza)), rc(irc) {};
	wall_cylinder(int iw_id,fpoint ixc,fpoint iyc,fpoint izc,fpoint ixa,fpoint iya,fpoint iza,fpoint irc)
		: w_id(iw_id), xc(ixc), yc(iyc), zc(izc), xa(ixa), ya(iya), za(iza),
		asi(1/(ixa*ixa+iya*iya+iza*iza)), rc(irc) {};
	bool point_inside(fpoint x,fpoint y,fpoint z);
	template<class n_option>
	inline void cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
	void cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
	void cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {cut_cell_base(c,x,y,z);}
};

#endif
