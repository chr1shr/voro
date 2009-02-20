// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

/** \file dynamic.hh
 * \brief Header file for the dynamic extension classes, which add
 * functionality for a variety of dynamic particle motions. */

#ifndef VOROPP_DYNAMIC_HH
#define VOROPP_DYNAMIC_HH

template<class r_option>
class container_dynamic_base : public container_base<r_option> {
	public:
		container_dynamic_base(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi) : container_base<r_option> (xa,xb,ya,yb,za,zb,xn,yn,zn,xper,yper,zper,memi);
		using container_base<r_option>::nx;
		using container_base<r_option>::ny;
		using container_base<r_option>::nz;
		using container_base<r_option>::co;
		using container_base<r_option>::p;
		using container_base<r_option>::sz;
		using container_base<r_option>::id;
		using container_base<r_option>::radius;
		using container_base<r_option>::wall_number;
		using container_base<r_option>::walls;
		void wall_diagnostic();
		void count() {};
		void spot() {};
		void gauss_spot() {};
		void relax() {};
		void full_relax() {};
	protected:

};

/** The basic dynamic container class. */
typedef container_dynamic_base<radius_mono> container_dynamic;

/** The polydisperse dynamic container class. */
typedef container_dynamic_base<radius_poly> container_dynamic_poly;

#endif
