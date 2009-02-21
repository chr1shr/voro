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

class velocity_internal {
	public:
		velocity_internal(fpoint **&ive) : track_ve(true), ve(ive) {};
		inline void vel(int ijk,int q,fpoint &x,fpoint &y,fpoint &z) {
			x+=ve[ijk][3*q];
			y+=ve[ijk][3*q+1];
			z+=ve[ijk][3*q+2];
		}
		const bool track_ve;
	private:
		fpoint **&ve;
};

class velocity_brownian {
	public:
		velocity_brownian() : track_ve(false), mag(0.05), tmag(2*mag) {}; 
		inline void vel(int ijk,int q,fpoint &x,fpoint &y,fpoint &z) {
			x+=tmag*rnd()-mag;
			y+=tmag*rnd()-mag;
			z+=tmag*rnd()-mag;
		}
		const bool track_ve;
	private:
		const fpoint mag,tmag;
		inline fpoint rnd() {return fpoint(rand())/RAND_MAX;}
};

template<class r_option>
class container_dynamic_base : public container_base<r_option> {
	public:
		container_dynamic_base(fpoint xa,fpoint xb,fpoint ya,fpoint yb,fpoint za,fpoint zb,int xn,int yn,int zn,bool xper,bool yper,bool zper,int memi);
		~container_dynamic_base();
		using container_base<r_option>::xperiodic;
		using container_base<r_option>::yperiodic;
		using container_base<r_option>::zperiodic;
		using container_base<r_option>::ax;
		using container_base<r_option>::ay;
		using container_base<r_option>::az;
		using container_base<r_option>::nx;
		using container_base<r_option>::ny;
		using container_base<r_option>::nz;
		using container_base<r_option>::nxyz;
		using container_base<r_option>::xsp;
		using container_base<r_option>::ysp;
		using container_base<r_option>::zsp;
		using container_base<r_option>::co;
		using container_base<r_option>::p;
		using container_base<r_option>::sz;
		using container_base<r_option>::id;
		using container_base<r_option>::radius;
		using container_base<r_option>::wall_number;
		using container_base<r_option>::walls;
		using container_base<r_option>::mem;
		void wall_diagnostic();
		int count(fpoint x,fpoint y,fpoint z,fpoint r);
		void spot() {};
		void gauss_spot() {};
		void relax() {};
		void full_relax(fpoint alpha);
		template<class v_class>
		inline void move() {v_class vcl;move(&vcl);}
		template<class v_class>
		void move(v_class *vcl);
		void add_particle_memory(int i);
	protected:
		int *gh;
		fpoint **ve;
		velocity_internal v_inter;
	private:
		inline int step_mod(int a,int b);
		inline int step_int(fpoint a);
};

/** The basic dynamic container class. */
typedef container_dynamic_base<radius_mono> container_dynamic;

/** The polydisperse dynamic container class. */
typedef container_dynamic_base<radius_poly> container_dynamic_poly;

#endif
