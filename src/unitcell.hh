#ifndef VOROPP_UNITCELL_HH
#define VOROPP_UNITCELL_HH

#include "config.hh"
#include "cell.hh"

class unitcell {
	public:
		const double bx;
		const double bxy,by;
		const double bxz,byz,bz;
		voronoicell unit_voro;
		unitcell(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_);
		void draw_domain_gnuplot(FILE *fp=stdout);
		void draw_domain_pov(FILE *fp=stdout);
		bool intersects_image(double dx,double dy,double dz,double &vol);
	private:
		inline void unit_voro_apply(int i,int j,int k);
		bool unit_voro_intersect(int l);
		inline bool unit_voro_test(int i,int j,int k);
};

#endif
