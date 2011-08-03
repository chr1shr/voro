#ifndef VOROPP_UNITCELL_HH
#define VOROPP_UNITCELL_HH

#include <vector>
using namespace std;

#include "config.hh"
#include "cell.hh"

class unitcell {
	public:
		const double bx;
		const double bxy,by;
		const double bxz,byz,bz;
		voronoicell unit_voro;
		unitcell(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_);
		/** Draws an outline of the domain in Gnuplot format.
		 * \param[in] filename the filename to write to. */
		inline void draw_domain_gnuplot(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_domain_gnuplot(fp);
			fclose(fp);
		}
		void draw_domain_gnuplot(FILE *fp=stdout);
		/** Draws an outline of the domain in Gnuplot format.
		 * \param[in] filename the filename to write to. */
		inline void draw_domain_pov(const char* filename) {
			FILE *fp(voropp_safe_fopen(filename,"w"));
			draw_domain_pov(fp);
			fclose(fp);
		}
		void draw_domain_pov(FILE *fp=stdout);
		bool intersects_image(double dx,double dy,double dz,double &vol);
		void images(vector<int> &vi,vector<double> &vd);
	protected:
		double max_uv_y;
		double max_uv_z;
	private:
		inline void unit_voro_apply(int i,int j,int k);
		bool unit_voro_intersect(int l);
		inline bool unit_voro_test(int i,int j,int k);
};

#endif
