#ifndef VOROPP_PRE_CONTAINER_HH
#define VOROPP_PRE_CONTAINER_HH

#include "container.hh"

class pre_container_base {
	public:
		const double ax;
		const double bx;
		const double ay;
		const double by;
		const double az;
		const double bz;
		const bool xperiodic;
		const bool yperiodic;
		const bool zperiodic;
		pre_container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				bool xperiodic_,bool yperiodic_,bool zperiodic_,int ps_) :
			ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
			xperiodic(xperiodic_), yperiodic(yperiodic_), zperiodic(zperiodic_), ps(ps_)	{}
	protected:
		const int ps;
		/** This array holds pointers to any wall objects that have
		 * been added to the container. */
		wall **walls;
		/** The current number of wall objects, initially set to zero. */
		int wall_number;
		/** The current amount of memory allocated for walls. */
		int current_wall_size;
}

class pre_container {
	public:
		pre_container(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				bool xperiodic_,bool yperiodic_,bool zperiodic_)
			: pre_container_base(ax_,bx_,ay_,by_,az_,bz_,xperiodic_,yperiodic,zperiodic_) {};
	private:
}

class pre_container_poly {

}

#endif
