// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// Set up constants for the container geometry
const fpoint x_min=-15,x_max=15;
const fpoint y_min=-7,y_max=7;
const fpoint z_min=-15,z_max=15;

// Set the computational grid size
const int n_x=7,n_y=7,n_z=7;

// Set the number of particles that are going to be randomly introduced
const int particles=2000;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}


/** \brief A class representing a spherical wall object.
 *
 * This class represents a spherical wall object. */
struct wall_sphere_acc : public wall {
	public:
		/** Constructs a spherical wall object.
		 * \param[in] iw_id an ID number to associate with the wall for
		 * neighbor tracking. 
		 * \param[in] (ixc,iyc,izc) a position vector for the sphere's
		 * center.
		 * \param[in] irc the radius of the sphere. */
		wall_sphere_acc(fpoint ixc,fpoint iyc,fpoint izc,fpoint irc,int iw_id=-99)
			: w_id(iw_id), xc(ixc), yc(iyc), zc(izc), rc(irc) {};
		bool point_inside(fpoint x,fpoint y,fpoint z);
		template<class n_option>
		inline bool cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z);
		bool cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {return cut_cell_base(c,x,y,z);}
		bool cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {return cut_cell_base(c,x,y,z);}
		void min_distance(fpoint x,fpoint y,fpoint z,fpoint &dx,fpoint &dy,fpoint &dz);
	private:
		const int w_id;
		const fpoint xc,yc,zc,rc;
};

/** Tests to see whether a point is inside the sphere wall object.
 * \param[in] (x,y,z) the vector to test.
 * \return true if the point is inside, false if the point is outside. */ 
bool wall_sphere_acc::point_inside(fpoint x,fpoint y,fpoint z) {
	return (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)<rc*rc;
}

/** Cuts a cell by the sphere wall object. The spherical wall is approximated by
 * a single plane applied at the point on the sphere which is closest to the center
 * of the cell. This works well for particle arrangements that are packed against
 * the wall, but loses accuracy for sparse particle distributions.
 * \param[in] &c the Voronoi cell to be cut.
 * \param[in] (x,y,z) the location of the Voronoi cell.
 * \return true if the cell still exists, false if the cell is deleted. */
template<class n_option>
inline bool wall_sphere_acc::cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
	fpoint xd=x-xc,yd=y-yc,zd=z-zc,dq,idq;
	fpoint xd2,yd2,zd2,xd3,yd3,zd3,sinp;
	fpoint xs,ys,zs,sd;
	dq=xd*xd+yd*yd+zd*zd;
	if (dq>tolerance) {
		idq=1/sqrt(dq);xd*=idq;yd*=idq;zd*=idq;
		sinp=sqrt(1-zd*zd);
		if(sinp<tolerance) {dq=2*(sqrt(dq)*rc-dq);return c.nplane(xd,yd,zd,dq,w_id);}
		xd2=-xd*zd/sinp;yd2=-yd*zd/sinp;zd2=sinp;
		xd3=-yd/sinp;yd3=xd/sinp;zd3=0;
	//	cout << xd*xd+yd*yd+zd*zd << endl;
	//	cout << xd*xd2+yd*yd2+zd*zd2 << endl;
	//	cout << xd*xd3+yd*yd3+zd*zd3 << endl;
	//	cout << xd2*xd2+yd2*yd2+zd2*zd2 << endl;
	//	cout << xd2*xd3+yd2*yd3+zd2*zd3 << endl;
	//	cout << xd3*xd3+yd3*yd3+zd3*zd3 << endl << endl;
		fpoint i,j;
		for(i=-0.4;i<=0.4;i+=0.04) {
			for(j=-0.4;j<=0.4;j+=0.04) {
				xs=xd+i*xd2+j*xd3;
				ys=yd+i*yd2+j*yd3;
				zs=zd+i*zd2+j*zd3;
				sd=1/sqrt(xs*xs+ys*ys+zs*zs);
				xs*=sd;ys*=sd;zs*=sd;sd=rc-sqrt(dq/(1+i*i+j*j));
				if(!c.nplane(xs,ys,zs,2*sd,w_id)) return false;
			}
		}
	}
	return true;
}

/** Finds the vector of minimum distance from a given position vector to the
 * sphere wall object.
 * \param[in] (x,y,z) the vector to test.
 * \param[in] (&dx,&dy,&dz) the vector of minimum distance to the wall. */ 
void wall_sphere_acc::min_distance(fpoint x,fpoint y,fpoint z,fpoint &dx,fpoint &dy,fpoint &dz) {
	fpoint rad=sqrt(x*x+y*y+z*z);
	if(rad>tolerance) {
		rad=rc/rad-1;
		dx=x*rad;dy=y*rad;dz=z*rad;
	} else {
		dx=rc;dy=dz=0;
	}
}

int main() {
	int i=0,j,k;
	fpoint x,y,z,theta;//char q[256];
	int *u[201];

	for(j=0;j<=200;j++) {
		u[j]=new int[400];
		for(k=0;k<400;k++) u[j][k]=0;
	}

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Add a cylindrical wall to the container
	//wall_sphere_acc sph(0,0,0,4);
	//con.add_wall(sph);
	
	// Randomly add particles into the container
	while(i<particles) {
		theta=3.1415926535897932384626433832795*rnd();
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		if (x*x+y*y+z*z>4*4) continue;
		//x=x*x*0.2*(x<0?-1:1);
		//y=y*y*0.2*(y<0?-1:1);
		//z=z*z*0.2*(z<0?-1:1);
		con.put(i,x+10*sin(theta),y,z+10*cos(theta));i++;
	}

	for(i=0;i<=200;i++) {
		cout << i << " " << con.packing_badness<cond_all>() << endl;
		con.full_relax(0.8);
	}
	
	con.draw_particles_pov("acc_sphere_p.pov");
	con.draw_cells_pov("acc_sphere_v.pov");
	con.draw_particles("acc_pack");

	fpoint svol=4*3.1415926535897932384626433832795/3*4*4*4;
	fpoint vvol=con.sum_cell_volumes();
	cout << vvol << " " << svol << " " << (vvol/svol-1)*100 << "%" << endl;
}
