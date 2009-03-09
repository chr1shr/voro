// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"
#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
#include <gsl/gsl_linalg.h>

const fpoint pi=3.1415926535897932384626433832795;

// Set up constants for the container geometry
const fpoint x_min=-16,x_max=16;
const fpoint y_min=-16,y_max=16;
const fpoint z_min=-40,z_max=32;

// Set the computational grid size
const int n_x=16,n_y=16,n_z=32;

// Helix parameters
const fpoint aradius=7.5,bradius=4,lam=3,upz=20,downz=-32,insbuf=0.4;

// Set the number of particles that are going to be randomly introduced
const int particles=pi*bradius*bradius*aradius*(upz-downz)/lam*1.2;

// This function returns a random fpoint between 0 and 1
fpoint rnd() {return fpoint(rand())/RAND_MAX;}

inline fpoint arg(fpoint x,fpoint y) {
	return x+y>0?(x>y?atan(y/x):pi*0.5-atan(x/y)):(x>y?-atan(x/y)-pi*0.5:atan(y/x)+(y>0?pi:-pi)); 
}

class wall_helix_end;

class wall_helix : public wall {
	public:
		wall_helix(fpoint ir,fpoint il,fpoint ib,int iw_id=-99) : w_id(iw_id), b(ib),
			r(ir), l(il), lilr(il/sqrt(ir*ir+il*il)), rilr(ir/sqrt(ir*ir+il*il)),
			og(gsl_matrix_view_array(o,3,3)),
			vg(gsl_vector_view_array(v,3)),
       			eg(gsl_vector_view_array(e,3)),
			pg(gsl_permutation_alloc(3)) {};
		bool point_inside(fpoint x,fpoint y,fpoint z) {
			const fpoint rmin=b>r?0:(b-r)*(b-r),rmax=(b+r)*(b+r);
			fpoint rad=x*x+y*y,a,t,p;
			if(rad<rmin||rad>rmax) return false;
			convert(x,y,z,a,t,p);
			return a<b;
		}
		void convert(fpoint x,fpoint y,fpoint z,fpoint &a,fpoint &t,fpoint &p) {
			fpoint temp1,temp2,rad=x*x+y*y;
			t=arg(x,y);
			if(t<z/l-pi) {
				do {
					t+=2*pi;
				} while(t<z/l-pi);
			} else {
				while(t>z/l+pi) t-=2*pi;
			}
			temp1=z-t*lam;
			temp2=sqrt(rad)-4;
			a=sqrt(temp1*temp1+temp2*temp2);
			p=arg(temp2,temp1);
			int s,i=0,j=0;
			fpoint cp,ct,sp,st,racp;
			do {
				cp=cos(p);ct=cos(t);
				sp=sin(p);st=sin(t);
				racp=r+a*cp;

				o[0]=cp*ct+lilr*st*sp;
				o[1]=-racp*st+a*lilr*ct*sp;
				o[2]=-a*sp*ct+a*lilr*cp*st;

				o[3]=cp*st-lilr*ct*sp;
				o[4]=racp*ct+a*lilr*st*sp;
				o[5]=-a*sp*st-a*lilr*cp*ct;

				o[6]=rilr*sp;
				o[7]=l;
				o[8]=a*rilr*cp;

				v[0]=x-(racp*ct+a*lilr*st*sp);
				v[1]=y-(racp*st-a*lilr*ct*sp);
				v[2]=z-(t*l+a*rilr*sp);

				if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]<tolerance*tolerance) j++;
				i++;if(i>100) {cout << "Newton-Raphson didn't converge" << endl;exit(0);}
				
				gsl_linalg_LU_decomp(&og.matrix,pg,&s);
				gsl_linalg_LU_solve(&og.matrix,pg,&vg.vector,&eg.vector);
				a+=e[0];t+=e[1];p+=e[2];
			} while (j<5);
			if(a<0) {p+=pi;a=-a;}
			//cout << a << " " << t << " c " << p << endl;
		}
		template<class n_option>
		inline bool cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
			fpoint dx,dy,dz,a,t,p;
			convert(x,y,z,a,t,p);
			dx=(cos(p)*cos(t)+lilr*sin(t)*sin(p));
			dy=(cos(p)*sin(t)-lilr*cos(t)*sin(p));
			dz=rilr*sin(p);
			return c.plane(dx,dy,dz,2*(b-a));
		}
		void min_distance(fpoint x,fpoint y,fpoint z,fpoint &dx,fpoint &dy,fpoint &dz) {
			fpoint a,t,p;
			convert(x,y,z,a,t,p);
			dx=(b-a)*(cos(p)*cos(t)+lilr*sin(t)*sin(p));
			dy=(b-a)*(cos(p)*sin(t)-lilr*cos(t)*sin(p));
			dz=(b-a)*rilr*sin(p);
		}
		void compute(fpoint a,fpoint t,fpoint p,fpoint &x,fpoint &y,fpoint &z) {
			fpoint racp=r+a*cos(p);
			x=racp*cos(t)+a*lilr*sin(t)*sin(p);
			y=racp*sin(t)-a*lilr*cos(t)*sin(p);
			z=t*l+a*rilr*sin(p);
		}
		bool cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {return cut_cell_base(c,x,y,z);}
		bool cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		const fpoint b,r,l,lilr,rilr;
		fpoint o[9],v[3],e[3];
		gsl_matrix_view og;
		gsl_vector_view vg;
		gsl_vector_view eg;
		gsl_permutation *pg;
		friend class wall_helix_end;
};

class wall_helix_end : public wall {
	public:
		wall_helix_end(wall_helix &iwh,fpoint z,fpoint upsign) :
			w_id(iwh.w_id), wh(iwh), radc((iwh.b+0.5)*(iwh.b+0.5)), walld(2) {
			fpoint t=z/wh.l;
			xc=-wh.rilr*sin(t)*upsign;yc=wh.rilr*cos(t)*upsign;zc=wh.lilr*upsign;  // Normal
			xa=wh.r*cos(t);ya=wh.r*sin(t);za=t*wh.l; // Pos vector of axis
		}
		template<class n_option>
		inline bool cut_cell_base(voronoicell_base<n_option> &c,fpoint x,fpoint y,fpoint z) {
			fpoint xp=x-xa,yp=y-ya,zp=z-za;
			fpoint dp=xp*xc+yp*yc+zp*zc;
			xp-=dp*xc;yp-=dp*yc;zp-=dp*zc;
//			return c.plane(xc,yc,zc,-2*dp);
			return xp*xp+yp*yp+zp*zp<radc&&dp>-walld?c.plane(xc,yc,zc,-2*dp):true;
		}
		bool point_inside(fpoint x,fpoint y,fpoint z) {
			fpoint xp=x-xa,yp=y-ya,zp=z-za;
			fpoint dp=xp*xc+yp*yc+zp*zc;
			xp-=dp*xc;yp-=dp*yc;zp-=dp*zc;
			return (xp*xp+yp*yp+zp*zp>radc||dp<0);
		}
		void min_distance(fpoint x,fpoint y,fpoint z,fpoint &dx,fpoint &dy,fpoint &dz) {
			fpoint xp=x-xa,yp=y-ya,zp=z-za;
			fpoint dp=xp*xc+yp*yc+zp*zc;
			xp-=dp*xc;yp-=dp*yc;zp-=dp*zc;
			if(xp*xp+yp*yp+zp*zp<radc&&dp>-walld) {
				dx=-xc*dp;dy=-yc*dp;dz=-zc*dp;
			} else {
				dx=1e4;dy=1e4;dz=1e4;
			}
		}
		bool cut_cell(voronoicell_base<neighbor_none> &c,fpoint x,fpoint y,fpoint z) {return cut_cell_base(c,x,y,z);}
		bool cut_cell(voronoicell_base<neighbor_track> &c,fpoint x,fpoint y,fpoint z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		wall_helix &wh;
		const fpoint radc,walld;
		fpoint xc,yc,zc;
		fpoint xa,ya,za;
		fpoint r;
};

int main() {
	int i=0;
	fpoint a,t,p,x,y,z;

	container_dynamic con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);

	wall_helix wh(aradius,lam,bradius);
	wall_helix_end whe1(wh,upz,1);
	wall_helix_end whe2(wh,downz,-1);
	con.add_wall(wh);
	con.add_wall(whe1);
	con.add_wall(whe2);
	
	while(i<particles) {
		a=(bradius-insbuf)*rnd();
		t=(downz/lam+insbuf/sqrt(lam*lam+aradius*aradius))+rnd()*((upz-downz)/lam-2*insbuf/sqrt(lam*lam+aradius*aradius));
		p=2*pi*rnd();
		wh.compute(a,t,p,x,y,z);
		//cout << a << " " << t << " " << p << endl;
		if(!con.point_inside(x,y,z)) {cout << i << " " << x << " " << y << " " << z << endl;return 0;}
		con.put(i,x,y,z);i++;
	}
	for(i=0;i<200;i++) {
		con.full_relax(0.05*(i<12?i:12));
		cout << i << " " << con.packing_badness<cond_all>() << endl;
	}

	con.draw_particles_pov("helix_p.pov");
	con.draw_cells_pov("helix_v.pov");
}
