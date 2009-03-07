// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

//#include "voro++.cc"
//#include "dynamic.cc"
#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
#include <gsl/gsl_linalg.h>

// Set up constants for the container geometry
/*const fpoint x_min=-5,x_max=5;
const fpoint y_min=-5,y_max=5;
const fpoint z_min=-0.5,z_max=0.5;

// Set the computational grid size
const int n_x=8,n_y=8,n_z=1;

// Set the number of particles that are going to be randomly introduced
const int particles=100;

// This function returns a random double between 0 and 1
*/

const double pi=3.1415926535897932384626433832795;
inline double arg(double x,double y) { 
return x+y>0?(x>y?atan(y/x):pi*0.5-atan(x/y)):(x>y?-atan(x/y)-pi*0.5:atan(y/x)+(y>0?pi:-pi)); 
}

double rnd() {return double(rand())/RAND_MAX;}

/*void output_all(container_dynamic &con,int i) {
	char q[256];
	sprintf(q,"output/%04d_p.pov",i);con.draw_particles_pov(q);
	sprintf(q,"gzip -f -9 output/%04d_p.pov",i);system(q);
	sprintf(q,"output/%04d_v.pov",i);con.draw_cells_pov(q);
	sprintf(q,"gzip -f -9 output/%04d_v.pov",i);system(q);
}*/

class wall_helix {
	public:
		wall_helix(double ir,double il) : r(ir), l(il),
			lilr(il/sqrt(ir*ir+il*il)), rilr(ir/sqrt(ir*ir+il*il)),
			og(gsl_matrix_view_array(o,3,3)),
			vg(gsl_vector_view_array(v,3)),
       			eg(gsl_vector_view_array(e,3)),
			pg(gsl_permutation_alloc(3)) {};
		void convert(double x,double y,double z,double &a,double &t,double &p) {
			int s,i=0,j=0;
			double cp,ct,sp,st,racp;
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

			/*	printf ("\no = \n");
  				gsl_matrix_fprintf(stdout,&og.matrix,"%g");
				printf ("\nv = \n");
  				gsl_vector_fprintf(stdout,&vg.vector,"%g");*/
				if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]<1e-20) j++;
				i++;if(i>100) {cout << "Newton-Raphson didn't converge" << endl;}

				gsl_linalg_LU_decomp(&og.matrix,pg,&s);
				gsl_linalg_LU_solve(&og.matrix,pg,&vg.vector,&eg.vector);
			//	cout << v[0]*v[0]+v[1]*v[1]+v[2]*v[2] << " " << a << " " << p << " " << t << " " << e[0] << " " << e[1] << " " << e[2] << endl;
				a+=e[0];t+=e[1];p+=e[2];
			} while (j<5);
		}
		void compute(double a,double t,double p,double &x,double &y,double &z) {
			double racp=r+a*cos(p);
			x=racp*cos(t)+a*lilr*sin(t)*sin(p);
			y=racp*sin(t)-a*lilr*cos(t)*sin(p);
			z=t*l+a*rilr*sin(p);
		}
	private:
		const double r,l,lilr,rilr;
		double o[9],v[3],e[3];
		gsl_matrix_view og;
		gsl_vector_view vg;
		gsl_vector_view eg;
		gsl_permutation *pg;
};

int main() {
	const double lam=0.5;
	int i=0;
	double x,y,z,dx,dy,dz,a,p,t,rad,temp1,temp2;
	wall_helix wh(4,lam);
	
	while(i<5000) {
		x=5*(2*rnd()-1);
		y=5*(2*rnd()-1);
		z=5*(2*rnd()-1);
		rad=x*x+y*y;
		if(rad>9&&rad<25) {
			t=arg(x,y);
			if(t<z/lam-pi) {
				do {
					t+=2*pi;
				} while(t<z/lam-pi);
			} else {
				while(t>z/lam+pi) t-=2*pi;
			}
			temp1=z-t*lam;
			temp2=sqrt(rad)-4;
			a=sqrt(temp1*temp1+temp2*temp2);
			p=arg(temp2,temp1);
		//	cout << x << " " << y << " " << z << " " << a << " " << t << " " << p << endl;
		//	cout << x << " " << y << " " << z << " " << cos(t)*(4+a*cos(p)) << " " << sin(t)*(4+a*cos(p)) << " " << t*lam+a*sin(p) << endl;
			wh.convert(x,y,z,a,t,p);
			if(a<1&&a>0) {
				a=1;
				wh.compute(a,t,p,dx,dy,dz);
				cout << x << " " << y << " " << z << " " << dx-x << " " << dy-y << " " << dz-z << endl;
			}
			i++;
		}
	}
}



/*
output_all(con,0);
for(i=0;i<=10;i++) con.full_relax(0.8);
	output_all(con,1);

	for(i=20;i<=1000;i++) {
		con.spot(0,0,0,0.005,0.005,0,2.5);
	//	con.relax(4,4,0,3.5,0.8);
		con.full_relax(0.8);
		if(i%10==0) {
			output_all(con,i/10);
		}
	}
}*/
