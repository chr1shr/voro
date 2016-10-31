// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include <cstdio>

#include "config.hh"
#include "voro++.cc"
#include "dynamic.cc"

#include "Magick++.h"

using namespace Magick;

// Set up constants for the container geometry
const fpoint x_min=0,x_max=20;
const fpoint y_min=0,y_max=20;
const fpoint z_min=-0.5,z_max=0.5;

const double pi=3.1415926535897932384626433832795;

// Set the computational grid size
const int n_x=8,n_y=8,n_z=1;

// Set the number of particles that are going to be randomly introduced
const int particles=1440;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

void output_all(container_dynamic_poly &con,int i) {
	char q[256];
	sprintf(q,"output/%04d_p.pov",i);con.draw_particles(q);
}

inline double rfunc(double x) {
	return x*x*(3-2*x);
}

inline double arg(double x,double y) {
	return x+y>0?(x>y?atan(y/x):pi*0.5-atan(x/y)):
                     (x>y?-atan(x/y)-pi*0.5:atan(y/x)+(y>0?pi:-pi));
}

void rainbow2(double arg,double &re,double &gr,double &bl) {
	re=0.5+0.5*cos(arg);
	gr=0.5+0.5*cos(arg-2.09439510239319549230842892219);
	bl=0.5+0.5*cos(arg-4.18879020478639098461685784436);
	re=rfunc(re);
	gr=rfunc(gr);
	bl=rfunc(bl);
	re+=0.1*gr;if(re>1) re=1;
	bl+=0.1*gr;if(bl>1) bl=1;
	gr*=0.7;
//	re=0.3+0.7*re;
//	gr=0.3+0.7*gr;
//	bl=0.3+0.7*bl;
}

int main() {
	int i=0,j;
	fpoint x,y,z,r;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic_poly con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			true,true,false,8);

	// Add a cylindrical wall to the container
	wall_sphere sph(10,10,0,6);
	con.add_wall(sph);
	wall_sphere sph2(10,10,0,9);
	con.add_wall(sph2);
	
	// Randomly add particles into the container
	while(i<particles) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		r=0.1+0.6*rnd()*rnd();
		con.put(i,x,y,0,r);
		i++;
	}

//	output_all(con,0);
	for(i=0;i<=4000;i++) {
		con.move<velocity_brownian>();
		con.full_relax(1.6);
	}
//	output_all(con,1);
	for(i=0;i<=1000;i++) con.full_relax(1.1);
	for(i=0;i<=8001;i++) con.full_relax(0.5);
	for(i=0;i<=8001;i++) con.full_relax(0.2);
//	output_all(con,2);

	double rr,re,gr,bl,scf=1200/20.0;
	double **p=con.p;int *co=con.co;
//	Image im(Geometry(1200,1200),Color("white"));
//	im.strokeWidth(0);
	for(i=0;i<n_x*n_y*n_z;i++) {
		for(j=0;j<co[i];j++) {
			x=p[i][4*j];y=p[i][4*j+1];r=p[i][4*j+3];
			rr=sqrt((x-10)*(x-10)+(y-10)*(y-10));
			if(rr<6||rr>9) continue;
			x*=scf;y*=scf;r*=scf;
			z=arg(x-600,y-600)+rnd()-0.5;
		//	z=tpi/8*(int(8.0/tpi*(tpi+z)));
		//	rainbow2(tpi*rnd(),re,gr,bl);
			rainbow2(z,re,gr,bl);

			printf("\\put(%d,%d){\\color[rgb]{%.4f,%.4f,%.4f}\\circle*{%d}\n",int(x),int(y),re,gr,bl,int(r*2));
	//		im.fillColor(Color(re*MaxRGB,gr*MaxRGB,bl*MaxRGB,0));
	//		im.draw(DrawableCircle(x,y,x+0.98*r,y));
		}
	}
//	im.write("png24:test.png");
}
