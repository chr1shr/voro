// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

void output_all(container_dynamic &con,int i) {
/*	char q[256];
	sprintf(q,"output/%04d_p.pov",i);con.draw_particles_pov(q);
	sprintf(q,"gzip -f -9 output/%04d_p.pov",i);system(q);
	sprintf(q,"output/%04d_v.pov",i);con.draw_cells_pov(q);
	sprintf(q,"gzip -f -9 output/%04d_v.pov",i);system(q);*/
}

class cond_center {
	public:
		inline bool test(fpoint cx,fpoint cy,fpoint cz) {return cx>-15&&cx<15&&cz>15&&cz<75;}
};

class velocity_compress {
	public:
		velocity_compress() : track_ve(false), zfac(exp(0.001*(log(110)-log(110)))) {}; 
		inline void vel(int ijk,int q,fpoint &x,fpoint &y,fpoint &z) {
			z*=zfac;
		}
		const bool track_ve;
	private:
		const fpoint zfac;
};

int main() {
	int i=0,j,k;
	fpoint x,y,z;
	int *u[11];

	for(j=0;j<=10;j++) {
		u[j]=new int[400];
		for(k=0;k<400;k++) u[j][k]=0;
	}

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(-25,25,-4,4,0,115,25,4,75,
			false,false,false,8);
	
	// Randomly add particles into the container
	while(i<54708) {
		x=-24.5+rnd()*49;
		y=-3.5+rnd()*7;
		z=0.5+rnd()*109;
		con.put(i,x,y,z);
		i++;
	}

	output_all(con,0);
//	con.neighbor_distribution<cond_center>(u[0],0.02,400);
	con.clear_velocities();
	output_all(con,0);
	con.neighbor_distribution<cond_center>(u[0],0.02,400);
	for(i=1;i<=1000;i++) {
		con.move<velocity_compress>();
		if(i<40) con.full_relax(i*0.025+0.0);
		else if (i<800) con.full_relax(1);
		else con.full_relax(0.0+(1000-i)*0.0005);
	//	cout << i << " " << con.full_count() << " " << con.packing_badness<cond_all>() << endl;
		if(i%10==0) output_all(con,i/10);
		if(i%100==0) con.neighbor_distribution<cond_center>(u[i/100],0.02,400);
	}

	for(k=0;k<400;k++) {
		cout << (double(k)+0.5)*0.02;
		for(i=0;i<=10;i++) cout << " " << u[i][k];
		cout << "\n";
	}
}
