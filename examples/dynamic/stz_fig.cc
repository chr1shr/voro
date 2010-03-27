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
	char q[256];
	sprintf(q,"%04d_p.pov",i);con.draw_particles_pov(q);
	sprintf(q,"%04d_v.pov",i);con.draw_cells_pov(q);
}

class cond_center {
	public:
		inline bool test(fpoint cx,fpoint cy,fpoint cz) {return cx>-15&&cx<15&&cz>15&&cz<75;}
};

class velocity_compress {
	public:
		velocity_compress() : track_ve(false), cfac(exp(0.001*(log(110)-log(100)))) {}; 
		inline void vel(int ijk,int q,fpoint &x,fpoint &y,fpoint &z) {
			x*=cfac;
			y/=cfac;
		}
		const bool track_ve;
	private:
		const fpoint cfac;
};

int main() {
	int i=0;
	fpoint x,y;
	const fpoint x1=-1.9,y1=0.6;
	const fpoint x2=2.5,y2=1.6;
	const fpoint x3=2.1,y3=-2.7;
	const fpoint x4=-1.9,y4=-3.0;
	const fpoint da=0.5,db=sqrt(3.0)/2;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(-4,4,-4,4,-0.5,0.5,4,4,1,
			true,true,false,8);

	con.put(i,x1+da,y1,0);i++;
	con.put(i,x1-da,y1,0);i++;
	con.put(i,x1,y1+db,0);i++;
	con.put(i,x1,y1-db,0);i++;
	con.put(i,x2+db,y2,0);i++;
	con.put(i,x2-db,y2,0);i++;
	con.put(i,x2,y2+da,0);i++;
	con.put(i,x2,y2-da,0);i++;
	con.put(i,x3+da,y3,0);i++;
	con.put(i,x3-da,y3,0);i++;
	con.put(i,x3,y3+db,0);i++;
	con.put(i,x3,y3-db,0);i++;
	con.put(i,x4+db,y4,0);i++;
	con.put(i,x4-db,y4,0);i++;
	con.put(i,x4,y4+da,0);i++;
	con.put(i,x4,y4-da,0);i++;

	// Randomly add particles into the container
	while(i<67) {
		x=4-8*rnd();
		y=4-8*rnd();
		con.put(i,x,y,0.1);
		i++;
	}

	output_all(con,0);
	for(i=0;i<1;i++) con.full_relax(2);
	for(i=0;i<5;i++) con.full_relax(1.2);
	for(i=0;i<100;i++) con.full_relax(0.8);

	output_all(con,1);

/*	con.neighbor_distribution<cond_center>(u[0],0.02,400);
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
	}*/
}
