// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main(int argc,char **argv) {
	if(argc!=8&&argc!=9) {
		cerr << "Syntax: ./gen_packing <x_size> <y_size> <z_size> { <num> | <frac>% }\n"
			"                      <num_relax> <relax_param> <filename> [<seed>]\n";
		return 1;
	}
	int i;

	// Read the box dimensions and check they make sense
	fpoint x(atof(argv[1])),y(atof(argv[2])),z(atof(argv[3])),px,py,pz;
	if(x<=0||y<=0||z<=0) {
		cerr << "The box dimensions must be positive\n";
		return 1;
	}

	// Read particle number, or compute it if a percent sign is detected
	int n=strlen(argv[4]);double packf;
	if(argv[4][n-1]=='%') {
		char *tmp=new char[n];
		for(i=0;i<n-1;i++) tmp[i]=argv[4][i];tmp[n-1]=0;
		packf=atof(argv[4]);	
		n=int (x*y*z*0.0190985931710274*packf+0.5);
		delete [] tmp;
	} else n=atoi(argv[4]);
	packf=double(n)/(0.0190985931710274*x*y*z);

	// Read in the relaxation number
	int nr(atoi(argv[5]));
	if(n<=0) {
		cerr << "There must be a strictly positive number of particles\n";
		return 1;
	}
	if(nr<0) {
		cerr << "There must be a positive number of particles\n";
		return 1;
	}
	
	// Check relaxation parameter
	double rp(atof(argv[6]));
	if(rp<=0||rp>=1)
		cerr << "Warning: relaxation value of " << rp << " is outside range of (0,1)\n";
	
	// Seed random number generator if requested
	if(argc==9) srand(atoi(argv[8]));

	// Compute optimal size of internal computational grid
	fpoint ilscale(pow(n/(5.0*x*y*z),1/3.0));
	int nx=int(x*ilscale+1),ny=int(y*ilscale+1),nz=int(z*ilscale+1);


	// Output diagnostic information 
	cout << "Using " << n << " particles, " << packf << "% packing fraction\n";
	cout << "Computational grid : " << nx << " by " << ny << " by " << nz << "\n\n";
	
	// Create particle container
	container_dynamic con(0,x,0,y,0,z,nx,ny,nz,false,false,false,8);
	
	// Randomly add particles into the container
	for(i=0;i<n;i++) {
		px=rnd()*x;
		py=rnd()*y;
		pz=rnd()*z;
		con.put(i,px,py,pz);
	}

	// Carry out relaxations
	cout << "Packing badnesses:\nInitial : " << con.packing_badness<cond_all>() << "\n";
	for(i=1;i<=nr;i++) {
		con.full_relax(rp);
		if(i%10==0||i==nr) cout << "Relaxation " << i << " : " << con.packing_badness<cond_all>() << "\n";
	}

	// Output particle configuration
	con.draw_particles(argv[7]);
}
