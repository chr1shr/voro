// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}
const double pi=3.1415926535897932384626433832795;



const int trialz=250;

int main() {
	int i=0,j,k,c=1,parts,tr;
	fpoint x,y,z;bool die=false;
	fpoint alpha[200][50],al,best,trial;
	ofstream os;

	for(j=0;j<50;j++) alpha[0][j]=0.2;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container_dynamic con(-25,25,-4,4,0,112,25,4,75,
			false,false,false,8);
	
	// Randomly add particles into the container
	best=0;parts=0;
	for(tr=0;tr<trialz;tr++) {
		
		con.clear();
		for(i=0;i<54708;i++) {
			x=-24.5+rnd()*49;
			y=-3.5+rnd()*7;
			z=0.5+rnd()*111;
			con.put(i,x,y,z);
		}
		for(j=0;j<50;j++) con.full_relax(alpha[0][j]);
		best+=con.packing_badness<cond_all>();
		parts+=con.full_count();
	}

	cout << "0 " << best/trialz << " " << double(parts)/trialz << endl;

	while(c<200) {
		for(j=0;j<50;j++) {
			al=alpha[c-1][j];
	//		k=rand()%9;
	//		if(k<2) {al*=1-0.02*rnd();} else if(k>6) {al*=1+0.02*rnd();}
			alpha[c][j]=al;
		}
		for(k=0;k<16;k++) {
			al=0.12*(1-0.05*k)*(2*rnd()-1);
			for(j=0;j<50;j++) alpha[c][j]*=(1+al*cos(pi*k*j/50.0));
		}
	
		trial=0;
		for(tr=0;tr<trialz;tr++) {
			con.clear();
			for(i=0;i<54708;i++) {
				x=-24.5+rnd()*49;
				y=-3.5+rnd()*7;
				z=0.5+rnd()*111;
				con.put(i,x,y,z);
			}
			for(j=0;j<50;j++) con.full_relax(alpha[c][j]);
			trial+=con.packing_badness<cond_all>();
			if(con.full_count()!=54708) {die=true;break;}
		}
		if(die) {die=false;continue;}

		cout << c << " " << trial/trialz << endl;

		if(trial<best) {
			best=trial;c++;
			os.open("pat",ofstream::out|ofstream::trunc);
			for(j=0;j<50;j++) {
				os << j;
				for(k=0;k<c;k++) os << " " << alpha[k][j];
				os << "\n";
			}
			os.close();
		}	
	}
}
