// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : July 1st 2008

#include "voro++.cc"
#include "dynamic.cc"

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

class cond_center {
	public:
		inline bool test(fpoint cx,fpoint cy,fpoint cz) {return cx>-15&&cx<15&&cz>15&&cz<75;}
};

const int trialz=8;

int main() {
	int i=0,j,k,c=1,parts,tr;
	fpoint x,y,z;
	fpoint alpha[200][50],al,best,trial;
	ofstream os;

	for(j=0;j<50;j++) alpha[0][j]=j<30?(j<7?0.5+fpoint(j)*0.1:1.2):1.2;

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
			k=rand()%9;
			if(k<2) {al*=1-0.2*rnd();} else if(k>6) {al*=1+0.2*rnd();}
			if(al>1.2) al=1.2;
			alpha[c][j]=al;
		}
	
		trial=0;
		parts=0;
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
			parts+=con.full_count();
		}

		cout << c << " " << trial/trialz << " " << double(parts)/trialz << endl;

		if(trial<best&&parts>=54708*trialz-2) {
			best=trial;c++;
			os.open("pat2",ofstream::out|ofstream::trunc);
			for(j=0;j<50;j++) {
				os << j;
				for(k=0;k<c;k++) os << " " << alpha[k][j];
				os << "\n";
			}
			os.close();
		}	
	}
}
