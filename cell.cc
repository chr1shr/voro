// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"

// Constructs a Voronoi cell and sets up the initial memory
voronoicell::voronoicell() :
	currentvertices(initvertices), currentvertexorder(initvertexorder),
	currentdeletesize(initdeletesize), currentdeletesize2(initdeletesize2) {
	int i;
	ds=new int[currentdeletesize];
	ds2=new int[currentdeletesize2];
	mem=new int[currentvertexorder];
	mec=new int[currentvertexorder];
	mep=new int*[currentvertexorder];
	ed=new int*[currentvertices];
#ifdef FACETS_NEIGHBOR
	mne=new int*[currentvertexorder];
	ne=new int*[currentvertices];
#endif
	nu=new int[currentvertices];
	pts=new f_point[3*currentvertices];
	sure.p=pts;
	for(i=0;i<3;i++) {
		mem[i]=initnvertices;
		mep[i]=new int[initnvertices*(2*i+1)];
		mec[i]=0;
#ifdef FACETS_NEIGHBOR
		mne[i]=new int[initnvertices*i];
#endif
	}
	mem[3]=init3vertices;
	mep[3]=new int[init3vertices*7];
	mec[3]=0;
#ifdef FACETS_NEIGHBOR
	mne[3]=new int[init3vertices*3];
#endif
	for(i=4;i<currentvertexorder;i++) {
		mem[i]=initnvertices;
		mep[i]=new int[initnvertices*(2*i+1)];
		mec[i]=0;
#ifdef FACETS_NEIGHBOR
		mne[i]=new int[initnvertices*i];
#endif
	}
};

voronoicell::~voronoicell() {
	delete [] ds;
	delete [] ds2;
	for(int i=0;i<currentvertexorder;i++) if (mem[i]>0) {
		delete [] mep[i];
#ifdef FACETS_NEIGHBOR
		delete [] mne[i];
#endif
	}
	delete [] mem;
	delete [] mec;
	delete [] mep;
	delete [] ed;
	delete [] nu;
	delete [] pts;
};

// Increases the memory storage for a particular vertex order
void voronoicell::addmemory(int i) {
	int s=2*i+1;
	if(mem[i]==0) {
#ifdef FACETS_NEIGHBOR
		mne[i]=new int[initnvertices*i];
#endif		
		mep[i]=new int[initnvertices*s];
		mem[i]=initnvertices;
		cerr << "Order " << i << " vertex memory created " << endl;
	} else {
		int j,k,*l;
		mem[i]*=2;
		if (mem[i]>maxnvertices) throw fatal_error("Point memory allocation exceeded absolute maximum");
		cerr << "Order " << i << " vertex memory scaled up to " << mem[i] << endl;
		l=new int[s*mem[i]];
		j=0;
		while(j<s*mec[i]) {
			k=mep[i][j+2*i];
			if(k>=0) {
				ed[k]=l+j;
			} else {
				int m;
				for(m=0;m<stack2;m++) {
					if(ed[ds2[m]]==mep[i]+j) {
						ed[ds2[m]]=l+j;
						break;
					}
				}
				if(m==stack2) throw fatal_error("Couldn't relocate dangling pointer");
				cerr << "Relocated dangling pointer" << endl;
			}
			for(k=0;k<=2*i;k++) {
				l[j]=mep[i][j];
				j++;
			}
		}
		delete [] mep[i];
		mep[i]=l;
#ifdef FACETS_NEIGHBOR
		int *q;
		q=new int[i*mem[i]];
		for(j=0;j<i*mec[i];j++) q[j]=mne[i][j];
		delete [] mne[i];
		mne[i]=q;
#endif
	}
};

void voronoicell::addmemory_vertices() {
	int i=2*currentvertices,j,**ped,*pnu;
	if (i>maxvertices) throw fatal_error("Vertex memory allocation exceeded absolute maximum");
	cerr << "Vertex memory scaled up to " << i << endl;
	f_point *ppts;
	ped=new int*[i];
	for(j=0;j<currentvertices;j++) ped[j]=ed[j];
	delete [] ed;ed=ped;
	pnu=new int[i];
	for(j=0;j<currentvertices;j++) pnu[j]=nu[j];
	delete [] nu;nu=pnu;
	ppts=new f_point[3*i];
	for(j=0;j<3*currentvertices;j++) ppts[j]=pts[j];
	delete [] pts;sure.p=pts=ppts;
	currentvertices=i;
};

void voronoicell::addmemory_vorder() {
	int i=2*currentvertexorder,j,*p1,**p2;
	if (i>maxvertexorder) throw fatal_error("Vertex order memory allocation exceeded absolute maximum");
	cerr << "Vertex order memory scaled up to " << i << endl;
	p1=new int[i];
	for(j=0;j<currentvertexorder;j++) p1[j]=mem[j];while(j<i) p1[j++]=0;
	delete [] mem;mem=p1;
	p2=new int*[i];
	for(j=0;j<currentvertexorder;j++) p2[j]=mep[j];
	delete [] mep;mep=p2;
	p1=new int[i];
	for(j=0;j<currentvertexorder;j++) p1[j]=mec[j];while(j<i) p1[j++]=0;
	delete [] mec;mec=p1;
#ifdef FACETS_RADIUS
	p2=new int*[i];
	for(j=0;j<currentvertexorder;j++) p2[j]=mne[j];
	delete [] mne;mne=p2;
#endif
	currentvertexorder=i;
};

void voronoicell::addmemory_ds() {
	int i=2*currentdeletesize,j,*pds;
	if (i>maxdeletesize) throw fatal_error("Delete stack 1 memory allocation exceeded absolute maximum");
	cerr << "Delete stack 1 memory scaled up to " << i << endl;
	pds=new int[i];
	for(j=0;j<currentdeletesize;j++) pds[j]=ds[j];
	delete [] ds;ds=pds;
	currentdeletesize=i;
};

void voronoicell::addmemory_ds2() {
	int i=2*currentdeletesize2,j,*pds2;
	if (i>maxdeletesize2) throw fatal_error("Delete stack 2 memory allocation exceeded absolute maximum");
	cerr << "Delete stack 2 memory scaled up to " << i << endl;
	pds2=new int[i];
	for(j=0;j<currentdeletesize2;j++) pds2[j]=ds2[j];
	delete [] ds2;ds2=pds2;
	currentdeletesize2=i;
};

// Initializes a Voronoi cell as a rectangular box with the given dimensions
inline void voronoicell::init(f_point xmin,f_point xmax,f_point ymin,f_point ymax,f_point zmin,f_point zmax) {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;
	mec[3]=p=8;xmin*=2;xmax*=2;ymin*=2;ymax*=2;zmin*=2;zmax*=2;
	pts[0]=xmin;pts[1]=ymin;pts[2]=zmin;
	pts[3]=xmax;pts[4]=ymin;pts[5]=zmin;
	pts[6]=xmin;pts[7]=ymax;pts[8]=zmin;
	pts[9]=xmax;pts[10]=ymax;pts[11]=zmin;
	pts[12]=xmin;pts[13]=ymin;pts[14]=zmax;
	pts[15]=xmax;pts[16]=ymin;pts[17]=zmax;
	pts[18]=xmin;pts[19]=ymax;pts[20]=zmax;
	pts[21]=xmax;pts[22]=ymax;pts[23]=zmax;
	int *q=mep[3];
	q[0]=1;q[1]=4;q[2]=2;q[3]=2;q[4]=1;q[5]=0;q[6]=0;
	q[7]=3;q[8]=5;q[9]=0;q[10]=2;q[11]=1;q[12]=0;q[13]=1;
	q[14]=0;q[15]=6;q[16]=3;q[17]=2;q[18]=1;q[19]=0;q[20]=2;
	q[21]=2;q[22]=7;q[23]=1;q[24]=2;q[25]=1;q[26]=0;q[27]=3;
	q[28]=6;q[29]=0;q[30]=5;q[31]=2;q[32]=1;q[33]=0;q[34]=4;
	q[35]=4;q[36]=1;q[37]=7;q[38]=2;q[39]=1;q[40]=0;q[41]=5;
	q[42]=7;q[43]=2;q[44]=4;q[45]=2;q[46]=1;q[47]=0;q[48]=6;
	q[49]=5;q[50]=3;q[51]=6;q[52]=2;q[53]=1;q[54]=0;q[55]=7;
	ed[0]=q;ed[1]=q+7;ed[2]=q+14;ed[3]=q+21;
	ed[4]=q+28;ed[5]=q+35;ed[6]=q+42;ed[7]=q+49;
#ifdef FACETS_NEIGHBOR
	q=mne[3];
	q[0]=-5;q[1]=-3;q[2]=-1;
	q[3]=-5;q[4]=-2;q[5]=-3;
	q[6]=-5;q[7]=-1;q[8]=-4;
	q[9]=-5;q[10]=-4;q[11]=-2;
	q[12]=-6;q[13]=-1;q[14]=-3;
	q[15]=-6;q[16]=-3;q[17]=-2;
	q[18]=-6;q[19]=-4;q[20]=-1;
	q[21]=-6;q[22]=-2;q[23]=-4;
	ne[0]=q;ne[1]=q+3;ne[2]=q+6;ne[3]=q+9;
	ne[4]=q+12;ne[5]=q+15;ne[6]=q+18;ne[7]=q+21;
#endif
	nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=nu[6]=nu[7]=3;
};

// Initializes a Voroni cell as a regular octahedron
inline void voronoicell::init_octahedron(f_point l) {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;
	mec[4]=p=6;l*=2;
	pts[0]=-l;pts[1]=0;pts[2]=0;
	pts[3]=l;pts[4]=0;pts[5]=0;
	pts[6]=0;pts[7]=-l;pts[8]=0;
	pts[9]=0;pts[10]=l;pts[11]=0;
	pts[12]=0;pts[13]=0;pts[14]=-l;
	pts[15]=0;pts[16]=0;pts[17]=l;
	int *q=mep[4];
	q[0]=2;q[1]=5;q[2]=3;q[3]=4;q[4]=0;q[5]=0;q[6]=0;q[7]=0;q[8]=0;
	q[9]=2;q[10]=4;q[11]=3;q[12]=5;q[13]=2;q[14]=2;q[15]=2;q[16]=2;q[17]=1;
	q[18]=0;q[19]=4;q[20]=1;q[21]=5;q[22]=0;q[23]=3;q[24]=0;q[25]=1;q[26]=2;
	q[27]=0;q[28]=5;q[29]=1;q[30]=4;q[31]=2;q[32]=3;q[33]=2;q[34]=1;q[35]=3;
	q[36]=0;q[37]=3;q[38]=1;q[39]=2;q[40]=3;q[41]=3;q[42]=1;q[43]=1;q[44]=4;
	q[45]=0;q[46]=2;q[47]=1;q[48]=3;q[49]=1;q[50]=3;q[51]=3;q[52]=1;q[53]=5;
	ed[0]=q;ed[1]=q+9;ed[2]=q+18;ed[3]=q+27;ed[4]=q+36;ed[5]=q+45;
#ifdef FACETS_NEIGHBOR
	q=mne[4];
	q[0]=-5;q[1]=-6;q[2]=-7;q[3]=-8;
	q[4]=-1;q[5]=-2;q[6]=-3;q[7]=-4;
	q[8]=-6;q[9]=-5;q[10]=-2;q[11]=-1;
	q[12]=-8;q[13]=-7;q[14]=-4;q[15]=-3;
	q[16]=-5;q[17]=-8;q[18]=-3;q[19]=-2;
	q[20]=-7;q[21]=-6;q[22]=-1;q[23]=-4;
	ne[0]=q;ne[1]=q+4;ne[2]=q+8;ne[3]=q+12;ne[4]=q+16;ne[5]=q+20;
#endif
	nu[0]=nu[1]=nu[2]=nu[3]=nu[4]=nu[5]=4;
};

// Initializes an arbitrary test object using the add_vertex() and
// relconstruct() routines
inline void voronoicell::init_test() {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;p=0;

	/*add_vertex(1,-2,-1,5,1,3);
	add_vertex(0,-1,1,2,0,5);
	add_vertex(0,1,0,4,6,3,1);
	add_vertex(1,2,-1,0,2,6);
	add_vertex(-1,2,-1,6,2,5);
	add_vertex(-1,-2,-1,1,0,4);
	add_vertex(0,3,0,3,2,4);*/

	/*add_vertex(-2,2,-1,1,4,3);
	add_vertex(2,2,-1,2,5,0);
	add_vertex(2,-2,-1,3,6,1);
	add_vertex(-2,-2,-1,0,7,2);
	add_vertex(-1,1,0,5,7,0);
	add_vertex(1,1,1,6,4,1);
	add_vertex(1,-1,1,5,2,7);
	add_vertex(-1,-1,1,6,3,4);*/

	/*add_vertex(1,-2,-1,4,3,1);
	add_vertex(-1,-2,-1,5,4,0,2);
	add_vertex(-1,2,-1,3,6,1);
	add_vertex(1,2,-1,0,5,6,2);
	add_vertex(0,-1,1,0,1,5);
	add_vertex(0,0,0,6,3,4,1);
	add_vertex(0,1,1,3,5,2);*/

	/*add_vertex(-1,-1,-1,1,3,4);
	add_vertex(1,-1,-1,5,2,0);
	add_vertex(1,1,-1,3,1,6);
	add_vertex(-1,1,-1,7,0,2);
	add_vertex(-1,-1,1,8,5,0,7);
	add_vertex(1,-1,1,8,6,1,4);
	add_vertex(1,1,1,8,7,2,5);
	add_vertex(-1,1,1,8,4,3,6);
	add_vertex(0,0,2,6,5,4,7);*/

	/*add_vertex(1,-3,-1,1,6,5);
	add_vertex(-1,-3,-1,2,6,0);
	add_vertex(-3,0,-1,3,8,7,1);
	add_vertex(-1,3,-1,4,9,2);
	add_vertex(1,3,-1,5,9,3);
	add_vertex(3,0,-1,0,7,8,4);
	add_vertex(0,-2,1,0,1,7);
	add_vertex(0,-1,0,5,6,2,8);
	add_vertex(0,1,0,5,7,2,9);
	add_vertex(0,2,1,4,8,3);*/

	/*add_vertex(-1,-3,-1,12,8,1,7);
	add_vertex(1,-3,-1,0,8,12,2);
	add_vertex(3,-1,-1,1,13,9,3);
	add_vertex(3,1,-1,2,9,13,4);
	add_vertex(1,3,-1,3,14,10,5);
	add_vertex(-1,3,-1,4,10,14,6);
	add_vertex(-3,1,-1,5,15,11,7);
	add_vertex(-3,-1,-1,6,11,15,0);
	add_vertex(0,-2,1,12,1,0);
	add_vertex(2,0,1,13,3,2);
	add_vertex(0,2,1,5,4,14);
	add_vertex(-2,0,1,6,15,7);
	add_vertex(0,-1,0.5,16,1,8,0);
	add_vertex(1,0,0.5,16,3,9,2);
	add_vertex(0,1,0.5,16,5,10,4);
	add_vertex(-1,0,0.5,16,7,11,6);
	add_vertex(0,0,0,14,13,12,15);*/

	/*add_vertex(2,-3,-1,1,4,3);
	add_vertex(-2,-3,-1,2,4,0);
	add_vertex(-2,3,-1,3,7,1);
	add_vertex(2,3,-1,0,6,2);
	add_vertex(0,-2,0,1,5,0);
	add_vertex(0,1,0,7,6,4);
	add_vertex(1,2,1,3,5,7);
	add_vertex(-1,2,1,2,6,5);*/

	/*add_vertex(3,-2,-1,1,3,2);
	add_vertex(-3,-2,-1,2,4,0);
	add_vertex(0,4,-1,0,5,1);
	add_vertex(1.5,-1,0,6,0,7);
	add_vertex(-1.5,-1,0,8,7,1);
	add_vertex(0,2,0,2,6,8);
	add_vertex(0.75,0.5,0,5,3,9);
	add_vertex(0,-1,0,9,3,4);
	add_vertex(-0.75,0.5,0,5,9,4);
	add_vertex(0,0,1,6,7,8);*/

	add_vertex(0,0,0,2,1,3);
	add_vertex(1,0,1,0,2,3);
	add_vertex(1,1,0,1,0,3);
	add_vertex(2,0,0,0,1,2,4,6);
	add_vertex(3,1,0,5,8,6,3);
	add_vertex(3,2,0,4);
	add_vertex(4,0,0,8,7,3,4);
	add_vertex(5,0,0,6);
	add_vertex(4,1,0,4,6);

	relconstruct();
};

// Adds an order 1 vertex to the memory structure, and specifies its edge
void voronoicell::add_vertex(f_point x,f_point y,f_point z,int a) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=1;
	if (mem[1]==mec[1]) addmemory(1);
	int *q=mep[1]+3*mec[1]++;ed[p]=q;
	q[0]=a;q[2]=p++;
};

// Adds an order 2 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(f_point x,f_point y,f_point z,int a,int b) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=2;
	if (mem[2]==mec[2]) addmemory(2);
	int *q=mep[2]+5*mec[2]++;ed[p]=q;
	q[0]=a;q[1]=b;q[4]=p++;
};

// Adds an order 3 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(f_point x,f_point y,f_point z,int a,int b,int c) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=3;
	if (mem[3]==mec[3]) addmemory(3);
	int *q=mep[3]+7*mec[3]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[6]=p++;
};

// Adds an order 4 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(f_point x,f_point y,f_point z,int a,int b,int c,int d) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=4;
	if (mem[4]==mec[4]) addmemory(4);
	int *q=mep[4]+9*mec[4]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[8]=p++;
};

// Adds an order 5 vertex to the memory structure, and specifies its edges
void voronoicell::add_vertex(f_point x,f_point y,f_point z,int a,int b,int c,int d,int e) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=5;
	if (mem[5]==mec[5]) addmemory(5);
	int *q=mep[5]+11*mec[5]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[4]=e;q[10]=p++;
};

// Checks that the relational table of the Voronoi cell is accurate, and prints
// out any errors. This algorithm is O(p), so running it every time the plane
// routine is called will result in a significant slowdown.
inline void voronoicell::relcheck() {
	int i,j;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if (ed[ed[i][j]][ed[i][nu[i]+j]]!=i) cout << "Relational error at point " << i << ", edge " << j << "." << endl;
		}
	}
};

// This routine checks for any two vertices that are connected by more than one
// edge. The plane algorithm is designed so that this should not happen, so any
// occurrences are most likely errors. Note that the routine is O(p), so
// running it every time the plane routine is called will result in a significant
// slowdown. 
inline void voronoicell::duplicatecheck() {
	int i,j,k;
	for(i=0;i<p;i++) {
		for(j=1;j<nu[i];j++) {
			for(k=0;k<j;k++) {
				if (ed[i][j]==ed[i][k]) cout << "Duplicate edges: (" << i << "," << j << ") and (" << i << "," << k << ") [" << ed[i][j] << "]" << endl;
			}
		}
	}
};

// Constructs the relational table if the edges have been specified
inline void voronoicell::relconstruct() {
	int i,j,k,l;
	for(i=0;i<p;i++) for(j=0;j<nu[i];j++) {
#ifdef FACETS_RADIUS
		mne[i][j]=-1;
#endif
		k=ed[i][j];
		l=0;
		while(ed[k][l]!=i) {
			l++;
			if (l==nu[k]) throw fatal_error("Relation table construction failed");
		}
		ed[i][nu[i]+j]=l;
	}
};

// Cuts the Voronoi cell by a particle whose center is at a separation of
// (x,y,z) from the cell center. The value of rsq should be initially set to
// x*x+y*y+z*z.
#ifdef FACETS_NEIGHBOR
bool voronoicell::nplane(f_point x,f_point y,f_point z,f_point rsq,int p_id) {
#else
bool voronoicell::plane(f_point x,f_point y,f_point z,f_point rsq) {
#endif
	int count=0,i,j,k,up=0,lp=0,tp,cp,qp=1,rp,stack=0;stack2=0;
	int us=0,ls=0,qs,iqs,cs,uw,qw=0,lw,tw;
	int *edp,*emp;
	f_point u,l,t,r,q;bool complicatedsetup=false,newdoubleedge=false,doubleedge=false;

	//Initialize the safe testing routine
	sure.init(x,y,z,rsq);

	//Test approximately sqrt(n)/4 points for their proximity to the plane
	//and keep the one which is closest
	uw=sure.test(up,u);t=abs(u);
	tw=qp=1;rp=p>>3;
	while(tw<rp) {
		qw=sure.test(qp,q);
		r=abs(q);
		if(r<t) {up=qp;u=q;t=r;uw=qw;}
		tw+=qp++;
	}
	lp=up;lw=uw;l=u;

	// Starting from an initial guess, we now move from vertex to vertex,
	// to try and find an edge which intersects the cutting plane,
	// or a vertex which is on the plane
	try {
		if(uw==1) {

			// The test point is within the cutting space
			do {

				// If we have been around this loop more times
				// than there are points, there's a floating
				// point problem, so we'll bail out 
				if (++count>=p) throw true;
				
				// Test all the neighbors of the current point
				// and find the one which is closest to the
				// plane
				u=l;up=lp;uw=lw;
				for(i=0;i<nu[up];i++) {
					tp=ed[up][i];
					tw=sure.test(tp,t);
					if(t<l) {l=t;lw=tw;lp=tp;us=i;}
				}

				// If we couldn't find a point and the object
				// is convex, then the whole cell must be
				// within the cutting space, so there's nothing
				// left
				if (lp==up) {
					cerr << "Failed to find intersection" << endl;
					return false;
				}
			} while (lw==1);
			ls=ed[up][nu[up]+us];

			// If the last point in the iteration is within the
			// plane, we need to do the complicated setup
			// routine. Otherwise, we use the regular iteration.
			if (lw==0) {
				up=lp;
				complicatedsetup=true;
			} else complicatedsetup=false;
		} else if (uw==-1) {

			// The test point is outside of the cutting space
			do {

				// If we have been around this loop more times
				// than there are points, there's a floating
				// point problem, so we'll bail out 
				if (++count>=p) throw true;
				
				// Test all the neighbors of the current point
				// and find the one which is closest to the
				// plane
				l=u;lp=up;lw=uw;
				for(i=0;i<nu[lp];i++) {
					tp=ed[lp][i];
					tw=sure.test(tp,t);
					if(t>u) {u=t;uw=tw;up=tp;ls=i;}
				}

				// If we couldn't find a point and the object
				// is convex, then the whole cell must be
				// outside the cutting space, so it's not
				// intersected at all 
				if (up==lp) return true;
			} while (uw==-1);
			us=ed[lp][nu[lp]+ls];
			complicatedsetup=(uw!=1);
		} else {

			// Our original test point was on the plane, so we
			// automatically head for the complicated setup
			// routine
			complicatedsetup=true;
		}
	}
	catch(bool except) {

		// This routine is a fall-back, in case floating point errors
		// cause the usual search routine to fail. In the fall-back
		// routine, we just test every edge to find one straddling
		// the plane.
		cerr << "Bailed out of convex calculation\n";
		for(qp=0;qp<p;qp++) {
			qw=sure.test(qp,q);
			if (qw==1) {

				// The point is inside the cutting space. Now
				// see if we can find a neighbor which isn't.
				for(us=0;us<nu[qp];us++) {
					lp=ed[qp][us];
					if(lp<qp) {
						lw=sure.test(lp,l);
						if (lw!=1) break;
					}
				}
				if(us<nu[qp]) {
					up=qp;
					if(lw==0) {
						complicatedsetup=true;
					} else {
						complicatedsetup=false;
						u=q;
						ls=ed[up][nu[up]+us];
					}
					break;
				}
			} else if (qw==-1) {

				// The point is outside the cutting space. See
				// if we can find a neighbor which isn't.
				for(ls=0;ls<nu[qp];ls++) {
					up=ed[qp][ls];
					if(up<qp) {
						uw=sure.test(up,u);
						if (uw!=-1) break;
					}
				}
				if(ls<nu[qp]) {
					if(uw==0) {
						up=qp;
						complicatedsetup=true;
					} else {
						complicatedsetup=false;
						lp=qp;l=q;
						us=ed[lp][nu[lp]+ls];
					}
					break;
				}
			} else {
				
				// The point is in the plane, so we just
				// proceed with the complicated setup routine
				up=qp;
				complicatedsetup=true;
				break;
			}
		}
		if(qp==p) return qw==-1?true:false;
	}

	// We're about to add the first point of the new facet. In either
	// routine, we have to add a point, so first check there's space for
	// it.
	if(p==currentvertices) addmemory_vertices();

	if (complicatedsetup) {

		// The search algorithm found a point which is on the cutting
		// plane. We leave that point in place, and create a new one at
		// the same location.
		pts[3*p]=pts[3*up];
		pts[3*p+1]=pts[3*up+1];
		pts[3*p+2]=pts[3*up+2];
		
		// Search for a collection of edges of the test vertex which
		// are outside of the cutting space. Begin by testing the
		// zeroth edge. 
		i=0;
		lp=ed[up][0];
		lw=sure.test(lp,l);
		if(lw!=-1) {

			// The first edge is either inside the cutting space,
			// or lies within the cutting plane. Test the edges
			// sequentially until we find one that is outside.
			rp=lw;
			do {
				i++;

				// If we reached the last edge with no luck
				// then all of the vertices are inside
				// or on the plane, so the cell is completely
				// deleted
				if (i==nu[up]) return false;
				lp=ed[up][i];
				lw=sure.test(lp,l);
			} while (lw!=-1);
			j=i+1;

			// We found an edge outside the cutting space. Keep
			// moving through these edges until we find one that's
			// inside or on the plane.
			while(j<nu[up]) {
				lp=ed[up][j];
				lw=sure.test(lp,l);
				if (lw!=-1) break;
				j++;
			}

			// Compute the number of edges for the new vertex. In
			// general it will be the number of outside edges
			// found, plus two. But we need to recognize the
			// special case when all but one edge is outside, and
			// the remaining one is on the plane. For that case we
			// have to reduce the edge count by one to prevent
			// doubling up.
			if(j==nu[up]&&i==1&&rp==0) {
				nu[p]=nu[up];
				doubleedge=true;
			} else nu[p]=j-i+2;
			k=1;

			// Add memory for the new vertex if needed, and
			// initialize
			while (nu[p]>=currentvertexorder) addmemory_vorder();
			if (mec[nu[p]]==mem[nu[p]]) addmemory(nu[p]);
#ifdef FACETS_NEIGHBOR
			ne[p]=mne[nu[p]]+nu[p]*mec[nu[p]];
#endif
			ed[p]=mep[nu[p]]+(2*nu[p]+1)*mec[nu[p]]++;
			ed[p][2*nu[p]]=p;

			// Copy the edges of the original vertex into the new
			// one. Delete the edges of the original vertex, and
			// update the relational table.
			us=vor_down(i,up);
			while(i<j) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			qs=i==nu[up]?0:i;

		} else {

			// In this case, the zeroth edge is outside the cutting
			// plane. Begin by searching backwards from the last
			// edge until we find an edge which isn't outside.
			i=nu[up]-1;
			lp=ed[up][i];
			lw=sure.test(lp,l);
			while(lw==-1) {
				i--;

				// If i reaches zero, then we have a point in
				// the plane all of whose edges are outside
				// the cutting space, so we just exit
				if (i==0) return true;
				lp=ed[up][i];
				lw=sure.test(lp,l);
			}

			// Now search forwards from zero
			j=1;
			qp=ed[up][j];
			qw=sure.test(qp,q);
			while(qw==-1) {
				j++;
				qp=ed[up][j];
				qw=sure.test(qp,l);
			}

			// Compute the number of edges for the new vertex. In
			// general it will be the number of outside edges
			// found, plus two. But we need to recognize the
			// special case when all but one edge is outside, and
			// the remaining one is on the plane. For that case we
			// have to reduce the edge count by one to prevent
			// doubling up.
			if (i==j&&qw==0) {
				doubleedge=true;
				nu[p]=nu[up];
			} else {
				nu[p]=nu[up]-i+j+1;
			}

			// Add memory to store the vertex if it doesn't exist
			// already
			k=1;
			while(nu[p]>=currentvertexorder) addmemory_vorder();
			if (mec[nu[p]]==mem[nu[p]]) addmemory(nu[p]);

			// Copy the edges of the original vertex into the new
			// one. Delete the edges of the original vertex, and
			// update the relational table.
#ifdef FACETS_NEIGHBOR
			ne[p]=mne[nu[p]]+nu[p]*mec[nu[p]];
#endif
			ed[p]=mep[nu[p]]+(2*nu[p]+1)*mec[nu[p]]++;
			ed[p][2*nu[p]]=p;
			us=i++;
			while(i<nu[up]) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			i=0;
			while(i<j) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			qs=j;
		}
		
		// Add this point to the auxiliary delete stack
		if (stack2==currentdeletesize2) addmemory_ds2();
		ds2[stack2++]=up;

		// Look at the edges on either side of the group that was
		// detected. We're going to commence facet computation by
		// moving along one of them. We are going to end up coming back
		// along the other one.
		cs=k;
		qp=up;q=u;
		i=ed[up][us];
		us=ed[up][nu[up]+us];
		up=i;
		ed[qp][2*nu[qp]]=-p;

	} else {

		// The search algorithm found an intersected edge between the
		// points lp and up. Create a new vertex between them which
		// lies on the cutting plane. Since u and l differ by at least
		// the tolerance, this division should never screw up.
		if (stack==currentdeletesize) addmemory_ds();
		ds[stack++]=up;
		r=1/(u-l);
		pts[3*p]=(pts[3*lp]*u-pts[3*up]*l)*r;
		pts[3*p+1]=(pts[3*lp+1]*u-pts[3*up+1]*l)*r;
		pts[3*p+2]=(pts[3*lp+2]*u-pts[3*up+2]*l)*r;

		// This point will always have three edges. Connect one of them
		// to lp.
		nu[p]=3;
		if (mec[3]==mem[3]) addmemory(3);
#ifdef FACETS_NEIGHBOR
		ne[p]=mne[3]+3*mec[3];
		ne[p][0]=p_id;
		ne[p][1]=ne[lp][ls];
		ne[p][2]=ne[qp][qs];
#endif
		ed[p]=mep[3]+7*mec[3]++;
		ed[p][6]=p;
		ed[up][us]=-1;
		ed[lp][ls]=p;
		ed[lp][nu[lp]+ls]=1;
		ed[p][1]=lp;
		ed[p][nu[p]+1]=ls;
		cs=2;

		// Set the direction to move in
		qs=vor_up(us,up);
		qp=up;q=u;
	}

	// When the code reaches here, we have initialized the first point, and
	// we have a direction for moving it to construct the rest of the facet
	cp=p;rp=p;p++;
	while(qp!=up||qs!=us) {

		// We're currently tracing round an intersected facet. Keep
		// moving around it until we find a point or edge which
		// intersects the plane. 
		lp=ed[qp][qs];
		lw=sure.test(lp,l);
		
		if (lw==1) {

			// The point is still in the cutting space. Just add it
			// to the delete stack and keep moving.
			if (stack==currentdeletesize) addmemory_ds();
			qs=vor_up(ed[qp][nu[qp]+qs],lp);
			qp=lp;
			q=l;
			ds[stack++]=qp;
		
		} else if (lw==-1) {

			// The point is outside of the cutting space, so we've
			// found an intersected edge. Introduce a regular point
			// at the point of intersection. Connect it to the
			// point we just tested. Also connect it to the previous
			// new point in the facet we're constructing.
			if(p==currentvertices) addmemory_vertices();
			r=1/(q-l);
			pts[3*p]=(pts[3*lp]*q-pts[3*qp]*l)*r;
			pts[3*p+1]=(pts[3*lp+1]*q-pts[3*qp+1]*l)*r;
			pts[3*p+2]=(pts[3*lp+2]*q-pts[3*qp+2]*l)*r;
			nu[p]=3;
			if (mec[3]==mem[3]) addmemory(3);
			ls=ed[qp][qs+nu[qp]];
#ifdef FACETS_NEIGHBOR
			ne[p]=mne[3]+3*mec[3];
			ne[p][0]=p_id;
			ne[p][1]=ne[lp][ls];
			ne[p][2]=ne[qp][qs];
#endif
			ed[p]=mep[3]+7*mec[3]++;
			ed[p][6]=p;
			ed[lp][ls]=p;
			ed[lp][nu[lp]+ls]=1;
			ed[p][1]=lp;
			ed[p][0]=cp;
			ed[p][nu[p]+1]=ls;
			ed[p][nu[p]]=cs;
			ed[cp][cs]=p;
			ed[cp][nu[cp]+cs]=0;
			ed[qp][qs]=-1;
			qs=vor_up(qs,qp);
			cp=p++;
			cs=2;
		} else {

			// We've found a point which is on the cutting plane.
			// We're going to introduce a new point right here, but
			// first we need to figure out the number of edges it
			// has.
			if(p==currentvertices) addmemory_vertices();
			
			// If the previous vertex detected a double edge, our
			// new vertex will have one less edge.
			k=doubleedge?0:1;
			qs=ed[qp][nu[qp]+qs];
			qp=lp;
			iqs=qs;

			// Start testing the edges of the current point until
			// we find one which isn't outside the cutting space
			do {
				k++;
				qs=vor_up(qs,qp);
				lp=ed[qp][qs];
				lw=sure.test(lp,l);
			} while (lw==-1);
			
			// Now we need to find out whether this marginal vertex
			// we are on has been visited before, because if that's
			// the case, we need to add vertices to the existing
			// new vertex, rather than creating a fresh one. We also
			// need to figure out whether we're in a case where we
			// might be creating a duplicate edge.
			j=-ed[qp][2*nu[qp]];
	 		if(qp==up&&qs==us) {

				// If we're heading into the final part of the
				// new facet, then we never worry about the
				// duplicate edge calculation.
				newdoubleedge=false;
				if(j>0) k+=nu[j];
			} else {
				if (j>0) {

					// This vertex was visited before, so
					// count those vertices to the ones we
					// already have.
					k+=nu[j];
					
					// The only time when we might make a
					// duplicate edge is if the point we're
					// going to move to next is also a
					// marginal point, so test for that
					// first.
					if(lw==0) {

						// Now see whether this marginal point
						// has been visited before.
						i=-ed[lp][2*nu[lp]];
						if(i>0) {

							// Now see if the last edge of that other
							// marginal point actually ends up here.
							if(ed[i][nu[i]-1]==j) {
								newdoubleedge=true;
								k-=1;
							} else newdoubleedge=false;
						} else {

							// That marginal point hasn't been visited
							// before, so we probably don't have to worry
							// about duplicate edges, except in the
							// case when that's the way into the end
							// of the facet, because that way always creates
							// an edge.
							if (j==rp&&lp==up&&ed[qp][nu[qp]+qs]==us) {
								newdoubleedge=true;
								k-=1;
							} else newdoubleedge=false;
						}
					} else newdoubleedge=false;
				} else {

					// The vertex hasn't been visited
					// before, but let's see if it's
					// marginal
					if(lw==0) {

						// If it is, we need to check
						// for the case that it's a
						// small branch, and that we're
						// heading right back to where
						// we came from
						i=-ed[lp][2*nu[lp]];
						if(i==cp) {
							newdoubleedge=true;
							k-=1;
						} else newdoubleedge=false;
					} else newdoubleedge=false;
				}
			}
			
			// k now holds the number of edges of the new vertex
			// we are forming. Add memory for it if it doesn't exist
			// already.
			while(k>=currentvertexorder) addmemory_vorder();
			if (mec[k]==mem[k]) addmemory(k);
			
			// Now create a new vertex with order k, or augment
			// the existing one.
			if(j>0) {

				// If we're augmenting a vertex but we don't
				// actually need any more edges, just skip this
				// routine to avoid memory confusion
				if(nu[j]!=k) {

					// Allocate memory and copy the edges
					// of the previous instance into it
					edp=mep[k]+(2*k+1)*mec[k]++;
					i=0;
					while(i<nu[j]) {
						edp[i]=ed[j][i];
						edp[k+i]=ed[j][nu[j]+i];
						i++;
					}
					edp[2*k]=j;

					// Remove the previous instance with
					// fewer vertices from the memory
					// structure
					mec[nu[j]]--;
					emp=mep[nu[j]]+(2*nu[j]+1)*mec[nu[j]];
					if(emp!=ed[j]) {
						for(lw=0;lw<=2*nu[j];lw++) ed[j][lw]=emp[lw];
						ed[emp[2*nu[j]]]=ed[j];
					}
					ed[j]=edp;
				} else i=nu[j];
			} else {

				// Allocate a new vertex of order k
				ed[p]=mep[k]+(2*k+1)*mec[k]++;
				ed[p][2*k]=p;
				if (stack2==currentdeletesize2) addmemory_ds2();
				ds2[stack2++]=qp;
				pts[3*p]=pts[3*qp];
				pts[3*p+1]=pts[3*qp+1];
				pts[3*p+2]=pts[3*qp+2];
				ed[qp][2*nu[qp]]=-p;
				j=p++;
				i=0;
			}
			nu[j]=k;

			// Unless the previous case was a double edge, connect
			// the first available edge of the new vertex to the 
			// last one in the facet
			if (!doubleedge) {
				ed[j][i]=cp;
				ed[j][nu[j]+i]=cs;
				ed[cp][cs]=j;
				ed[cp][nu[cp]+cs]=i;
				i++;
			}

			// Copy in the edges of the underlying vertex,
			// and do one less if this was a double edge 
			qs=iqs;
			while(i<(newdoubleedge?k:k-1)) {
				qs=vor_up(qs,qp);
				lp=ed[qp][qs];ls=ed[qp][nu[qp]+qs];
				ed[j][i]=lp;
				ed[j][nu[j]+i]=ls;
				ed[lp][ls]=j;
				ed[lp][nu[lp]+ls]=i;
				ed[qp][qs]=-1;
				i++;
			}
			qs=vor_up(qs,qp);
			cs=i;
			cp=j;

			// Update the doubleedge flag, to pass it
			// to the next instance of this routine
			doubleedge=newdoubleedge;
		}
	}

	// Connect the final created vertex to the initial one
	ed[cp][cs]=rp;
	ed[rp][0]=cp;
	ed[cp][nu[cp]+cs]=0;
	ed[rp][nu[rp]+0]=cs;

	// Delete points: first, remove any duplicates
	i=0;
	while(i<stack) {
		j=ds[i];
		if(ed[j][nu[j]]!=-1) {
			ed[j][nu[j]]=-1;
			i++;
		} else ds[i]=ds[--stack];
	}
	
	// Add the points in the auxiliary delete stack,
	// and reset their back pointers
	for(i=0;i<stack2;i++) {
		j=ds2[i];
		ed[j][2*nu[j]]=j;
		if(ed[j][nu[j]]!=-1) {
			ed[j][nu[j]]=-1;
			if (stack==currentdeletesize) addmemory_ds();
			ds[stack++]=j;
		}
	}
	
	// Scan connections and add in extras
	for(i=0;i<stack;i++) {
		cp=ds[i];
		for(j=0;j<nu[cp];j++) {
			qp=ed[cp][j];
			if(qp!=-1) {
				if (ed[qp][nu[qp]]!=-1) {
					if (stack==currentdeletesize) addmemory_ds();
					ds[stack++]=qp;
					ed[qp][nu[qp]]=-1;
				}
			}
		}
	}

	// Delete them from the array structure
	while(stack>0) {
		while(ed[--p][nu[p]]==-1) {
			j=nu[p];
			mec[j]--;
			for(i=0;i<=2*j;i++) ed[p][i]=(mep[j]+(2*j+1)*mec[j])[i];
			ed[ed[p][2*j]]=ed[p];
		}
		qp=ds[--stack];
		if (qp<p) {

			// Vertex management
			pts[3*qp]=pts[3*p];
			pts[3*qp+1]=pts[3*p+1];
			pts[3*qp+2]=pts[3*p+2];

			// Memory management
			j=nu[qp];
			mec[j]--;
			for(i=0;i<=2*j;i++) ed[qp][i]=(mep[j]+(2*j+1)*mec[j])[i];
			ed[ed[qp][2*j]]=ed[qp];

			// Edge management
			ed[qp]=ed[p];
			nu[qp]=nu[p];
			for(i=0;i<nu[qp];i++) {
				if (ed[qp][i]==-1) throw fatal_error("fishy");
				ed[ed[qp][i]][ed[qp][nu[qp]+i]]=qp;
			}
			ed[qp][2*nu[qp]]=qp;
		} else p++;
	}

	// Check for any vertices of zero order
	if (mec[0]>0) throw fatal_error("Zero order vertex formed");

	// Collapse any order 2 vertices and exit	
	return collapseorder2();
};

// During the creation of a new facet in the plane routine, it is possible
// that some order 2 vertices may arise. This routine removes them.
// Suppose an order 2 vertex joins c and d. If there's a edge between
// c and d already, then the order 2 vertex is just removed; otherwise,
// the order 2 vertex is removed and c and d are joined together directly.
// It is possible this process will create order 2 or order 1 vertices,
// and the routine is continually run until all of them are removed.
inline bool voronoicell::collapseorder2() {
	if(!collapseorder1()) return false;
	int a,b,i,j,k,l;
	while(mec[2]>0) {

		// Pick a order 2 vertex and read in its edges
		i=--mec[2];
		j=mep[2][5*i];k=mep[2][5*i+1];
		if (j==k) {
			cerr << "Order two vertex joins itself" << endl;
			return false;
		}

		// Scan the edges of j to see if joins k
		for(l=0;l<nu[j];l++) {
			if(ed[j][l]==k) break;
		}

		// If j doesn't already join k, join them together.
		// Otherwise delete the connection to the current
		// vertex from j and k.
		a=mep[2][5*i+2];b=mep[2][5*i+3];i=mep[2][5*i+4];
		if(l==nu[j]) {
			ed[j][a]=k;
			ed[k][b]=j;
			ed[j][nu[j]+a]=b;
			ed[k][nu[k]+b]=a;
		} else {
			if (!delete_connection(j,a)) return false;
			if (!delete_connection(k,b)) return false;
		}

		// Compact the memory
		--p;
		if(p!=i) {
			pts[3*i]=pts[3*p];
			pts[3*i+1]=pts[3*p+1];
			pts[3*i+2]=pts[3*p+2];
			for(k=0;k<nu[p];k++) ed[ed[p][k]][ed[p][nu[p]+k]]=i;
			ed[i]=ed[p];
			nu[i]=nu[p];
			ed[i][2*nu[i]]=i;
		}

		// Collapse any order 1 vertices if they were created
		if(!collapseorder1()) return false;
	}
	return true;
};

// Order 1 vertices can potentially be created during the order 2 collapse
// routine. This routine removes them.
inline bool voronoicell::collapseorder1() {
	int i,j,k;
	while(mec[1]>0) {
		cerr << "Order one collapse" << endl;
		i=--mec[1];
		j=mep[1][3*i];k=mep[1][3*i+1];
		i=mep[1][3*i+2];
		if(!delete_connection(j,k)) return false;
		--p;
		if(p!=i) {
			pts[3*i]=pts[3*p];
			pts[3*i+1]=pts[3*p+1];
			pts[3*i+2]=pts[3*p+2];
			for(k=0;k<nu[p];k++) ed[ed[p][k]][ed[p][nu[p]+k]]=i;
			ed[i]=ed[p];
			nu[i]=nu[p];
			ed[i][2*nu[i]]=i;
		}
	}
	return true;
};

// This routine deletes the kth edge of vertex j and reorganizes the memory
inline bool voronoicell::delete_connection(int j,int k) {
	int i=nu[j]-1,l,*edp,*edd,m;
	if(i<1) {
		cout << "Zero order vertex formed" << endl;
		return false;
	}
	if(mec[i]==mem[i]) addmemory(i);
	edp=mep[i]+(2*i+1)*mec[i]++;
	edp[2*i]=j;
	for(l=0;l<k;l++) {
		edp[l]=ed[j][l];
		edp[l+i]=ed[j][l+nu[j]];
	}
	while(l<i) {
		m=ed[j][l+1];
		edp[l]=m;
		k=ed[j][l+nu[j]+1];
		edp[l+i]=k;
		ed[m][nu[m]+k]--;
		l++;
	}
	edd=mep[nu[j]]+(2*nu[j]+1)*--mec[nu[j]];
	for(l=0;l<=2*nu[j];l++) ed[j][l]=edd[l];
	ed[edd[2*nu[j]]]=edd;
	ed[j]=edp;
	nu[j]=i;
	return true;
};

// Cuts a Voronoi cell using the influence of a particle at (x,y,z), first
// calculating the modulus squared of this vector before passing it to the
// routine above
inline bool voronoicell::plane(f_point x,f_point y,f_point z) {
	f_point rsq=x*x+y*y+z*z;
	return plane(x,y,z,rsq);
};

// For the neighbor-tracking version of the code, an extra version of the plane
// routine is provided that doesn't require passing a plane ID. It just makes
// up the plane ID to be zero. An nplane routine that works without passing
// the modulus squared is also provided.
#ifdef FACETS_NEIGHBOR
inline bool voronoicell::plane(f_point x,f_point y,f_point z,f_point rsq) {
	return nplane(x,y,z,rsq,0);
};
inline bool voronoicell::nplane(f_point x,f_point y,f_point z,int p_id) {
	f_point rsq=x*x+y*y+z*z;
	return nplane(x,y,z,rsq,p_id);
};	
#endif

// Simple functions for moving around the edges of a given Voronoi vertex
inline int voronoicell::vor_up(int a,int p) {
	return a==nu[p]-1?0:a+1;
};
inline int voronoicell::vor_down(int a,int p) {
	return a==0?nu[p]-1:a-1;
};

// Calculates the volume of a Voronoi cell
inline f_point voronoicell::volume() {
	const f_point fe=1/48.0;
	f_point vol=0;
	int i,j,k,l,m,n;
	f_point ux,uy,uz,vx,vy,vz,wx,wy,wz;
	for(i=1;i<p;i++) {
		ux=pts[0]-pts[3*i];
		uy=pts[1]-pts[3*i+1];
		uz=pts[2]-pts[3*i+2];
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				ed[i][j]=-1-k;
				l=vor_up(ed[i][nu[i]+j],k);
				vx=pts[3*k]-pts[0];
				vy=pts[3*k+1]-pts[1];
				vz=pts[3*k+2]-pts[2];
				m=ed[k][l];ed[k][l]=-1-m;
				while(m!=i) {
					n=vor_up(ed[k][nu[k]+l],m);
					wx=pts[3*m]-pts[0];
					wy=pts[3*m+1]-pts[1];
					wz=pts[3*m+2]-pts[2];
					vol+=ux*vy*wz+uy*vz*wx+uz*vx*wy-uz*vy*wx-uy*vx*wz-ux*vz*wy;
					k=m;l=n;vx=wx;vy=wy;vz=wz;
					m=ed[k][l];ed[k][l]=-1-m;
				}
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Volume routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}

	return vol*fe;
};

// Computes the maximum radius squared
inline f_point voronoicell::maxradsq() {
	int i;f_point r,s;
	r=pts[0]*pts[0]+pts[1]*pts[1]+pts[2]*pts[2];
	for(i=3;i<3*p;i+=3) {
		s=pts[i]*pts[i]+pts[i+1]*pts[i+1]+pts[i+2]*pts[i+2];
		if(s>r) r=s;
	}
	return r;
};

// Outputs the edges of the Voronoi cell (in POV-Ray format) to an open file
// stream, displacing the cell by an amount (x,y,z)
inline void voronoicell::dumppov(ostream &of,f_point x,f_point y,f_point z) {
	int i,j,k;f_point ux,uy,uz;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];uz=z+0.5*pts[3*i+2];
		of << "sphere{<" << ux << "," << uy << "," << uz << ">,r}" << endl;
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k<i) of << "cylinder{<" << ux << "," << uy << "," << uz << ">,<" << x+0.5*pts[3*k] << "," << y+0.5*pts[3*k+1] << "," << z+0.5*pts[3*k+2] << ">,r}" << endl;
		}
	}
};

// Outputs the edges of the Voronoi cell (in gnuplot format) to an open file
// stream, displacing the cell by an amount (x,y,z)
inline void voronoicell::dumpgnuplot(ostream &of,f_point x,f_point y,f_point z) {
	int i,j,k;f_point ux,uy,uz;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];uz=z+0.5*pts[3*i+2];
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (ed[i][j]<i) of << ux << " " << uy << " " << uz << endl << x+0.5*pts[3*k] << " " << y+0.5*pts[3*k+1] << " " << z+0.5*pts[3*k+2] << endl << endl << endl;
		}
	}
};

// Outputs the Voronoi cell in the POV mesh2 format
void voronoicell::dumppovmesh(ostream &of,f_point x,f_point y,f_point z) {
	int i,j,k,l,m,n;
	of << "mesh2 {" << endl << "vertex_vectors {" << endl << p << "," << endl;
	for(i=0;i<p;i++) {
		of << "<" << x+0.5*pts[3*i] << "," << y+0.5*pts[3*i+1] << "," << z+0.5*pts[3*i+2] << ">," << endl;
	}
	of << "}" << endl << "face_indices {" << endl << 2*(p-2) << "," << endl;
	for(i=1;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				ed[i][j]=-1-k;
				l=vor_up(ed[i][nu[i]+j],k);
				m=ed[k][l];ed[k][l]=-1-m;
				while(m!=i) {
					n=vor_up(ed[k][nu[k]+l],m);
					of << "<" << i << "," << k << "," << m << ">" << endl;
					k=m;l=n;
					m=ed[k][l];ed[k][l]=-1-m;
				}
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Mesh routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
	of << "}" << endl << "inside_vector <0,0,1>" << endl << "}" << endl;
};

// Randomly perturbs the points in the Voronoi cell by an amount r
inline void voronoicell::perturb(f_point r) {
	for(int i=0;i<3*p;i++) {
		pts[i]+=(2*double(rand())/RAND_MAX-1)*r;
	}
};

//Initialises the suretest class and creates a buffer for dubious points
suretest::suretest() : currentdubious(initdubious) {
	sn=new int[2*currentdubious];
};

// Suretest destructor to free memory allocation
suretest::~suretest() {
	delete [] sn;
};

// Sets up the suretest class with a particular test plane, and removes
// any special cases from the table
inline void suretest::init(f_point x,f_point y,f_point z,f_point rsq) {
	sc=0;px=x;py=y;pz=z;prsq=rsq;
};

inline int suretest::test(int n,f_point &ans) {
	ans=px*p[3*n]+py*p[3*n+1]+pz*p[3*n+2]-prsq;
	if(ans>tolerance2) {
		return 1;
	} else if(ans<-tolerance2) {
		return -1;
	} else {
		int i;
		for(i=0;i<sc;i+=2) if(sn[i]==n) return sn[i+1];
		if (sc==2*currentdubious) {
			i=2*currentdubious;
			if (i>maxdubious) throw fatal_error("Dubious case buffer allocation exceeded absolute maximum");
			cerr << "Dubious cases buffer scaled up to " << i << endl;
			int *psn=new int[2*i];
			for(int j=0;j<2*currentdubious;j++) psn[j]=sn[j];
			delete [] sn;sn=psn;
		}
		sn[sc++]=n;
		sn[sc++]=ans>tolerance?1:(ans<-tolerance?-1:0);
		return sn[sc-1];
	}
};

// Prints the vertices, their edges, the relation table,
// and also notifies if any glaring memory errors are visible.
void voronoicell::edgeprint() {
	int j;
	for(int i=0;i<p;i++) {
		cout << i << " " << nu[i] << "  ";
		for(j=0;j<nu[i];j++) cout << " " << ed[i][j];
		cout << "    ";
		while(j<2*nu[i]) cout << " " << ed[i][j++];
		cout << "     " << ed[i][j];
		cout << " " << pts[3*i] << " " << pts[3*i+1] << " " << pts[3*i+2];
		cout << " " << ed[i];
		if (ed[i]>=mep[nu[i]]+mec[nu[i]]*(2*nu[i]+1)) cout << " Memory error";
		cout << endl;
	}
};

// Prints out a list of all the facets and their vertices. If the neighbor option
// is defined, it lists each cutting plane.
void voronoicell::facets(ostream &of) {
	int i,j,k,l,m;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
#ifdef FACETS_NEIGHBOR
				of << "(" << i << "," << ne[i][j] << ")";
#else
				of << i;
#endif
				ed[i][j]=-1-k;
				l=vor_up(ed[i][nu[k]+j],k);
				do {
#ifdef FACETS_NEIGHBOR
					of << " (" << k << "," << ne[k][l] << ")";
#else
					of << " " << k;
#endif
					m=ed[k][l];
					ed[k][l]=-1-m;
					l=vor_up(ed[k][nu[k]+l],m);
					k=m;
				} while (k!=i);
				of << endl;
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Facet statistics routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
}

// Overloaded versions of facets
inline void voronoicell::facets() {
	facets(cout);
};
inline void voronoicell::facets(char *filename) {
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	facets(of);
	of.close();
};

// Examines all the facets, and evaluates them by the number of vertices that
// they have
void voronoicell::facet_statistics(ostream &of) {
	int *stat,*pstat,currentfacetsize=initfacetsize,newc,maxf=0;
	stat=new int[currentfacetsize];
	int i,j,k,l,m,q;
	for(i=0;i<currentfacetsize;i++) stat[i]=0;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				q=1;
				ed[i][j]=-1-k;
				l=vor_up(ed[i][nu[k]+j],k);
				do {
					q++;
					m=ed[k][l];
					ed[k][l]=-1-m;
					l=vor_up(ed[k][nu[k]+l],m);
					k=m;
				} while (k!=i);
				if (q>=currentfacetsize) {
					newc=currentfacetsize*2;
					pstat=new int[newc];
					for(k=0;k<currentfacetsize;k++) pstat[k]=stat[k];
					while(k<newc) pstat[k]=0;
					delete [] stat;
					currentfacetsize=newc;
					stat=pstat;
				}
				stat[q]++;
				if (q>maxf) maxf=q;
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Facet statistics routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
	for(i=0;i<=maxf;i++) of << i << " " << stat[i] << endl; 
	delete [] stat;
};

// Overloaded versions of facet_statistics
inline void voronoicell::facet_statistics() {
	facet_statistics(cout);
};
inline void voronoicell::facet_statistics(char *filename) {
	ofstream of;
	of.open(filename,ofstream::out|ofstream::trunc);
	facet_statistics(of);
	of.close();
};
