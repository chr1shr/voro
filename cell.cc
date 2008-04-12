// Voronoi calculation code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : February 27th 2008

#include "cell.hh"

/** Constructs a Voronoi cell and sets up the initial memory. */
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
	pts=new fpoint[3*currentvertices];
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

/** The voronoicell destructor deallocates all the dynamic memory. */
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
#ifdef FACETS_NEIGHBOR
	delete [] mne;
	delete [] ne;
#endif
	delete [] ed;
	delete [] nu;
	delete [] pts;
};

/** Increases the memory storage for a particular vertex order, by increasing
 * the size of the of the mep and mne <em>(Neighbor version only)</em> arrays.
 * If the arrays already exist, their size is doubled; if they don't exist,
 * then new ones of size initnvertices are allocated. The routine also ensures
 * that the pointers in the ed array are updated, by making use of the back
 * pointers. For the cases where the back pointer has been temporarily
 * overwritten in the marginal vertex code, the auxiliary delete stack is
 * scanned to find out how to update the ed value.
 * \param[in] i The order of the vertex memory to be increased.*/
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
		int j=0,k,*l;
		mem[i]*=2;
		if (mem[i]>maxnvertices) throw fatal_error("Point memory allocation exceeded absolute maximum");
		cerr << "Order " << i << " vertex memory scaled up to " << mem[i] << endl;
		l=new int[s*mem[i]];
#ifdef FACETS_NEIGHBOR
		int *q,m=0;
		q=new int[i*mem[i]];
#endif
		while(j<s*mec[i]) {
			k=mep[i][j+2*i];
			if(k>=0) {
				ed[k]=l+j;
#ifdef FACETS_NEIGHBOR
				ne[k]=q+m;
#endif
			} else {
				int m;
				for(m=0;m<stack2;m++) {
					if(ed[ds2[m]]==mep[i]+j) {
						ed[ds2[m]]=l+j;
#ifdef FACETS_NEIGHBOR
						ne[ds2[m]]=q+m;
#endif
						break;
					}
				}
				if(m==stack2) throw fatal_error("Couldn't relocate dangling pointer");
				cerr << "Relocated dangling pointer" << endl;
			}
			for(k=0;k<s;k++,j++) l[j]=mep[i][j];
#ifdef FACETS_NEIGHBOR
			for(k=0;k<i;k++,m++) q[m]=mne[i][m];
#endif
		}
		delete [] mep[i];
		mep[i]=l;
#ifdef FACETS_NEIGHBOR
		delete [] mne[i];
		mne[i]=q;
#endif
	}
};

/** Doubles the maximum number of vertices allowed, by reallocating the ed, nu,
 * pts, and ne <em>(Neighbor version only)</em> arrays. If the allocation
 * exceeds the absolute maximum set in maxvertices, then the routine throws a
 * fatal error. */
void voronoicell::addmemory_vertices() {
	int i=2*currentvertices,j,**pp,*pnu;
	if (i>maxvertices) throw fatal_error("Vertex memory allocation exceeded absolute maximum");
	cerr << "Vertex memory scaled up to " << i << endl;
	fpoint *ppts;
	pp=new int*[i];
	for(j=0;j<currentvertices;j++) pp[j]=ed[j];
	delete [] ed;ed=pp;
#ifdef FACETS_NEIGHBOR
	pp=new int*[i];
	for(j=0;j<currentvertices;j++) pp[j]=ne[j];
	delete [] ne;ne=pp;
#endif
	pnu=new int[i];
	for(j=0;j<currentvertices;j++) pnu[j]=nu[j];
	delete [] nu;nu=pnu;
	ppts=new fpoint[3*i];
	for(j=0;j<3*currentvertices;j++) ppts[j]=pts[j];
	delete [] pts;sure.p=pts=ppts;
	currentvertices=i;
};

/** Doubles the maximum allowed vertex order, by reallocating mem, mep, mec,
 * and mnu <em>(Neighbor version only)</em> arrays. If the allocation exceeds
 * the absolute maximum set in maxvertexorder, then the routine causes a fatal
 * error. */
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
#ifdef FACETS_NEIGHBOR
	p2=new int*[i];
	for(j=0;j<currentvertexorder;j++) p2[j]=mne[j];
	delete [] mne;mne=p2;
#endif
	currentvertexorder=i;
};

/** Doubles the size allocation of the main delete stack. If the allocation
 * exceeds the absolute maximum set in maxdeletesize, then routine causes a
 * fatal error. */
void voronoicell::addmemory_ds() {
	int i=2*currentdeletesize,j,*pds;
	if (i>maxdeletesize) throw fatal_error("Delete stack 1 memory allocation exceeded absolute maximum");
	cerr << "Delete stack 1 memory scaled up to " << i << endl;
	pds=new int[i];
	for(j=0;j<currentdeletesize;j++) pds[j]=ds[j];
	delete [] ds;ds=pds;
	currentdeletesize=i;
};

/** Doubles the size allocation of the auxiliary delete stack. If the
 * allocation exceeds the absolute maximum set in maxdeletesize2, then the
 * routine causes a fatal error. */
void voronoicell::addmemory_ds2() {
	int i=2*currentdeletesize2,j,*pds2;
	if (i>maxdeletesize2) throw fatal_error("Delete stack 2 memory allocation exceeded absolute maximum");
	cerr << "Delete stack 2 memory scaled up to " << i << endl;
	pds2=new int[i];
	for(j=0;j<currentdeletesize2;j++) pds2[j]=ds2[j];
	delete [] ds2;ds2=pds2;
	currentdeletesize2=i;
};

/** Initializes a Voronoi cell as a rectangular box with the given dimensions */
void voronoicell::init(fpoint xmin,fpoint xmax,fpoint ymin,fpoint ymax,fpoint zmin,fpoint zmax) {
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

/** Initializes a Voroni cell as a regular octahedron.
 * \param l The distance from the octahedron center to a vertex. Six vertices
 * are initialized at (-l,0,0), (l,0,0), (0,-l,0), (0,l,0), (0,0,-l), and (0,0,l).*/
inline void voronoicell::init_octahedron(fpoint l) {
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


/** Initializes an arbitrary test object using the add_vertex() and
 * relconstruct() routines. See the source code for information about
 * the specific objects.
 * \param[in] n the number of the test object (from 0 to 9)
 */
inline void voronoicell::init_test(int n) {
	for(int i=0;i<initvertexorder;i++) mec[i]=0;p=0;
	switch(n) {
		case 0:
			// A peaked object, with a high vertex 6, and a ridge
			// at z=0 from vertex 1 to 2. This can be used to test
			// the order 1 vertex collapse routine.
			add_vertex(1,-2,-1,3,1,5);
			add_vertex(0,-1,1,5,0,2);
			add_vertex(0,1,0,1,3,6,4);
			add_vertex(1,4,-1,4,6,2,0);
			add_vertex(-1,4,-1,3,5,2,6);
			add_vertex(-1,-2,-1,4,0,1);
			add_vertex(0,3,0,4,2,3);
			break;
		case 1:
			// A truncated pyramid shape, with vertex 4 in the z=0
			// plane. This can be used to test order 4 vertex
			// generation.
			add_vertex(-2,2,-1,3,4,1);
			add_vertex(2,2,-1,0,5,2);
			add_vertex(2,-2,-1,1,6,3);
			add_vertex(-2,-2,-1,2,7,0);
			add_vertex(-1,1,0,0,7,5);
			add_vertex(1,1,1,1,4,6);
			add_vertex(1,-1,1,7,2,5);
			add_vertex(-1,-1,1,4,3,6);
			break;
		case 2:
			// An object with two peaks at vertices 4 and 6,
			// connected with a trough at vertex 5. It can be used
			// to test the part of the routine that deals with
			// augmenting existing vertices.
			add_vertex(1,-2,-1,1,3,4);
			add_vertex(-1,-2,-1,2,0,4,5);
			add_vertex(-1,2,-1,1,6,3);
			add_vertex(1,2,-1,2,6,5,0);
			add_vertex(0,-1,1,5,1,0);
			add_vertex(0,0,0,1,4,3,6);
			add_vertex(0,1,1,2,5,3);
			break;
		case 3:
			// A box with a pyramid on top of it. This a good
			// test object to make sure that the code can handle
			// object with vertices of different orders.
			add_vertex(-1,-1,-1,4,3,1);
			add_vertex(1,-1,-1,0,2,5);
			add_vertex(1,1,-1,6,1,3);
			add_vertex(-1,1,-1,2,0,7);
			add_vertex(-1,-1,1,7,0,5,8);
			add_vertex(1,-1,1,4,1,6,8);
			add_vertex(1,1,1,5,2,7,8);
			add_vertex(-1,1,1,6,3,4,8);
			add_vertex(0,0,2,7,4,5,6);
			break;
		case 4:
			// A shape with two peaks (vertices 6 and 9) connected
			// with a trough (vertices 7 and 8). It can be used as
			// a basic test of the double edge skipping in the plane
			// generation routine, whereby an edge is omitted if
			// the routine is tracing along a part it has already
			// been down.
			add_vertex(1,-3,-1,5,6,1);
			add_vertex(-1,-3,-1,0,6,2);
			add_vertex(-3,0,-1,1,7,8,3);
			add_vertex(-1,3,-1,2,9,4);
			add_vertex(1,3,-1,3,9,5);
			add_vertex(3,0,-1,4,8,7,0);
			add_vertex(0,-2,1,7,1,0);
			add_vertex(0,-1,0,8,2,6,5);
			add_vertex(0,1,0,9,2,7,5);
			add_vertex(0,2,1,3,8,4);
			break;
		case 5:
			// An object with four peaks (vertices 8 to 11)
			// connected by a trough (vertex 16). It can be used to
			// the test multiple augmentation of edges of a
			// particular vertex.
			add_vertex(-1,-3,-1,7,1,8,12);
			add_vertex(1,-3,-1,2,12,8,0);
			add_vertex(3,-1,-1,3,9,13,1);
			add_vertex(3,1,-1,4,13,9,2);
			add_vertex(1,3,-1,5,10,14,3);
			add_vertex(-1,3,-1,6,14,10,4);
			add_vertex(-3,1,-1,7,11,15,5);
			add_vertex(-3,-1,-1,0,15,11,6);
			add_vertex(0,-2,1,0,1,12);
			add_vertex(2,0,1,2,3,13);
			add_vertex(0,2,1,14,4,5);
			add_vertex(-2,0,1,7,15,6);
			add_vertex(0,-1,0.5,0,8,1,16);
			add_vertex(1,0,0.5,2,9,3,16);
			add_vertex(0,1,0.5,4,10,5,16);
			add_vertex(-1,0,0.5,6,11,7,16);
			add_vertex(0,0,0,15,12,13,14);
			break;
		case 6:
			// An object with four peaks (vertices 8 to 11)
			// connected by a sequence of ridges. It can be used to
			// the test multiple cases of double edge skipping
			// on the same vertex.
			add_vertex(-1,-3,-1,7,1,8,12);
			add_vertex(1,-3,-1,2,12,8,0);
			add_vertex(3,-1,-1,3,9,13,1);
			add_vertex(3,1,-1,4,13,9,2);
			add_vertex(1,3,-1,5,10,14,3);
			add_vertex(-1,3,-1,6,14,10,4);
			add_vertex(-3,1,-1,7,11,15,5);
			add_vertex(-3,-1,-1,0,15,11,6);
			add_vertex(0,-2,1,0,1,12);
			add_vertex(2,0,1,2,3,13);
			add_vertex(0,2,1,14,4,5);
			add_vertex(-2,0,1,7,15,6);
			add_vertex(0,-1,0,0,8,1,16);
			add_vertex(1,0,0,2,9,3,16);
			add_vertex(0,1,0,4,10,5,16);
			add_vertex(-1,0,0,6,11,7,16);
			add_vertex(0,0,0,15,12,13,14);
			break;
		case 7:
			// A variation on the zeroth test shape, with a peak
			// (vertices 7 and 8) connected to a ridge at z=0
			// (vertices 4 and 5)
			add_vertex(2,-3,-1,3,4,1);
			add_vertex(-2,-3,-1,0,4,2);
			add_vertex(-2,3,-1,1,7,3);
			add_vertex(2,3,-1,2,6,0);
			add_vertex(0,-2,0,0,5,1);
			add_vertex(0,1,0,4,6,7);
			add_vertex(1,2,1,7,5,3);
			add_vertex(-1,2,1,5,6,2);
			break;
		case 8:
			// A triangular object with a skewed peak, that can be used to
			// test the order two vertex removal routine
			add_vertex(3,-2,-1,2,3,1);
			add_vertex(-3,-2,-1,0,4,2);
			add_vertex(0,4,-1,1,5,0);
			add_vertex(1.5,-1,0,7,0,6);
			add_vertex(-1.5,-1,0,1,7,8);
			add_vertex(0,2,0,8,6,2);
			add_vertex(0.75,0.5,0,9,3,5);
			add_vertex(0,-1,0,4,3,9);
			add_vertex(-0.75,0.5,0,4,9,5);
			add_vertex(0,0,1,8,7,6);
			break;
		case 9:
			// This a tetrahedron with some low-order extraneous edges, and can
			// be used to test the order 1 and order 2 removal routines
			add_vertex(0,0,0,3,1,2);
			add_vertex(1,0,1,3,2,0);
			add_vertex(1,1,0,3,0,1);
			add_vertex(2,0,0,6,4,2,1,0);
			add_vertex(3,1,0,3,6,8,5);
			add_vertex(3,2,0,4);
			add_vertex(4,0,0,4,3,7,8);
			add_vertex(5,0,0,6);
			add_vertex(4,1,0,6,4);
	}

	relconstruct();
#ifdef FACETS_NEIGHBOR
	label_facets();
#endif
};

/** Adds an order one vertex to the memory structure, and specifies its edge.
 * \param[in] (x,y,z) are the coordinates of the vertex
 * \param[in] a is the first and only edge of this vertex
 */
void voronoicell::add_vertex(fpoint x,fpoint y,fpoint z,int a) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=1;
	if (mem[1]==mec[1]) addmemory(1);
#ifdef FACETS_NEIGHBOR
	ne[p]=mne[1]+mec[1];
#endif
	int *q=mep[1]+3*mec[1]++;ed[p]=q;
	q[0]=a;q[2]=p++;
};

/** Adds an order 2 vertex to the memory structure, and specifies its edges. */
void voronoicell::add_vertex(fpoint x,fpoint y,fpoint z,int a,int b) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=2;
	if (mem[2]==mec[2]) addmemory(2);
#ifdef FACETS_NEIGHBOR
	ne[p]=mne[2]+2*mec[2];
#endif
	int *q=mep[2]+5*mec[2]++;ed[p]=q;
	q[0]=a;q[1]=b;q[4]=p++;
};

/** Adds an order 3 vertex to the memory structure, and specifies its edges. */
void voronoicell::add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=3;
	if (mem[3]==mec[3]) addmemory(3);
#ifdef FACETS_NEIGHBOR
	ne[p]=mne[3]+3*mec[3];
#endif
	int *q=mep[3]+7*mec[3]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[6]=p++;
};

/** Adds an order 4 vertex to the memory structure, and specifies its edges. */
void voronoicell::add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c,int d) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=4;
	if (mem[4]==mec[4]) addmemory(4);
#ifdef FACETS_NEIGHBOR
	ne[p]=mne[4]+4*mec[4];
#endif
	int *q=mep[4]+9*mec[4]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[8]=p++;
};

/** Adds an order 5 vertex to the memory structure, and specifies its edges. */
void voronoicell::add_vertex(fpoint x,fpoint y,fpoint z,int a,int b,int c,int d,int e) {
	pts[3*p]=x;pts[3*p+1]=y;pts[3*p+2]=z;nu[p]=5;
	if (mem[5]==mec[5]) addmemory(5);
#ifdef FACETS_NEIGHBOR
	ne[p]=mne[5]+5*mec[5];
#endif
	int *q=mep[5]+11*mec[5]++;ed[p]=q;
	q[0]=a;q[1]=b;q[2]=c;q[3]=d;q[4]=e;q[10]=p++;
};

/** Checks that the relational table of the Voronoi cell is accurate, and prints
 * out any errors. This algorithm is O(p), so running it every time the plane
 * routine is called will result in a significant slowdown. */
inline void voronoicell::relcheck() {
	int i,j;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if (ed[ed[i][j]][ed[i][nu[i]+j]]!=i) cout << "Relational error at point " << i << ", edge " << j << "." << endl;
		}
	}
};

/** This routine checks for any two vertices that are connected by more than one
 * edge. The plane algorithm is designed so that this should not happen, so any
 * occurrences are most likely errors. Note that the routine is O(p), so
 * running it every time the plane routine is called will result in a significant
 * slowdown. */
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

/** Constructs the relational table if the edges have been specified. */
inline void voronoicell::relconstruct() {
	int i,j,k,l;
	for(i=0;i<p;i++) for(j=0;j<nu[i];j++) {
		k=ed[i][j];
		l=0;
		while(ed[k][l]!=i) {
			l++;
			if (l==nu[k]) throw fatal_error("Relation table construction failed");
		}
		ed[i][nu[i]+j]=l;
	}
};

/** Cuts the Voronoi cell by a particle whose center is at a separation of
 * (x,y,z) from the cell center. The value of rsq should be initially set to
 * \f$x^2+y^2+z^2\f$. */
#ifdef FACETS_NEIGHBOR
bool voronoicell::nplane(fpoint x,fpoint y,fpoint z,fpoint rsq,int p_id) {
	int *nep,*ned;
#else
bool voronoicell::plane(fpoint x,fpoint y,fpoint z,fpoint rsq) {
#endif
	int count=0,i,j,k,up=0,lp=0,tp,cp,qp=1,rp,stack=0;stack2=0;
	int us=0,ls=0,qs,iqs,cs,uw,qw=0,lw,tw;
	int *edp,*edd;
	fpoint u,l,t,r,q;bool complicatedsetup=false,newdoubleedge=false,doubleedge=false;

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
			us=cycle_down(i,up);
			while(i<j) {
				qp=ed[up][i];
				qs=ed[up][nu[up]+i];
#ifdef FACETS_NEIGHBOR
				ne[p][k]=ne[up][i];
#endif
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
#ifdef FACETS_NEIGHBOR
				ne[p][k]=ne[up][i];
#endif
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
#ifdef FACETS_NEIGHBOR
				ne[p][k]=ne[up][i];
#endif
				ed[p][k]=qp;
				ed[p][nu[p]+k]=qs;
				ed[qp][qs]=p;
				ed[qp][nu[qp]+qs]=k;
				ed[up][i]=-1;
				i++;k++;
			}
			qs=j;
		}
#ifdef FACETS_NEIGHBOR
		if (!doubleedge) {
			ne[p][k]=ne[up][qs];
			ne[p][0]=p_id;
		} else ne[p][0]=ne[up][qs];
#endif
		
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
		ne[p][1]=ne[up][us];
		ne[p][2]=ne[lp][ls];
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
		qs=cycle_up(us,up);
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
			qs=cycle_up(ed[qp][nu[qp]+qs],lp);
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
			ne[p][1]=ne[qp][qs];
			ne[p][2]=ne[lp][ls];
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
			qs=cycle_up(qs,qp);
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
				qs=cycle_up(qs,qp);
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
			// the existing one
			if(j>0) {

				// If we're augmenting a vertex but we don't
				// actually need any more edges, just skip this
				// routine to avoid memory confusion
				if(nu[j]!=k) {
					// Allocate memory and copy the edges
					// of the previous instance into it
#ifdef FACETS_NEIGHBOR
					nep=mne[k]+k*mec[k];
#endif
					edp=mep[k]+(2*k+1)*mec[k]++;
					i=0;
					while(i<nu[j]) {
#ifdef FACETS_NEIGHBOR
						nep[i]=ne[j][i];
#endif
						edp[i]=ed[j][i];
						edp[k+i]=ed[j][nu[j]+i];
						i++;
					}
					edp[2*k]=j;

					// Remove the previous instance with
					// fewer vertices from the memory
					// structure
					edd=mep[nu[j]]+(2*nu[j]+1)*--mec[nu[j]];
					if(edd!=ed[j]) {
						for(lw=0;lw<=2*nu[j];lw++) ed[j][lw]=edd[lw];
#ifdef FACETS_NEIGHBOR
						ned=mne[nu[j]]+nu[j]*mec[nu[j]];
						for(lw=0;lw<nu[j];lw++) ne[j][lw]=ned[lw];
						ne[edd[2*nu[j]]]=ne[j];
#endif
						ed[edd[2*nu[j]]]=ed[j];
					}
#ifdef FACETS_NEIGHBOR
					ne[j]=nep;
#endif
					ed[j]=edp;
				} else i=nu[j];
			} else {

				// Allocate a new vertex of order k
#ifdef FACETS_NEIGHBOR
				ne[p]=mne[k]+k*mec[k];
#endif
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
#ifdef FACETS_NEIGHBOR
				ne[j][i]=p_id;
#endif
				ed[cp][cs]=j;
				ed[cp][nu[cp]+cs]=i;
				i++;
			}

			// Copy in the edges of the underlying vertex,
			// and do one less if this was a double edge 
			qs=iqs;
			while(i<(newdoubleedge?k:k-1)) {
				qs=cycle_up(qs,qp);
				lp=ed[qp][qs];ls=ed[qp][nu[qp]+qs];
#ifdef FACETS_NEIGHBOR
				ne[j][i]=ne[qp][qs];
#endif
				ed[j][i]=lp;
				ed[j][nu[j]+i]=ls;
				ed[lp][ls]=j;
				ed[lp][nu[lp]+ls]=i;
				ed[qp][qs]=-1;
				i++;
			}
			qs=cycle_up(qs,qp);
			cs=i;
			cp=j;
#ifdef FACETS_NEIGHBOR
			ne[j][newdoubleedge?0:cs]=ne[qp][qs];
#endif

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
#ifdef FACETS_NEIGHBOR
			for(i=0;i<j;i++) ne[p][i]=(mne[j]+j*mec[j])[i];
			ne[ed[p][2*j]]=ne[p];
#endif
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
#ifdef FACETS_NEIGHBOR
			for(i=0;i<j;i++) ne[qp][i]=(mne[j]+j*mec[j])[i];
			ne[ed[qp][2*j]]=ne[qp];
			ne[qp]=ne[p];
#endif
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

/** During the creation of a new facet in the plane routine, it is possible
 * that some order two vertices may arise. This routine removes them.
 * Suppose an order two vertex joins c and d. If there's a edge between
 * c and d already, then the order two vertex is just removed; otherwise,
 * the order two vertex is removed and c and d are joined together directly.
 * It is possible this process will create order two or order one vertices,
 * and the routine is continually run until all of them are removed.
 * \return false if the vertex removal was unsuccessful, indicative of
 * the cell reducing to zero volume and disappearing; true if the vertex
 * removal was successful. */
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
#ifdef FACETS_NEIGHBOR
			if (!delete_connection(j,a,false)) return false;
			if (!delete_connection(k,b,true)) return false;
#else
			if (!delete_connection(j,a)) return false;
			if (!delete_connection(k,b)) return false;
#endif
		}

		// Compact the memory
		--p;
		if(p!=i) {
			pts[3*i]=pts[3*p];
			pts[3*i+1]=pts[3*p+1];
			pts[3*i+2]=pts[3*p+2];
			for(k=0;k<nu[p];k++) ed[ed[p][k]][ed[p][nu[p]+k]]=i;
#ifdef FACETS_NEIGHBOR
			ne[i]=ne[p];
#endif
			ed[i]=ed[p];
			nu[i]=nu[p];
			ed[i][2*nu[i]]=i;
		}

		// Collapse any order 1 vertices if they were created
		if(!collapseorder1()) return false;
	}
	return true;
};

/** Order one vertices can potentially be created during the order two collapse
 * routine. This routine keeps removing them until there are none left.
 * \return false if the vertex removal was unsuccessful, indicative of
 * the cell having zero volume and disappearing; true if the vertex removal
 * was successful. */
inline bool voronoicell::collapseorder1() {
	int i,j,k;
	while(mec[1]>0) {
		cerr << "Order one collapse" << endl;
		i=--mec[1];
		j=mep[1][3*i];k=mep[1][3*i+1];
		i=mep[1][3*i+2];
#ifdef FACETS_NEIGHBOR
		if(!delete_connection(j,k,false)) return false;
#else
		if(!delete_connection(j,k)) return false;
#endif
		--p;
		if(p!=i) {
			pts[3*i]=pts[3*p];
			pts[3*i+1]=pts[3*p+1];
			pts[3*i+2]=pts[3*p+2];
			for(k=0;k<nu[p];k++) ed[ed[p][k]][ed[p][nu[p]+k]]=i;
#ifdef FACETS_NEIGHBOR
			ne[i]=ne[p];
#endif
			ed[i]=ed[p];
			nu[i]=nu[p];
			ed[i][2*nu[i]]=i;
		}
	}
	return true;
};

/** This routine deletes the kth edge of vertex j and reorganizes the memory.
 * If the neighbor computation is enabled, we also have to supply an
 * handedness flag to decide whether to preserve the plane on the left
 * or right of the connection.
 * \return false if a zero order vertex was formed, indicative of the cell
 * disappearing; true if the vertex removal was successful. */
#ifdef FACETS_NEIGHBOR
inline bool voronoicell::delete_connection(int j,int k,bool hand) {
	int *nep,*ned,q=hand?k:cycle_up(k,j);
#else
inline bool voronoicell::delete_connection(int j,int k) {
#endif
	int i=nu[j]-1,l,*edp,*edd,m;
	if(i<1) {
		cout << "Zero order vertex formed" << endl;
		return false;
	}
	if(mec[i]==mem[i]) addmemory(i);
#ifdef FACETS_NEIGHBOR
	nep=mne[i]+i*mec[i];
	for(l=0;l<q;l++) nep[l]=ne[j][l];
	while(l<i) {
		nep[l]=ne[j][l+1];
		l++;
	}
#endif
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
#ifdef FACETS_NEIGHBOR
	ned=mne[nu[j]]+nu[j]*mec[nu[j]];
	for(l=0;l<nu[j];l++) ne[j][l]=ned[l];
	ne[edd[2*nu[j]]]=ned;
	ne[j]=nep;
#endif
	ed[edd[2*nu[j]]]=edd;
	ed[j]=edp;
	nu[j]=i;
	return true;
};

/** Cuts a Voronoi cell using the influence of a particle at (x,y,z), first
 * calculating the modulus squared of this vector before passing it to the
 * routine above. */
inline bool voronoicell::plane(fpoint x,fpoint y,fpoint z) {
	fpoint rsq=x*x+y*y+z*z;
	return plane(x,y,z,rsq);
};

/** For the neighbor-tracking version of the code, an extra version of the plane
 * routine is provided that doesn't require passing a plane ID. It just makes
 * up the plane ID to be zero. An nplane routine that works without passing
 * the modulus squared is also provided. */
#ifdef FACETS_NEIGHBOR
inline bool voronoicell::plane(fpoint x,fpoint y,fpoint z,fpoint rsq) {
	return nplane(x,y,z,rsq,0);
};
inline bool voronoicell::nplane(fpoint x,fpoint y,fpoint z,int p_id) {
	fpoint rsq=x*x+y*y+z*z;
	return nplane(x,y,z,rsq,p_id);
};	
#endif

/** This is a simple inline function for picking out the index of the next edge
 * counterclockwise at the current vertex.
 * \param[in] a The index of an edge of the current vertex.
 * \param[in] p The number of the vertex.
 * \return 0 if a=nu[p]-1, or a+1 otherwise. */
inline int voronoicell::cycle_up(int a,int p) {
	return a==nu[p]-1?0:a+1;
};

/** This is a simple inline function for picking out the index of the next edge
 * clockwise from the current vertex.
 * \param[in] a The index of an edge of the current vertex.
 * \param[in] p The number of the vertex.
 * \return nu[p]-1 if a=0, or a-1 otherwise. */
inline int voronoicell::cycle_down(int a,int p) {
	return a==0?nu[p]-1:a-1;
};

/** Calculates the volume of the Voronoi cell, by decomposing the cell into
 * tetrahedra extending outward from the zeroth vertex, which are evaluated
 * using a scalar triple product.
 * \return A floating point number holding the calculated volume. */
fpoint voronoicell::volume() {
	const fpoint fe=1/48.0;
	fpoint vol=0;
	int i,j,k,l,m,n;
	fpoint ux,uy,uz,vx,vy,vz,wx,wy,wz;
	for(i=1;i<p;i++) {
		ux=pts[0]-pts[3*i];
		uy=pts[1]-pts[3*i+1];
		uz=pts[2]-pts[3*i+2];
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				ed[i][j]=-1-k;
				l=cycle_up(ed[i][nu[i]+j],k);
				vx=pts[3*k]-pts[0];
				vy=pts[3*k+1]-pts[1];
				vz=pts[3*k+2]-pts[2];
				m=ed[k][l];ed[k][l]=-1-m;
				while(m!=i) {
					n=cycle_up(ed[k][nu[k]+l],m);
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

/** Computes the maximum radius squared of a vertex from the center
 * of the cell. It can be used to determine when enough particles have
 * been testing an all planes that could cut the cell have been considered.
 * \return The maximum radius squared of a vertex.*/
fpoint voronoicell::maxradsq() {
	int i;fpoint r,s;
	r=pts[0]*pts[0]+pts[1]*pts[1]+pts[2]*pts[2];
	for(i=3;i<3*p;i+=3) {
		s=pts[i]*pts[i]+pts[i+1]*pts[i+1]+pts[i+2]*pts[i+2];
		if(s>r) r=s;
	}
	return r;
};

/** Outputs the edges of the Voronoi cell (in POV-Ray format) to an open file
 * stream, displacing the cell by an amount (x,y,z).
 * \param[in] of A output stream to write to.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position. */
void voronoicell::dumppov(ostream &os,fpoint x,fpoint y,fpoint z) {
	int i,j,k;fpoint ux,uy,uz;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];uz=z+0.5*pts[3*i+2];
		os << "sphere{<" << ux << "," << uy << "," << uz << ">,r}" << endl;
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k<i) os << "cylinder{<" << ux << "," << uy << "," << uz << ">,<" << x+0.5*pts[3*k] << "," << y+0.5*pts[3*k+1] << "," << z+0.5*pts[3*k+2] << ">,r}" << endl;
		}
	}
};

/** An overloaded version of the dumppov routine, that outputs the edges of
 * the Voronoi cell (in POV-Ray format) to a file.
 * \param[in] filename The file to write to.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumppov(char *filename,fpoint x,fpoint y,fpoint z) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	dumppov(os,x,y,z);
	os.close();
};

/** An overloaded version of the dumppov routine, that outputs the edges of the
 * Voronoi cell (in POV-Ray format) to standard output. 
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumppov(fpoint x,fpoint y,fpoint z) {
	dumppov(cout,x,y,z);
};

/** Outputs the edges of the Voronoi cell (in gnuplot format) to an output stream.
 * \param[in] of An output stream to write to.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
void voronoicell::dumpgnuplot(ostream &os,fpoint x,fpoint y,fpoint z) {
	int i,j,k;fpoint ux,uy,uz;
	for(i=0;i<p;i++) {
		ux=x+0.5*pts[3*i];uy=y+0.5*pts[3*i+1];uz=z+0.5*pts[3*i+2];
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (ed[i][j]<i) os << ux << " " << uy << " " << uz << endl << x+0.5*pts[3*k] << " " << y+0.5*pts[3*k+1] << " " << z+0.5*pts[3*k+2] << endl << endl << endl;
		}
	}
};

/** An overloaded version of the dumpgnuplot routine that writes directly to
 * a file.
 * \param[in] filename The name of the file to write to. 
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumpgnuplot(char *filename,fpoint x,fpoint y,fpoint z) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	dumpgnuplot(os,x,y,z);
	os.close();
};

/** An overloaded version of the dumpgnuplot routine, that prints to the
 * standard output.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumpgnuplot(fpoint x,fpoint y,fpoint z) {
	dumpgnuplot(cout,x,y,z);
};

/** Outputs the Voronoi cell in the POV mesh2 format, described in section
 * 1.3.2.2 of the POV-Ray documentation. The mesh2 output consists of a list of
 * vertex vectors, followed by a list of triangular faces. The routine also
 * makes use of the optional inside_vector specification, which makes the mesh
 * object solid, so the the POV-Ray Constructive Solid Geometry (CSG) can be
 * applied.
 * \param[in] of An output stream to write to.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumppovmesh(ostream &os,fpoint x,fpoint y,fpoint z) {
	int i,j,k,l,m,n;
	os << "mesh2 {" << endl << "vertex_vectors {" << endl << p << "," << endl;
	for(i=0;i<p;i++) {
		os << "<" << x+0.5*pts[3*i] << "," << y+0.5*pts[3*i+1] << "," << z+0.5*pts[3*i+2] << ">," << endl;
	}
	os << "}" << endl << "face_indices {" << endl << 2*(p-2) << "," << endl;
	for(i=1;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				ed[i][j]=-1-k;
				l=cycle_up(ed[i][nu[i]+j],k);
				m=ed[k][l];ed[k][l]=-1-m;
				while(m!=i) {
					n=cycle_up(ed[k][nu[k]+l],m);
					os << "<" << i << "," << k << "," << m << ">" << endl;
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
	os << "}" << endl << "inside_vector <0,0,1>" << endl << "}" << endl;
};

/** An overloaded version of the dumppovmesh routine, that writes directly to a
 * file.
 * \param[in] filename A filename to write to.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumppovmesh(char *filename,fpoint x,fpoint y,fpoint z) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	dumppovmesh(os,x,y,z);
	os.close();
};

/** An overloaded version of the dumppovmesh routine, that prints to the
 * standard output.
 * \param[in] (x,y,z) A displacement vector to be added to the cell's position.
 */
inline void voronoicell::dumppovmesh(fpoint x,fpoint y,fpoint z) {
	dumppovmesh(cout,x,y,z);
};

/** Randomly perturbs the points in the Voronoi cell by an amount r. */
inline void voronoicell::perturb(fpoint r) {
	for(int i=0;i<3*p;i++) {
		pts[i]+=(2*double(rand())/RAND_MAX-1)*r;
	}
};

/** Initialises the suretest class and creates a buffer for dubious points. */
suretest::suretest() : currentdubious(initdubious) {
	sn=new int[2*currentdubious];
};

/** Suretest destructor to free memory allocation. */
suretest::~suretest() {
	delete [] sn;
};

/** Sets up the suretest class with a particular test plane, and removes
 * any special cases from the table. */
inline void suretest::init(fpoint x,fpoint y,fpoint z,fpoint rsq) {
	sc=0;px=x;py=y;pz=z;prsq=rsq;
};

/** */
inline int suretest::test(int n,fpoint &ans) {
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

/** Prints the vertices, their edges, the relation table, and also notifies if
 * any glaring memory errors are visible. */
void voronoicell::edgeprint() {
	int j;
	for(int i=0;i<p;i++) {
		cout << i << " " << nu[i] << "  ";
		for(j=0;j<nu[i];j++) cout << " " << ed[i][j];
		cout << "    ";
		while(j<2*nu[i]) cout << " " << ed[i][j++];
		cout << "     " << ed[i][j];
#ifdef FACETS_NEIGHBOR
		cout << "    (";
		for(j=0;j<nu[i];j++) {
			cout << ne[i][j] << (j==nu[i]-1?")    ":",");

		}
#endif
		cout << " " << pts[3*i] << " " << pts[3*i+1] << " " << pts[3*i+2];
		cout << " " << ed[i];
		if (ed[i]>=mep[nu[i]]+mec[nu[i]]*(2*nu[i]+1)) cout << " Memory error";
		cout << endl;
	}
};

/** Prints out a list of all the facets and their vertices. If the neighbor option
 * is defined, it lists each cutting plane. */
void voronoicell::facets(ostream &os) {
	int i,j,k,l,m;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
#ifdef FACETS_NEIGHBOR
				os << "(" << i << "," << ne[i][j] << ")";
#else
				os << i;
#endif
				ed[i][j]=-1-k;
				l=cycle_up(ed[i][nu[i]+j],k);
				do {
#ifdef FACETS_NEIGHBOR
					os << " (" << k << "," << ne[k][l] << ")";
#else
					os << " " << k;
#endif
					m=ed[k][l];
					ed[k][l]=-1-m;
					l=cycle_up(ed[k][nu[k]+l],m);
					k=m;
				} while (k!=i);
				os << endl;
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Facet evaluation routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
};

/** <em>(Neighbor version only.)</em> This routine labels the facets
 * in an arbitrary order, starting from one. */
#ifdef FACETS_NEIGHBOR
void voronoicell::label_facets() {
	int i,j,k,l,m,q=1;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				ed[i][j]=-1-k;
				ne[i][j]=q;
				l=cycle_up(ed[i][nu[i]+j],k);
				do {
					m=ed[k][l];
					ed[k][l]=-1-m;
					ne[k][l]=q;
					l=cycle_up(ed[k][nu[k]+l],m);
					k=m;
				} while (k!=i);
				q++;
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Facet labeling routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
};
#endif

/** <em>(Neighbor version only.)</em> This routine checks to make sure the
 * neighbor information of each facets is consistent.*/
#ifdef FACETS_NEIGHBOR
void voronoicell::facet_check() {
	int i,j,k,l,m,q;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				ed[i][j]=-1-k;
				q=ne[i][j];
				l=cycle_up(ed[i][nu[i]+j],k);
				do {
					m=ed[k][l];
					ed[k][l]=-1-m;
					if (ne[k][l]!=q) cerr << "Facet error at (" << k << "," << l << ")=" << ne[k][l] << ", started from (" << i << "," << j << ")=" << q << endl;
					l=cycle_up(ed[k][nu[k]+l],m);
					k=m;
				} while (k!=i);
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Facet labeling routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
};
#endif

/** For the neighbor-tracking version of the code, this routine just provides a list
 * of plane IDs. */
#ifdef FACETS_NEIGHBOR
void voronoicell::neighbors(ostream &os) {
	int i,j,k,l,m;
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			k=ed[i][j];
			if (k>=0) {
				os << " " << ne[i][j];
				ed[i][j]=-1-k;
				l=cycle_up(ed[i][nu[i]+j],k);
				do {
					m=ed[k][l];
					ed[k][l]=-1-m;
					l=cycle_up(ed[k][nu[k]+l],m);
					k=m;
				} while (k!=i);
			}
		}
	}
	for(i=0;i<p;i++) {
		for(j=0;j<nu[i];j++) {
			if(ed[i][j]>=0) throw fatal_error("Neighbor routine didn't look everywhere");
			ed[i][j]=-1-ed[i][j];
		}
	}
}
#endif

/** An overloaded version of facets() which output the results to the standard
 * output. */
inline void voronoicell::facets() {
	facets(cout);
};

/** An overloaded version of facets(), which outputs the results to <filename>. */
inline void voronoicell::facets(char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	facets(os);
	os.close();
};

/** Examines all the facets, and evaluates them by the number of vertices that
 * they have.*/
void voronoicell::facet_statistics(ostream &os) {
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
				l=cycle_up(ed[i][nu[i]+j],k);
				do {
					q++;
					m=ed[k][l];
					ed[k][l]=-1-m;
					l=cycle_up(ed[k][nu[k]+l],m);
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
	for(i=0;i<=maxf;i++) os << i << " " << stat[i] << endl; 
	delete [] stat;
};

/** An overloaded version of facet_statistics() which outputs the results to
 * standard output. */
inline void voronoicell::facet_statistics() {
	facet_statistics(cout);
};

/** An overloaded version of facet_statistics() which outputs the results to
 * <filename>. */
inline void voronoicell::facet_statistics(char *filename) {
	ofstream os;
	os.open(filename,ofstream::out|ofstream::trunc);
	facet_statistics(os);
	os.close();
};
