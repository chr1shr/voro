#include "neighbor.hh"

/** This constructs the neighbor_track class, within a current
 * voronoicell_neighbor class. It allocates memory for neighbor storage in a
 * similar way to the voronoicell constructor.
 * \param[in] ivc a pointer to the parent voronoicell_neighbor class. */
neighbor_track::neighbor_track(voronoicell_base<neighbor_track> *vc_) : vc(vc_) {
	int i;
	mne=new int*[vc->current_vertex_order];
	ne=new int*[vc->current_vertices];
	for(i=0;i<3;i++) mne[i]=new int[init_n_vertices*i];
	mne[3]=new int[init_3_vertices*3];
	for(i=4;i<vc->current_vertex_order;i++) mne[i]=new int[init_n_vertices*i];
}

/** The destructor for the neighbor_track class deallocates the arrays
 * for neighbor tracking. */
neighbor_track::~neighbor_track() {
	for(int i=0;i<vc->current_vertex_order;i++) if(vc->mem[i]>0) delete [] mne[i];
	delete [] mne;
	delete [] ne;
}

/** This allocates a single array for neighbor tracking.
 * \param[in] i the vertex order of the array to be extended.
 * \param[in] m the size of the array to be extended. */
inline void neighbor_track::allocate(int i,int m) {
	mne[i]=new int[m*i];
}

/** This increases the size of the ne[] array.
 * \param[in] i the new size of the array. */
inline void neighbor_track::add_memory_vertices(int i) {
	int **pp;
	pp=new int*[i];
	for(int j=0;j<vc->current_vertices;j++) pp[j]=ne[j];
	delete [] ne;ne=pp;
}

/** This increases the size of the maximum allowable vertex order in the
 * neighbor tracking.
 * \param[in] i the new size of the neighbor vertex order array. */
inline void neighbor_track::add_memory_vorder(int i) {
	int **p2;
	p2=new int*[i];
	for(int j=0;j<vc->current_vertex_order;j++) p2[j]=mne[j];
	delete [] mne;mne=p2;
}

/** This initializes the neighbor information for a rectangular box and is
 * called during the initialization routine for the voronoicell class. */
inline void neighbor_track::init() {
	int *q;
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
}

/** This initializes the neighbor information for an octahedron. The eight
 * initial faces are assigned ID numbers from -1 to -8. */
inline void neighbor_track::init_octahedron() {
	int *q;
	q=mne[4];
	q[0]=-5;q[1]=-6;q[2]=-7;q[3]=-8;
	q[4]=-1;q[5]=-2;q[6]=-3;q[7]=-4;
	q[8]=-6;q[9]=-5;q[10]=-2;q[11]=-1;
	q[12]=-8;q[13]=-7;q[14]=-4;q[15]=-3;
	q[16]=-5;q[17]=-8;q[18]=-3;q[19]=-2;
	q[20]=-7;q[21]=-6;q[22]=-1;q[23]=-4;
	ne[0]=q;ne[1]=q+4;ne[2]=q+8;ne[3]=q+12;ne[4]=q+16;ne[5]=q+20;
}

/** This initializes the neighbor information for an tetrahedron. The four
 * initial faces are assigned ID numbers from -1 to -4.*/
inline void neighbor_track::init_tetrahedron() {
	int *q;
	q=mne[3];
	q[0]=-4;q[1]=-3;q[2]=-2;
	q[3]=-3;q[4]=-4;q[5]=-1;
	q[6]=-4;q[7]=-2;q[8]=-1;
	q[9]=-2;q[10]=-3;q[11]=-1;
	ne[0]=q;ne[1]=q+3;ne[2]=q+6;ne[3]=q+9;
}

/** This is a basic operation to set a new pointer in the ne[] array.
 * \param[in] p the index in the ne[] array to set.
 * \param[in] n the order of the vertex. */
inline void neighbor_track::set_pointer(int p,int n) {
	ne[p]=mne[n]+n*vc->mec[n];
}

/** This is a basic operation to copy ne[c][d] to ne[a][b]. */
inline void neighbor_track::copy(int a,int b,int c,int d) {
	ne[a][b]=ne[c][d];
}

/** This is a basic operation to carry out ne[a][b]=c. */
inline void neighbor_track::set(int a,int b,int c) {
	ne[a][b]=c;
}

/** This is a basic operation to set the auxiliary pointer paux1.
 * \param[in] k the order of the vertex to point to. */
inline void neighbor_track::set_aux1(int k) {
	paux1=mne[k]+k*vc->mec[k];
}

/** This is a basic operation to copy a neighbor into paux1.*/
inline void neighbor_track::copy_aux1(int a,int b) {
	paux1[b]=ne[a][b];
}

/** This is a basic operation to copy a neighbor into paux1 with a shift. It is
 * used in the delete_connection() routine of the voronoicell class. */
inline void neighbor_track::copy_aux1_shift(int a,int b) {
	paux1[b]=ne[a][b+1];
}

/** This routine sets the second auxiliary pointer to a new section of memory,
 * and then copies existing neighbor information into it. */
inline void neighbor_track::set_aux2_copy(int a,int b) {
	paux2=mne[b]+b*vc->mec[b];
	for(int i=0;i<b;i++) ne[a][i]=paux2[i];
}

/** This is a basic routine to copy ne[b] into ne[a]. */
inline void neighbor_track::copy_pointer(int a,int b) {
	ne[a]=ne[b];
}

/** This sets ne[j] to the first auxiliary pointer. */
inline void neighbor_track::set_to_aux1(int j) {
	ne[j]=paux1;
}

/** This sets ne[j] to the second auxiliary pointer. */
inline void neighbor_track::set_to_aux2(int j) {
	ne[j]=paux2;
}

/** This prints out the neighbor information for vertex i. */
inline void neighbor_track::print_edges(int i) {
	int num(vc->nu[i]-1);
	if(num>0) {
		int j;
		printf("     (");
		for(j=0;j<vc->nu[i]-1;j++) {
			printf("%d,",ne[i][j]);
		}
		printf("%d)",ne[i][j]);
	} else printf("     ()");
}

/** This allocates a new array and sets the auxiliary pointer to it. */
inline void neighbor_track::allocate_aux1(int i) {
	paux1=new int[i*vc->mem[i]];
}

/** This deletes a particular neighbor array and switches the pointer to the
 * auxiliary pointer. */
inline void neighbor_track::switch_to_aux1(int i) {
	delete [] mne[i];
	mne[i]=paux1;
}

/** This routine copies neighbor information into the auxiliary pointer. */
inline void neighbor_track::copy_to_aux1(int i,int m) {
	paux1[m]=mne[i][m];
}

/** This sets ne[k] to the auxiliary pointer with an offset. */
inline void neighbor_track::set_to_aux1_offset(int k,int m) {
	ne[k]=paux1+m;
}

/** This routine checks to make sure the neighbor information of each face is
 * consistent.*/
void neighbor_track::check_facets() {
	int **edp(vc->ed),*nup(vc->nu);
	int i,j,k,l,m,q;
	for(i=1;i<vc->p;i++) for(j=0;j<nup[i];j++) {
		k=edp[i][j];
		if(k>=0) {
			edp[i][j]=-1-k;
			q=ne[i][j];
			l=vc->cycle_up(edp[i][nup[i]+j],k);
			do {
				m=edp[k][l];
				edp[k][l]=-1-m;
				fprintf(stderr,"Facet error at (%d,%d)=%d, started from (%d,%d)=%d\n",k,l,ne[k][l],i,j,q);
				l=vc->cycle_up(edp[k][nup[k]+l],m);
				k=m;
			} while (k!=i);
		}
	}
	vc->reset_edges();
}

/** This routine saves a list of plane IDs to a vector.
 * \param[out] v the vector to store the results in. */
void neighbor_track::neighbors(vector<int> &v) {
	int **edp(vc->ed),*nup(vc->nu);
	v.clear();
	int i,j,k,l,m;
	for(i=1;i<vc->p;i++) for(j=0;j<nup[i];j++) {
		k=edp[i][j];
		if(k>=0) {
			v.push_back(ne[i][j]);
			edp[i][j]=-1-k;
			l=vc->cycle_up(edp[i][nup[i]+j],k);
			do {
				m=edp[k][l];
				edp[k][l]=-1-m;
				l=vc->cycle_up(edp[k][nup[k]+l],m);
				k=m;
			} while (k!=i);
		}
	}
	vc->reset_edges();
}
