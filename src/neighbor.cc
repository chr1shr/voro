#include "neighbor.hh"

/** The destructor for the neighbor_track class deallocates the arrays
 * for neighbor tracking. */
neighbor_track::~neighbor_track() {
}

/** This allocates a single array for neighbor tracking.
 * \param[in] i the vertex order of the array to be extended.
 * \param[in] m the size of the array to be extended. */
inline void neighbor_track::allocate(int i,int m) {
	;
}

/** This increases the size of the ne[] array.
 * \param[in] i the new size of the array. */
inline void neighbor_track::add_memory_vertices(int i) {
	}

/** This increases the size of the maximum allowable vertex order in the
 * neighbor tracking.
 * \param[in] i the new size of the neighbor vertex order array. */
inline void neighbor_track::add_memory_vorder(int i) {

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
	
}

/** This is a basic operation to copy a neighbor into paux1.*/
inline void neighbor_track::copy_aux1(int a,int b) {
	
}

/** This is a basic operation to copy a neighbor into paux1 with a shift. It is
 * used in the delete_connection() routine of the voronoicell class. */
inline void neighbor_track::copy_aux1_shift(int a,int b) {
	
}

/** This routine sets the second auxiliary pointer to a new section of memory,
 * and then copies existing neighbor information into it. */
inline void neighbor_track::set_aux2_copy(int a,int b) {

}

/** This is a basic routine to copy ne[b] into ne[a]. */
inline void neighbor_track::copy_pointer(int a,int b) {
	
}

/** This sets ne[j] to the first auxiliary pointer. */
inline void neighbor_track::set_to_aux1(int j) {
	
}

/** This sets ne[j] to the second auxiliary pointer. */
inline void neighbor_track::set_to_aux2(int j) {
	
}



/** This allocates a new array and sets the auxiliary pointer to it. */
inline void neighbor_track::allocate_aux1(int i) {
	
}

/** This deletes a particular neighbor array and switches the pointer to the
 * auxiliary pointer. */
inline void neighbor_track::switch_to_aux1(int i) {

}

/** This routine copies neighbor information into the auxiliary pointer. */
inline void neighbor_track::copy_to_aux1(int i,int m) {
	
}

/** This sets ne[k] to the auxiliary pointer with an offset. */
inline void neighbor_track::set_to_aux1_offset(int k,int m) {
	
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

tore the results in. */
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
