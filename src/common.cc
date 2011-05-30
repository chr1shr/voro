// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file common.cc
 * \brief Implementations of the small helper functions. */

#include "common.hh"

/** \brief Prints a vector of integers.
 *
 * Prints a vector of integers.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voropp_print_vector(vector<int> &v,FILE *fp) {
	int k(0),s(v.size());
	while(k+4<s) {
		fprintf(fp,"%d %d %d %d ",v[k],v[k+1],v[k+2],v[k+3]);
		k+=4;
	}
	if(k+3<=s) {
		if(k+4==s) fprintf(fp,"%d %d %d %d",v[k],v[k+1],v[k+2],v[k+3]);
		else fprintf(fp,"%d %d %d",v[k],v[k+1],v[k+2]);
	} else {
		if(k+2==s) fprintf(fp,"%d %d",v[k],v[k+1]);
		else fprintf(fp,"%d",v[k]);
	}
}

/** \brief Prints a vector of doubles.
 *
 * Prints a vector of doubles.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voropp_print_vector(vector<double> &v,FILE *fp) {
	int k(0),s(v.size());
	while(k+4<s) {
		fprintf(fp,"%g %g %g %g ",v[k],v[k+1],v[k+2],v[k+3]);
		k+=4;
	}
	if(k+3<=s) {
		if(k+4==s) fprintf(fp,"%g %g %g %g",v[k],v[k+1],v[k+2],v[k+3]);
		else fprintf(fp,"%g %g %g",v[k],v[k+1],v[k+2]);
	} else {
		if(k+2==s) fprintf(fp,"%g %g",v[k],v[k+1]);
		else fprintf(fp,"%g",v[k]);
	}
}
