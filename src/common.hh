// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : May 18th 2011

/** \file cell.hh
 * \brief Header file for the small helper functions. */

#ifndef VOROPP_COMMON_HH
#define VOROPP_COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

#include "config.hh"

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
inline void voropp_fatal_error(const char *p,int status) {
	fprintf(stderr,"voro++: %s\n",p);
	exit(status);
}

/** \brief Prints a vector of integers.
 *
 * Prints a vector of integers.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
inline void voropp_print_vector(vector<int> &v,FILE *fp=stdout) {
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
inline void voropp_print_vector(vector<double> &v,FILE *fp=stdout) {
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

/** \brief Prints a vector of positions.
 *
 * Prints a vector of positions as bracketed triplets.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
inline void voropp_print_positions(vector<double> &v,FILE *fp=stdout) {
	if(v.size()>0) {
		fprintf(fp,"(%g,%g,%g)",v[0],v[1],v[2]);
		for(int k=3;(unsigned int) k<v.size();k+=3) {
			fprintf(fp," (%g,%g,%g)",v[k],v[k+1],v[k+2]);
		}
	}
}

/** \brief Prints a vector a face vertex information.
 *
 * Prints a vector of face vertex information. A value is read, which
 * corresponds to the number of vertices in the next face. The routine reads
 * this number of values and prints them as a bracked list. This is repeated
 * until the end of the vector is reached.
 * \param[in] v the vector to interpret and print.
 * \param[in] fp the file stream to print to. */
inline void voropp_print_face_vertices(vector<int> &v,FILE *fp=stdout) {
	int j,k=0,l;
	if(v.size()>0) {
		l=v[k++];
		if(l<=1) {
			if(l==1) fprintf(fp,"(%d)",v[k++]);
			else fputs("()",fp);
		} else {
			j=k+l;
			fprintf(fp,"(%d",v[k++]);
			while(k<j) fprintf(fp,",%d",v[k++]);
			fputs(")",fp);
		}
		while((unsigned int) k<v.size()) {
			l=v[k++];
			if(l<=1) {
				if(l==1) fprintf(fp," (%d)",v[k++]);
				else fputs(" ()",fp);
			} else {
				j=k+l;
				fprintf(fp," (%d",v[k++]);
				while(k<j) fprintf(fp,",%d",v[k++]);
				fputs(")",fp);
			}
		}
	}
}

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
inline FILE* voropp_safe_fopen(const char *filename,const char *mode) {
	FILE *fp(fopen(filename,mode));
	if(fp==NULL) {
		fprintf(stderr,"voro++: Unable to open file '%s'\n",filename);
		exit(VOROPP_FILE_ERROR);
	}
	return fp;
}

#endif
