// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

/** \file common.cc
 * \brief Header file for the small helper functions. */

#ifndef VOROPP_COMMON_HH
#define VOROPP_COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "config.hh"

namespace voro {

void check_duplicate(int n,double x,double y,double z,int id,double *qp);

bool voro_read_precision(FILE *fp,char *&fmp,int &pr);
void voro_fatal_error(const char *p,int status);
void voro_print_positions_2d(std::vector<double> &v,FILE *fp=stdout);
void voro_print_positions_2d(int pr,std::vector<double> &v,FILE *fp=stdout);
void voro_print_positions_3d(std::vector<double> &v,FILE *fp=stdout);
void voro_print_positions_3d(int pr,std::vector<double> &v,FILE *fp=stdout);
FILE* safe_fopen(const char *filename,const char *mode);
void voro_print_vector(std::vector<int> &v,FILE *fp=stdout);
void voro_print_vector(std::vector<double> &v,FILE *fp=stdout);
void voro_print_vector(int pr,std::vector<double> &v,FILE *fp=stdout);
void voro_print_face_vertices(std::vector<int> &v,FILE *fp=stdout);

}

#endif
