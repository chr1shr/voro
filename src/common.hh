// Voro++, a cell-based Voronoi library
// By Chris H. Rycroft and the Rycroft Group

#ifndef VOROPP_COMMON_HH
#define VOROPP_COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "config.hh"

namespace voro {

void check_duplicate(int n,double x,double y,double z,int id,double *qp);

void voro_fatal_error(const char *p,int status);
void voro_print_positions(std::vector<double> &v,FILE *fp=stdout);
FILE* safe_fopen(const char *filename,const char *mode);
void voro_print_vector(std::vector<int> &v,FILE *fp=stdout);
void voro_print_vector(std::vector<double> &v,FILE *fp=stdout);
void voro_print_face_vertices(std::vector<int> &v,FILE *fp=stdout);

}

#endif
