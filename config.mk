# Voro++, a 3D cell-based Voronoi library
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# This a common configuration file that includes definitions used by all
# the Makefiles.

# C++ compiler
CXX=g++-mp-4.8

# C++ compiler flags
CFLAGS=-Wall -ansi -pedantic -O3

# Include and library mat
E_INC=-I../../src
E_LIB=-L../../src
