# Voro++, a 3D cell-based Voronoi library
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# This a common configuration file that includes definitions used by all
# the Makefiles.

# C++ compiler
CC=g++

# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic -ggdb

# Relative include and library paths for compilation of the examples
E_INC=-I../../src
E_LIB=-L../../src
