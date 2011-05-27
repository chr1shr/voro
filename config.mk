# Voro++, a 3D cell-based Voronoi library
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# This a common configuration file that includes definitions used by all
# the Makefiles.

# Compiler location and compiler switches. Add "-DFACETS_SINGLE_PRECISION" to
# switch the code from using double floating point arithmetic to just single.
# On a slower computer, it may be worth lowering the level of optimization to
# "-O2" or switching it off completely, to speed up the compilation. See
# config.hh for more details.
CC=g++
CFLAGS=-Wall -ansi -pedantic -O3

E_INC=-I../../src
E_LIB=-L../../src
