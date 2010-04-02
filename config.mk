#Voro++, a 3D cell-based Voronoi library
#
#Author : Chris H. Rycroft (LBL / UC Berkeley)
#Email  : chr@alum.mit.edu
#Date   : July 1st 2008

#This a common configuration file that includes definitions used by all
#the Makefiles.

#Compiler location and compiler switches. Add "-DFACETS_SINGLE_PRECISION" to
#switch the code from using double floating point arithmetic to just single.
#On a slower computer, it may be worth lowering the level of optimization to
#"-O2" or switching it off completely, to speed up the compilation. See
#config.hh for more details.
CC=g++
CFLAGS=-Wall -ansi -pedantic -O3

#These optional compiler flags for the GNU C++ compiler force more function
#inlining to take place. Since the code makes heavy use of inlined functions,
#this can result in marginally better performance. On Apple systems, including
#the "-fast" option can also create marginal improvements, but this option
#will break on some Linux and Cygwin systems.
#CFLAGS=-Wall -pedantic -O3 --param large-function-growth=1000 --param max-inline-insns-single=2000 -Winline
