#Facets makefile
#
#Author : Chris H. Rycroft (LBL / UC Berkeley)
#Email  : chr@alum.mit.edu
#Date   : February 27th 2008

#Compiler location and compiler switches. Add "-DFACETS_SINGLE_PRECISION" to
#switch the code from using double floating point arithmetic to just single.
#See config.hh for more details.
CC=g++
CFLAGS=-Wall -O3

#Optional switches for the radical tessellation code, and the neighbor
#tracking code
RAD=-DFACETS_RADICAL
NEI=-DFACETS_NEIGHBOR

#List of the common source files
SOURCE=container.cc container.hh config.hh cell.hh cell.cc

#Lists of the command line utilities and the demo scripts 
FACETS=facets facets-radical facets-neighbor facets-radical-neighbor
DEMO=platonic higher_test higher_test2 cell_test voronoi_test

#Makefile rules
all: facets-all demo-all doc

facets-all: $(FACETS)

facets: $(SOURCE) facets.cc
	$(CC) $(CFLAGS) -o facets facets.cc

facets-radical: $(SOURCE) facets.cc
	$(CC) $(CFLAGS) $(RAD) -o facets-radical facets.cc

facets-neighbor: $(SOURCE) facets.cc
	$(CC) $(CFLAGS) $(NEI) -o facets-neighbor facets.cc

facets-radical-neighbor: $(SOURCE) facets.cc
	$(CC) $(CFLAGS) $(RAD) $(NEI) -o facets-radical-neighbor facets.cc

demo-all: $(DEMO)

platonic: $(SOURCE) platonic.cc
	$(CC) $(CFLAGS) -o platonic platonic.cc

higher_test: $(SOURCE) higher_test.cc
	$(CC) $(CFLAGS) -o higher_test higher_test.cc

higher_test2: $(SOURCE) higher_test2.cc
	$(CC) $(CFLAGS) -o higher_test2 higher_test2.cc

cell_test: $(SOURCE) cell_test.cc
	$(CC) $(CFLAGS) -o cell_test cell_test.cc

voronoi_test: $(SOURCE) voronoi_test.cc
	$(CC) $(CFLAGS) -o voronoi_test voronoi_test.cc

clean:
	rm $(FACETS) $(DEMO)

doc: Doxyfile $(SOURCE)
	doxygen Doxyfile
	cd latex ; make
