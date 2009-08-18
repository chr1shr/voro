# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Load the common configuration file
include ../../config.mk

# The relative of path of the main library source files
INCLUDE=-I../../src

# List of the common source files
SOURCE=../../src/container.cc ../../src/container.hh ../../src/config.hh \
       ../../src/cell.hh ../../src/cell.cc ../../src/wall.cc ../../src/wall.hh \
       ../../src/worklist.cc ../../src/worklist.hh ../../src/voro++.cc \
       ../../src/voro++.hh

# List of executables
EXECUTABLES=single_cell platonic random_points import 

# Makefile rules
all: single_cell platonic random_points import 

single_cell: $(SOURCE) single_cell.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o single_cell single_cell.cc

platonic: $(SOURCE) platonic.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o platonic platonic.cc

random_points: $(SOURCE) random_points.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o random_points random_points.cc

import: $(SOURCE) import.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o import import.cc

clean:
	rm -f $(EXECUTABLES)
