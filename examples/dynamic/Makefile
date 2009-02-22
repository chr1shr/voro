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
       ../../src/voro++.hh ../../src/dynamic.hh ../../src/dynamic.cc

# List of executables
EXECUTABLES=dyn_test

# Makefile rules
all: dyn_test square_test

dyn_test: $(SOURCE) dyn_test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o dyn_test dyn_test.cc

square_test: $(SOURCE) square_test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o square_test square_test.cc

clean:
	rm -f $(EXECUTABLES)
