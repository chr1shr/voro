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
EXECUTABLES=radical_test

# Makefile rules
all: radical_test

radical_test: $(SOURCE) radical_test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o radical_test radical_test.cc

clean:
	rm -f $(EXECUTABLES)
