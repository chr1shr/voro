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
EXECUTABLES=network p_test

# Makefile rules
all: p_test

network: $(SOURCE) network.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o network network.cc

p_test: $(SOURCE) p_test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o p_test p_test.cc

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
