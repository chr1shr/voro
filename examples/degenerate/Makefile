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
EXECUTABLES=degenerate degenerate2

# Makefile rules
all: degenerate degenerate2

degenerate: $(SOURCE) degenerate.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o degenerate degenerate.cc

degenerate2: $(SOURCE) degenerate2.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o degenerate2 degenerate2.cc

clean:
	rm $(EXECUTABLES)
