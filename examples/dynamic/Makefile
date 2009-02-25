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
EXECUTABLES=dyn_test square_test packing_gen best_gen yeast

# Makefile rules
all: dyn_test square_test packing_gen best_gen yeast

dyn_test: $(SOURCE) dyn_test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o dyn_test dyn_test.cc

square_test: $(SOURCE) square_test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o square_test square_test.cc

packing_gen: $(SOURCE) packing_gen.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o packing_gen packing_gen.cc

best_gen: $(SOURCE) best_gen.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o best_gen best_gen.cc

yeast: $(SOURCE) yeast.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o yeast yeast.cc

clean:
	rm -f $(EXECUTABLES)
