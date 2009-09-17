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
EXECUTABLES=radical cell_statistics custom_output

# Makefile rules
all: radical cell_statistics custom_output

radical: $(SOURCE) radical.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o radical radical.cc

cell_statistics: $(SOURCE) cell_statistics.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o cell_statistics cell_statistics.cc

custom_output: $(SOURCE) custom_output.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o custom_output custom_output.cc

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
