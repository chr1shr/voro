# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Load the common configuration file
include ../config.mk

# The relative of path of the main library source files
INCLUDE=-I../src

# List of the common source files
SOURCE=../src/container_2d.cc ../src/container_2d.hh ../src/config.hh \
       ../src/cell_2d.hh ../src/cell_2d.cc 

# List of executables
EXECUTABLES=test_cell test_ctr test_import test_lloyd

# Makefile rules
all: $(EXECUTABLES) 

test_cell: $(SOURCE) test_cell.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o test_cell test_cell.cc

test_ctr: $(SOURCE) test_ctr.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o test_ctr test_ctr.cc

test_import: $(SOURCE) test_import.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o test_import test_import.cc

test_lloyd: $(SOURCE) test_lloyd.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o test_lloyd test_lloyd.cc

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
