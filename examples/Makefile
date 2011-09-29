# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Load the common configuration file
include ../config.mk

# List of executables
EXECUTABLES=test_cell test_ctr test_import test_lloyd test_intersect test_circle

# Makefile rules
all: $(EXECUTABLES) 

test_cell: $(SOURCE) test_cell.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o test_cell test_cell.cc -lvoro++_2d

test_ctr: test_ctr.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o test_ctr test_ctr.cc -lvoro++_2d

test_import: test_import.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o test_import test_import.cc -lvoro++_2d

test_lloyd: test_lloyd.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o test_lloyd test_lloyd.cc -lvoro++_2d

test_intersect: test_intersect.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o test_intersect test_intersect.cc -lvoro++_2d

test_circle: test_circle.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o test_circle test_circle.cc -lvoro++_2d

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
