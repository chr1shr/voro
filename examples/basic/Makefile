# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=single_cell platonic random_points import 

# Makefile rules
all: single_cell platonic random_points import 

single_cell: single_cell.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o single_cell single_cell.cc -lvoro++

platonic: platonic.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o platonic platonic.cc -lvoro++

random_points: random_points.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o random_points random_points.cc -lvoro++

import: import.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o import import.cc -lvoro++

#nearest: $(SOURCE) nearest.cc
#	$(CC) $(CFLAGS) $(INCLUDE) -o nearest nearest.cc

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
