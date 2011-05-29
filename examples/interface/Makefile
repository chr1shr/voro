# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=loops polygons sub_loop_test

# Makefile rules
all: $(EXECUTABLES)

loops: loops.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o loops loops.cc -lvoro++

polygons: polygons.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o polygons polygons.cc -lvoro++

sub_loop_test: sub_loop_test.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o sub_loop_test sub_loop_test.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
