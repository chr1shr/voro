# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 28th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=loops polygons odd_even

# Makefile rules
all: $(EXECUTABLES)

loops: loops.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o loops loops.cc -lvoro++

polygons: polygons.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o polygons polygons.cc -lvoro++

odd_even: odd_even.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o odd_even odd_even.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
