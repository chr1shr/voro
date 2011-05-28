# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=cylinder tetrahedron frustum torus

# Makefile rules
all: $(EXECUTABLES)

cylinder: cylinder.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o cylinder cylinder.cc -lvoro++

tetrahedron: tetrahedron.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o tetrahedron tetrahedron.cc -lvoro++

frustum: frustum.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o frustum frustum.cc -lvoro++

torus: torus.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o torus torus.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
