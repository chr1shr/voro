# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=cylinder tetrahedron frustum torus cylinder_inv

# Makefile rules
all: $(EXECUTABLES)

cylinder: cylinder.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o cylinder cylinder.cc -lvoro++

cylinder_inv: cylinder_inv.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o cylinder_inv cylinder_inv.cc -lvoro++

tetrahedron: tetrahedron.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o tetrahedron tetrahedron.cc -lvoro++

frustum: frustum.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o frustum frustum.cc -lvoro++

torus: torus.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o torus torus.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
