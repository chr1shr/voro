# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=rad_test finite_sys cylinder_inv single_cell_2d period

# Makefile rules
all: $(EXECUTABLES)

rad_test: rad_test.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o rad_test rad_test.cc -lvoro++

finite_sys: finite_sys.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o finite_sys finite_sys.cc -lvoro++

cylinder_inv: cylinder_inv.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o cylinder_inv cylinder_inv.cc -lvoro++

single_cell_2d: single_cell_2d.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o single_cell_2d single_cell_2d.cc -lvoro++

period: period.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o period period.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
