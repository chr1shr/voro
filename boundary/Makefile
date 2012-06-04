# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Load the common configuration file
include ../config.mk

# List of executables
EXECUTABLES=voro++_bdry

# Makefile rules
all: $(EXECUTABLES) 

voro++_bdry: voro++_bdry.cc
	$(CXX) $(CFLAGS) -I../src -L../src -o voro++_bdry voro++_bdry.cc -lvoro++_2d

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
