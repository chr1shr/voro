# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 28th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=box_cut cut_region superellipsoid irregular

# Makefile rules
all: $(EXECUTABLES)

box_cut: box_cut.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o box_cut box_cut.cc -lvoro++

cut_region: cut_region.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o cut_region cut_region.cc -lvoro++

superellipsoid: superellipsoid.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o superellipsoid superellipsoid.cc -lvoro++

irregular: irregular.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o irregular irregular.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
