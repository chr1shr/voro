# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011

# Tell make that these are phony targets
.PHONY: all help clean

# Build all of the executable files
all:
	cd examples && $(MAKE)

# Build the help files (with doxygen)
help:
	cd src && $(MAKE) help

# Clean up the executable files
clean:
	cd examples && $(MAKE) clean
