# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Tell make that these are phony targets
.PHONY: all help clean

# Build all of the executable files
all:
	cd src && $(MAKE)
	cd examples && $(MAKE)

# Build the help files (with doxygen)
help:
	cd src && $(MAKE) help

# Clean up the executable files
clean:
	cd examples && $(MAKE) clean
	cd src && $(MAKE) clean
