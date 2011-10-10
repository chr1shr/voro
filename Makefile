# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011

# Tell make that these are phony targets
.PHONY: all help clean

# Build all of the executable files
all:
	$(MAKE) -C src
	$(MAKE) -C boundary
	$(MAKE) -C examples

# Build the help files (with doxygen)
help:
	$(MAKE) -C src help

# Clean up the executable files
clean:
	$(MAKE) -C src clean
	$(MAKE) -C boundary clean
	$(MAKE) -C examples clean
