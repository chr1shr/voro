# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=loops

# Makefile rules
all:loops

loops: loops.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o loops loops.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
