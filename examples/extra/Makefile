# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=intersect lloyd

# Makefile rules
all: $(EXECUTABLES) 

intersect: $(SOURCE) intersect.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o intersect intersect.cc -lvoro++_2d

lloyd: lloyd.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o lloyd lloyd.cc -lvoro++_2d

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
