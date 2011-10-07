# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=circle polygon

# Makefile rules
all: $(EXECUTABLES) 

circle: $(SOURCE) circle.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o circle circle.cc -lvoro++_2d

polygon: polygon.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o polygon polygon.cc -lvoro++_2d

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
