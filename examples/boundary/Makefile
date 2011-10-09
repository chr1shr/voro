# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=nonconvex_cell container_bd

# Makefile rules
all: $(EXECUTABLES) 

nonconvex_cell: $(SOURCE) nonconvex_cell.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o nonconvex_cell nonconvex_cell.cc -lvoro++_2d

container_bd: $(SOURCE) container_bd.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o container_bd container_bd.cc -ggdb

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
