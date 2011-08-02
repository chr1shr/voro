# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Load the common configuration file
include ../config.mk

# List of executables
EXECUTABLES=network

# Makefile rules
all: $(EXECUTABLES) 

network: network.cc v_network.cc v_network.hh r_table.cc
	$(CC) $(CFLAGS) -I../src -L../src -o v_network v_network.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
