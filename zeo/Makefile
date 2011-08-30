# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../config.mk

# List of executables
EXECUTABLES=network images

# Makefile rules
all: $(EXECUTABLES)

# List of the common source files
objs=v_network.o
src=$(patsubst %.o,%.cc,$(objs))

%.o: %.cc
	$(CXX) $(CFLAGS) -I../src -c $<

network: network.cc v_network.o v_network.hh r_table.cc
	$(CC) $(CFLAGS) -I../src -L../src -o network network.cc v_network.o -lvoro++

images: images.cc
	$(CC) $(CFLAGS) -I../src -L../src -o images images.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
