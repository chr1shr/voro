# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=cell_statistics custom_output radical boundaries

# Makefile rules
all: $(EXECUTABLES)

cell_statistics: cell_statistics.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o cell_statistics cell_statistics.cc -lvoro++

custom_output: custom_output.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o custom_output custom_output.cc -lvoro++

radical: radical.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o radical radical.cc -lvoro++

boundaries: boundaries.cc
	$(CC) $(CFLAGS) $(E_INC) $(E_LIB) -o boundaries boundaries.cc

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
