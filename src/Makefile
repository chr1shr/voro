# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : February 27th 2008

# Load the common configuration file
include ../config.mk

# List of the common source files
SRC_CELL=config.hh cell_2d.cc cell_2d.hh

SRC_CTR=container_2d.cc container_2d.hh

# Lists of the command line utilities
EXECUTABLES=test_cell test_ctr

# Makefile rules
all: $(EXECUTABLES)

test_cell: $(SRC_CELL) test_cell.cc
	$(CC) $(CFLAGS) -o test_cell test_cell.cc 

test_ctr: $(SRC_CELL) $(SRC_CTR) test_ctr.cc
	$(CC) $(CFLAGS) -o test_ctr test_ctr.cc 

#help: Doxyfile $(SRC)
#	doxygen Doxyfile
#	cd ../latex && make

clean:
	rm -f $(EXECUTABLES)

.PHONY: all help clean
