# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : February 27th 2008

# Load the common configuration file
include ../config.mk

# List of the common source files
SOURCE=config.hh cell_2d.cc cell_2d.hh

# Lists of the command line utilities
EXECUTABLES=test_2d

# Makefile rules
all: $(EXECUTABLES)

test_2d: $(SOURCE) test_2d.cc
	$(CC) $(CFLAGS) -o test_2d test_2d.cc 

#help: Doxyfile $(SOURCE)
#	doxygen Doxyfile
#	cd ../latex && make

clean:
	rm -f $(EXECUTABLES)

.PHONY: all help clean
