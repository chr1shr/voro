# Voro++ makefile
#
# Author : Chris H. Rycroft (Harvard University / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../../config.mk

# List of executables
EXECUTABLES=timing_test

# Makefile rules
all: $(EXECUTABLES)

timing_test: timing_test.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o timing_test timing_test.cc -lvoro++

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
