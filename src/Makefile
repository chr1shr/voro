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

# Makefile rules
all: help

help: Doxyfile $(SRC_CELL) $(SRC_CTR)
	doxygen Doxyfile
	cd ../latex && make

.PHONY: all help
