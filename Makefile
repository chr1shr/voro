# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Makefile rules
all:
	cd examples && make
	cd src && make

help:
	cd src && make help

clean:
	cd examples && make clean
	cd src && make clean
