# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Makefile rules
all:
	cd examples && make
	cd src && make

doc:
	cd src && make doc

clean:
	cd examples && make clean
	cd src && make clean
