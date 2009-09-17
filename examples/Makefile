# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : July 1st 2008

# Makefile rules
all:
	cd basic && make
	cd walls && make
	cd custom && make
	cd extra && make
	cd degenerate && make

clean:
	cd basic && make clean
	cd walls && make clean
	cd custom && make clean
	cd extra && make clean
	cd degenerate && make clean

.PHONY: all clean
