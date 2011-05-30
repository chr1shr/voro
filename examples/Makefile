# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Makefile rules
all: ex_basic ex_walls ex_custom ex_extra ex_degenerate ex_interface

ex_basic:
	cd basic && $(MAKE)

ex_walls:
	cd walls && $(MAKE)

ex_custom:
	cd custom && $(MAKE)

ex_extra:
	cd extra && $(MAKE)

ex_degenerate:
	cd degenerate && $(MAKE)

ex_interface:
	cd interface && $(MAKE)

clean:
	cd basic && $(MAKE) clean
	cd walls && $(MAKE) clean
	cd custom && $(MAKE) clean
	cd extra && $(MAKE) clean
	cd degenerate && $(MAKE) clean
	cd interface && $(MAKE) clean

.PHONY: all ex_basic ex_walls ex_custom ex_extra ex_degenerate clean
