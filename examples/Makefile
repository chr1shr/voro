# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : May 18th 2011 

# Makefile rules
all: ex_basic ex_walls ex_custom ex_extra ex_degenerate

ex_basic:
	cd basic && $(MAKE) $(AM_MAKEFLAGS)

ex_walls:
	cd walls && $(MAKE) $(AM_MAKEFLAGS)

ex_custom:
	cd custom && $(MAKE) $(AM_MAKEFLAGS)

ex_extra:
	cd extra && $(MAKE) $(AM_MAKEFLAGS)

ex_degenerate:
	cd degenerate && $(MAKE) $(AM_MAKEFLAGS)

clean:
	cd basic && $(MAKE) $(AM_MAKEFLAGS) clean
	cd walls && $(MAKE) $(AM_MAKEFLAGS) clean
	cd custom && $(MAKE) $(AM_MAKEFLAGS) clean
	cd extra && $(MAKE) $(AM_MAKEFLAGS) clean
	cd degenerate && $(MAKE) $(AM_MAKEFLAGS) clean

.PHONY: all ex_basic ex_walls ex_custom ex_extra ex_degenerate clean
