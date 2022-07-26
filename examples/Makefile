# Voro++ makefile
# By Chris H. Rycroft and the Rycroft Group

# Makefile rules
all: ex_basic ex_walls ex_custom ex_extra ex_degenerate ex_interface ex_timing

ex_basic:
	$(MAKE) -C basic

ex_walls:
	$(MAKE) -C walls

ex_custom:
	$(MAKE) -C custom

ex_extra:
	$(MAKE) -C extra

ex_degenerate:
	$(MAKE) -C degenerate

ex_interface:
	$(MAKE) -C interface

ex_timing:
	$(MAKE) -C timing

clean:
	$(MAKE) -C basic clean
	$(MAKE) -C walls clean
	$(MAKE) -C custom clean
	$(MAKE) -C extra clean
	$(MAKE) -C degenerate clean
	$(MAKE) -C interface clean
	$(MAKE) -C timing clean

.PHONY: all ex_basic ex_walls ex_custom ex_extra ex_degenerate clean
