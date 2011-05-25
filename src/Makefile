# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : February 27th 2008

# Load the common configuration file
include ../config.mk

# List of the common source files
objs=cell.o container.o unitcell.o v_compute.o v_loops.o v_base.o wall.o
src=$(patsubst %.o,%.cc,$(objs))

# Lists of the command line utilities
execs=single_cell

# Makefile rules
all: libvoro++.a

depend:
	$(CXX) -MM $(src) >Makefile.dep

execs: $(execs)

include Makefile.dep

libvoro++.a: $(objs)
	ar rs libvoro++.a $^

%.o: %.cc
	$(CXX) $(CFLAGS) -c $<

single_cell: single_cell.cc libvoro++.a
	g++-mp-4.5 $(CFLAGS) -funit-at-a-time -L. -o single_cell single_cell.cc -lvoro++

help: Doxyfile $(SOURCE)
	doxygen Doxyfile
	cd ../latex && make

clean:
	rm -f $(objs)

.PHONY: all help execs depend
