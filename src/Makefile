# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : February 27th 2008

# Load the common configuration file
include ../config.mk

# List of the common source files
objs=cell.o common.o container.o unitcell.o v_compute.o v_loops.o v_base.o wall.o pre_container.o
src=$(patsubst %.o,%.cc,$(objs))

# Lists of the command line utilities
execs=single_cell

# Makefile rules
all: libvoro++.a ../bin/voro++

depend:
	$(CXX) -MM $(src) >Makefile.dep

execs: $(execs)

include Makefile.dep

libvoro++.a: $(objs)
	ar rs libvoro++.a $^

../bin/voro++: libvoro++.a cmd_line.cc
	$(CXX) $(CFLAGS) -L. -o ../bin/voro++ cmd_line.cc -lvoro++


%.o: %.cc
	$(CXX) $(CFLAGS) -c $<

help: Doxyfile $(SOURCE)
	doxygen Doxyfile
	cd ../latex && $(MAKE)

clean:
	rm -f $(objs)

.PHONY: all help execs depend
