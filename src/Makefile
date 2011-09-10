# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : February 27th 2008

# Load the common configuration file
include ../config.mk

# List of the common source files
objs=common.o cell_2d.o container_2d.o
src=$(patsubst %.o,%.cc,$(objs))

# Makefile rules
all: libvoro++_2d.a

depend:
	$(CXX) -MM $(src) >Makefile.dep

execs: $(execs)

include Makefile.dep

libvoro++_2d.a: $(objs)
	rm -f libvoro++_2d.a
	ar rs libvoro++_2d.a $^

%.o: %.cc
	$(CXX) $(CFLAGS) -c $<

help: Doxyfile $(SOURCE)
	doxygen Doxyfile
	cd ../latex && $(MAKE)

clean:
	rm -f $(objs) libvoro++_2d.a

.PHONY: all help execs depend
