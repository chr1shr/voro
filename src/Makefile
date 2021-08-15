# Voro++, a cell-based Voronoi library
# By Chris H. Rycroft and the Rycroft Group

# Load the common configuration file
include ../config.mk

# List of the common source files
objs=cell_3d.o common.o container_3d.o container_tri.o unitcell.o v_compute_3d.o \
     c_loops.o v_base_3d.o wall.o wall_3d.o particle_list.o 
src=$(patsubst %.o,%.cc,$(objs))

# Makefile rules
all: libvoro++.a voro++

depend:
	$(CXX) -MM $(src) >Makefile.dep

include Makefile.dep

libvoro++.a: $(objs)
	rm -f libvoro++.a
	ar rs libvoro++.a $^

voro++: libvoro++.a cmd_line.cc
	$(CXX) $(CFLAGS) -L. -o voro++ cmd_line.cc -lvoro++

%.o: %.cc
	$(CXX) $(CFLAGS) -c $<

help: Doxyfile $(SOURCE)
	doxygen Doxyfile

clean:
	rm -f $(objs) voro++ libvoro++.a

.PHONY: all help execs depend
