# Voro++ makefile
# By Chris H. Rycroft and the Rycroft Group

# Tell GNU Make that these are phony targets
.PHONY: all help clean install uninstall

include config.mk

# Build all of the executable files
all:
	$(MAKE) -C src
	$(MAKE) -C examples

# Build the help files (with Doxygen)
help:
	$(MAKE) -C src help

# Clean up the executable files
clean:
	$(MAKE) -C src clean
	$(MAKE) -C examples clean

# Install the executable, man page, and shared library
install:
	$(MAKE) -C src
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/bin
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/lib
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/man
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/man/man1
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/include
	$(INSTALL) -d $(IFLAGS_EXEC) $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) man/voro++.1 $(PREFIX)/man/man1
	$(INSTALL) $(IFLAGS) src/libvoro++.a $(PREFIX)/lib
	$(INSTALL) $(IFLAGS_EXEC) src/voro++ $(PREFIX)/bin
	$(INSTALL) $(IFLAGS) src/c_info.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/cell_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/cell_3d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/common.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/config.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/container_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/container_3d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/container_tri.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/iter_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/iter_3d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/particle_list.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/particle_order.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/rad_option.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/unitcell.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_base_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_base_3d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_compute_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/v_compute_3d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/voro++.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/wall.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/wall_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/wall_3d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/worklist_2d.hh $(PREFIX)/include/voro++
	$(INSTALL) $(IFLAGS) src/worklist_3d.hh $(PREFIX)/include/voro++

# Uninstall the executable, man page, and shared library
uninstall:
	rm -f $(PREFIX)/bin/voro++
	rm -f $(PREFIX)/man/man1/voro++.1
	rm -f $(PREFIX)/lib/libvoro++.a
	rm -f $(PREFIX)/include/voro++/c_info.hh
	rm -f $(PREFIX)/include/voro++/cell_2d.hh
	rm -f $(PREFIX)/include/voro++/cell_3d.hh
	rm -f $(PREFIX)/include/voro++/common.hh
	rm -f $(PREFIX)/include/voro++/config.hh
	rm -f $(PREFIX)/include/voro++/container_2d.hh
	rm -f $(PREFIX)/include/voro++/container_3d.hh
	rm -f $(PREFIX)/include/voro++/container_tri.hh
	rm -f $(PREFIX)/include/voro++/iter_2d.hh
	rm -f $(PREFIX)/include/voro++/iter_3d.hh
	rm -f $(PREFIX)/include/voro++/particle_list.hh
	rm -f $(PREFIX)/include/voro++/particle_order.hh
	rm -f $(PREFIX)/include/voro++/rad_option.hh
	rm -f $(PREFIX)/include/voro++/unitcell.hh
	rm -f $(PREFIX)/include/voro++/v_base_2d.hh
	rm -f $(PREFIX)/include/voro++/v_base_3d.hh
	rm -f $(PREFIX)/include/voro++/v_compute_2d.hh
	rm -f $(PREFIX)/include/voro++/v_compute_3d.hh
	rm -f $(PREFIX)/include/voro++/voro++.hh
	rm -f $(PREFIX)/include/voro++/wall.hh
	rm -f $(PREFIX)/include/voro++/wall_2d.hh
	rm -f $(PREFIX)/include/voro++/wall_3d.hh
	rm -f $(PREFIX)/include/voro++/worklist_2d.hh
	rm -f $(PREFIX)/include/voro++/worklist_3d.hh
	rmdir $(PREFIX)/include/voro++
