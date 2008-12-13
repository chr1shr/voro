#Voro++ makefile
#
#Author : Chris H. Rycroft (LBL / UC Berkeley)
#Email  : chr@alum.mit.edu
#Date   : February 27th 2008

#Load the common configuration file
include ../config.mk

#List of the common source files
SOURCE=container.cc container.hh config.hh cell.hh cell.cc wall.cc wall.hh worklist.cc worklist.hh

#Lists of the command line utilities
EXECUTABLES="../bin/voro++"

#Makefile rules
all: $(EXECUTABLES)

"../bin/voro++": $(SOURCE) cmd_line.cc
	$(CC) $(CFLAGS) -o ../bin/voro++ cmd_line.cc

clean:
	rm -f $(EXECUTABLES)

doc: Doxyfile $(SOURCE)
	doxygen Doxyfile
	cd ../latex && make
