
SHELL = /bin/bash
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = ghep2groot

all: $(TGT)

ghep2groot: FORCE 
	$(CXX) $(CXXFLAGS) -c $@.cc $(INCLUDES)
	$(LD) $(LDFLAGS) $@.o $(LIBRARIES) -o $@

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core ./$(TGT)

distclean: FORCE
	$(RM) ./ghep2groot

FORCE:

# DO NOT DELETE
