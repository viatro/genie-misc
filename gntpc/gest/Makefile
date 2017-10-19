SHELL = /bin/bash
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs) -lEG
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = ghep2gest

all: $(TGT)

ghep2gest: FORCE
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $(LIBRARIES) $@.cc -o ./$@

clean: FORCE
	$(RM) *.o ./ghep2gest

FORCE:

# DO NOT DELETE
