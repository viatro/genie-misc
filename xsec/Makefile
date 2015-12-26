SHELL = /bin/bash
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gXSecCalc

all: $(TGT)

gXSecCalc: FORCE 
	$(CXX) $(CXXFLAGS) -c $@.cc $(INCLUDES)
	$(LD) $(LDFLAGS) $@.o $(LIBRARIES) -o ./gxscalc

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core ./gxscalc

distclean: FORCE
	$(RM) ./gxscalc

FORCE:

# DO NOT DELETE
