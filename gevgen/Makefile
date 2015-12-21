SHELL = /bin/bash
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gEvGen_mix

all: $(TGT)

# gevgen event generation driver
#
gEvGen_mix: FORCE
	$(CXX) $(CXXFLAGS) -c $@.cc $(INCLUDES)
	$(LD) $(LDFLAGS) $@.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_mix 

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gevgen_mix

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_mix

FORCE:

# DO NOT DELETE
