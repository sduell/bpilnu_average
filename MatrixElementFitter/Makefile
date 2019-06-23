#
# Makefile to compile the Project
#

SHELL = /bin/sh

# ROOT stuff
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTLDFLAGS   = $(shell $(ROOTSYS)/bin/root-config --ldflags)

CXXFLAGS     += $(ROOTCFLAGS) -lMinuit -lMinuit2 -lRooFitCore -lRooFit -lRooStats

FitterMain: FitterMain.C
	$(CXX) ../utils/AtlasStyle.C ../utils/AtlasLabels.C ../utils/AtlasUtils.C ../utils/Utils.C \
	FitClass.C PlotClass.C TheoryClass.C \
	$(CXXFLAGS) $(ROOTLIBS) $@.C -o $@

# Clean up: remove executables and outdated files.
clean:
	rm -rf *.o
	rm -rf FitterMain
	
