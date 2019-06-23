#
# Makefile to compile the HyyLimitPlotter fitter
#

SHELL = /bin/sh

# ROOT stuff
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTLDFLAGS   = $(shell $(ROOTSYS)/bin/root-config --ldflags)

CXXFLAGS     += $(ROOTCFLAGS) -lMinuit -lMinuit2 -lRooFitCore -lRooFit -lRooStats

main/HiggsCombiner: main/HiggsCombiner.C
	$(CXX) utils/AtlasStyle.C utils/AtlasLabels.C utils/AtlasUtils.C utils/Utils.C \
	main/HiggsDifferentialCombiner.C main/HiggsDifferentialInput.C main/HiggsDifferentialWorkspace.C \
	main/HiggsDifferentialPlot.C \
	$(CXXFLAGS) $(ROOTLIBS) $@.C -o $@

# Clean up: remove executables and outdated files.
clean:
	rm -rf *.o
	rm -rf main/HiggsCombiner
	
