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

calc_covmatrix: calc_covmatrix.C
	$(CXX) ../utils/AtlasStyle.C ../utils/AtlasLabels.C ../utils/AtlasUtils.C ../utils/Utils.C \
	$(CXXFLAGS) $(ROOTLIBS) $@.C -o $@

calc_tot_matrices: calc_tot_matrices.C
	$(CXX) ../utils/AtlasStyle.C ../utils/AtlasLabels.C ../utils/AtlasUtils.C ../utils/Utils.C \
	$(CXXFLAGS) $(ROOTLIBS) $@.C -o $@

# Clean up: remove executables and outdated files.
clean:
	rm -rf *.o
	rm -rf calc_covmatrix
	rm -rf calc_tot_matrices
	
