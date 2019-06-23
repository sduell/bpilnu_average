/*
 *  HiggsCombiner: Florian Bernlochner
 */

#ifndef _HiggsPlot
#define _HiggsPlot


#include "../utils/Utils.h"

#include "HiggsDifferentialCommon.h"
#include "HiggsDifferentialInput.h"
#include "HiggsDifferentialWorkspace.h"

#include "RooFitResult.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TBox.h"
class HiggsPlot {
    
public:
    
    // Constructors
    HiggsPlot();
    HiggsPlot(settings set, HiggsInput *input, HiggsWorkspace *workspace,TFile *fout);
    
    // Plot the combination result
    void Plot(RooFitResult *result, RooFormulaVar *nll, vector<TGraph*>* scans);
    void PlotDifferentialCrossSections(RooFitResult *result);
    void PlotNormalizedCrossSections(RooFitResult *result);
    void PlotNuisanceParameters(RooFitResult *result);
    void Projections(RooFitResult* result, RooFormulaVar* nll);
    void PlotNuisanceParameterProjection(RooFitResult* result, RooFormulaVar* nll);
    void PlotNuisanceParameterProjection_Dummy(RooFitResult* result, RooFormulaVar* nll);
    void PlotScans(RooFitResult* result, RooFormulaVar* nll, vector<TGraph*>* scans);
    void PlotNuisanceScans(RooFitResult* result, RooFormulaVar* nll, vector<TGraph*>* scans);
    void PlotNPScans(RooFitResult* result, RooFormulaVar* nll, vector<TGraph*>* scans);
    
    void PlotComparison(RooFitResult *result);
    
    // Just some random plot functions for now
    void PlotMeasurements();
    void PlotDifferentialMeasurements();
    void PlotCovariances();
    void PlotCorrelations();
    
    void Finalize() { _can->Print(_ps+"]"); };
    
    
private:
    
    settings _set;
    TFile *_fout;
    
    HiggsInput *_input; HiggsWorkspace *_workspace;
    
    TCanvas *_can;
    Str _ps;
    
    VecD Unfold(input in);
    VecD UnfoldStatError(input in, int index);
    VecD UnfoldSysError(input in, int index);
    VecD UnfoldError(input in, int index);
    
    VecD UnfoldNormalized(input in);
    VecD UnfoldNormalizedError(input in, int index);
    
    VecD ShapeErrors(VecD xsec, VecD errors, TMatrixD Cov);
    
    
};

#endif
