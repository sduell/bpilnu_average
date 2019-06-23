/*
 *  HiggsCombiner: Florian Bernlochner
 */

#ifndef _HiggsWorkspace
#define _HiggsWorkspace


#include "../utils/Utils.h"
#include "HiggsDifferentialCommon.h"
#include "HiggsDifferentialInput.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooMultiVarGaussian.h"
#include "RooDataHist.h"
#include "RooAbsReal.h"
#include "RooHistPdf.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TCanvas.h"

#include <map>

using namespace RooFit;

class HiggsWorkspace {
    
public:
    
    // Constructors
    HiggsWorkspace();
    HiggsWorkspace(settings set, HiggsInput *input);
    
    // Return pointer to workspace
    RooWorkspace* GetWorkspace() { return _workspace; }
    
    TMatrixDSym* GetSyscovMatrix() { return _SyscovMatrix; }
    TMatrixDSym* GetStatcovMatrix() { return _StatcovMatrix; }
    TMatrixDSym* GetTotcovMatrix() { return _TotcovMatrix; }
    
    TMatrixDSym* GetSyscorMatrix() { return _SyscorMatrix; }
    TMatrixDSym* GetStatcorMatrix() { return _StatcorMatrix; }
    TMatrixDSym* GetTotcorMatrix() { return _TotcorMatrix; }
    
    
    int GetCovPos(int index) { return index > _start_bin.size() ? -99 :  _start_bin[index]; };
    
    int GetNMaxBins() { return _max_bins; }
    VecD GetMaxBins() { return _input->GetInput()[_finest_index].Bins; }
    
private:
    
    void SetupWorkspace();
    void SetupMultiWorkspace();
    void SetupSimplifiedWorkspace();
    void SetupDemoWorkspace();
    
    void ConstructCovariance();
    
    TMatrixDSym MergeCovariance(TMatrixDSym, TMatrixDSym);
    
    settings _set;
    HiggsInput *_input;
    
    RooWorkspace *_workspace;
    
    std::map<Str,TMatrixDSym *> _cov_mult;
    StrV _sv_source_mul, _sv_source_const;
    
    int _n_meas; VecD _start_bin; int _finest_index, _max_bins;
    
    
    TMatrixDSym *_SyscovMatrix, *_StatcovMatrix, *_TotcovMatrix;
    TMatrixDSym *_SyscorMatrix, *_StatcorMatrix, *_TotcorMatrix;
    
};

#endif