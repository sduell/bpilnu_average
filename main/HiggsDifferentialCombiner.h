/*
 *  HiggsCombiner: Florian Bernlochner
 */

#ifndef _HiggsCombiner
#define _HiggsCombiner

#include "../utils/Utils.h"
#include "HiggsDifferentialCommon.h"
#include "HiggsDifferentialInput.h"
#include "HiggsDifferentialWorkspace.h"
#include "HiggsDifferentialPlot.h"

#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooFitResult.h"

class HiggsCombiner {

   public:

      // Constructors
      HiggsCombiner();
      HiggsCombiner(settings set);
      
      void Combine();
      void ManualScan();
      void NLL_yy_ZZ();
      
      void SetStatusNPs(bool status = true);
    
      RooFitResult *GetCombinationResult() { return _result; };
                  
   private:

      TFile* _fout; // output root file to save results to
   
      settings _set; 
      
      HiggsInput *_input; 
      HiggsWorkspace *_workspace;
      HiggsPlot *_plot;
      
      RooFitResult *_result, *_stat_result;
      RooFormulaVar* _nll;

      RooFitResult *_result_yy, *_stat_result_yy;
      RooFormulaVar* _nll_yy;
      RooFitResult *_result_zz, *_stat_result_zz;
      RooFormulaVar* _nll_zz;

      vector<TGraph*>* _scans;

};

#endif
