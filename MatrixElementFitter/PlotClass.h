#ifndef _PlotClass
#define _PlotClass

// -------------------------------------------------------------------------------------------
// PlotClass: a simple c++ interface to minuit
// Stephan and Florian
// -------------------------------------------------------------------------------------------

#include <vector>
#include "TString.h"

#include "TMatrixD.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "THStack.h"
#include "TRandom3.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TheoryClass.h"
#include "TLegend.h"

#include "../utils/Utils.h"


class PlotClass{
     public: 
      PlotClass();
      PlotClass(TEnv* set);

	private:
    
    TEnv *_set;
    TFile *_fout;
    
    TCanvas *_can;
    Str _ps;

    VecD _pull;
    TFile *_file;

    TheoryClass* _theory;

    public:
      void PlotFitResult(VecD pars, TMatrixDSym Cov, double chi2, int ndf);
      void PlotFitResult(VecD pars, TMatrixDSym Cov, double chi2, int ndf, double exerr, double theoerr);
      void PlotFitResultMarg(VecD pars, TMatrixDSym Cov, double chi2, int ndf);
      void calcVariatedFF(int variedPar, double NsigmaFF);
      void calcVariatedFF(int variedPar, double NsigmaFF, TVectorD* FFval, TMatrixDSym Cov);

      void setSettings(TEnv* set){_set=set;};

      void PlotHistogram(VecD vals){
        gStyle->SetOptStat(1101);
        TH1D *hist =new TH1D("hist","hist", vals.size()/100.,-0.003,0.003);
        for(int i=0; i<vals.size(); i++) hist->Fill(vals[i]);
        hist->Draw();
        _can->Print(_ps);
        _can->Clear();
        delete hist;
        };

      void PlotToyResult(TVectorD datvec, VecD pars, TMatrixDSym Cov, double chi2, int ndf);

      void Drawpval(VecD chi2);
      void Plotfpvalues(double mean, VecD pars);
      void PlotFitSteps(TVectorD bin, TVectorD data, TVectorD theo,TVectorD errvec,TVectorD erryvec,double *x, double chi2);
      void PlotTheoryPred(TVectorD bin, double *x);
      void PlotRatePred(TVectorD bin, double *x);
      void PlotMatrix(TMatrixD* mat);
      void PlotMatrix(TMatrixDSym* mat);
      void PlotToyPulls(VecD pars, VecD parerrs, VecVecD toypars);
      void PlotToyPulls(VecD pars, VecVecD parerrs, VecVecD toypars);
      void PlotToyPulls(VecD pars, TMatrixDSym* Cov, VecVecD toypars);
      void Finalize(){cout << "Finalizing plot"<<endl; _can->Print(_ps+"]"); delete _file;};
};
#endif