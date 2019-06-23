#ifndef _FitClass
#define _FitClass

// -------------------------------------------------------------------------------------------
// FitClass: a simple c++ interface to minuit
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

//#include "TheoryClass.h"
#include "PlotClass.h"

#include "../utils/Utils.h"


class FitClass{
     public: 
      FitClass(TEnv* set);
       
       // Global instances
       static FitClass* getGlobalInstance() { return _gblInstance; };
       void setAsGlobalInstance() { _gblInstance = this; }
       
       VecD getPars(){return _pars;};
       VecD getParErrs(){return _pars_err;};

       double GetCovarianceMatrixElement(int i, int j){return _tfitter->GetCovarianceMatrixElement(i,j);};

       double getChi2(){return chi2;};

       double getfp0(){return fplus_0;};

       double getFitStatus(){return fitstatus;};

       int getNdf(){return _tfitter->GetNumberFreeParameters();};

       TMatrixDSym* getCovMat(){return _Cov;};

       // Chisq function
       double chisq(double *x);

       double toychisq(double *x);

       // Execute the fit
       void Fit();

       // Evalute the Chi Square
       void Fcn();

       //Create Pulls for signal and background
       //void CreatePulls();

       void setToyBool(bool b){toybool=b;};
       bool getToyBool(){return toybool;};

       void Close(){_plot->Finalize();};

       // Return the private fitter object if requested
       TFitter* getFitter() { return _tfitter; } 
       
       VecD _pars, _pars_err;
              
       // Global instance pointer for class
       static FitClass* _gblInstance; 
     private:
      double chi2; double fplus_0; int fitstatus;
       // TEnv object to read out settings file
 	   TEnv *_set;
           TMatrixDSym *_Cov;
          // Define Covariance Matrix
           TMatrixDSym *Cov, *Covinvert;
          // Define Covariance Matrix of theory
           TMatrixDSym *TheoCor, *TheoCov, *TheoCovinvert;
 	   // Input class pointer
 	   VecD _values; VecD _bins; VecD _low; VecD _high; VecVecD _covvalues; VecD _theocovvalues; StrV FitPars; VecVecD FitParVals; VecVecD FitParErrs;
 	   // TFitter object to interface with minuit
       TFitter *_tfitter;

       bool toybool=false;

       PlotClass *_plot;

       // values to store current likelihood and minuit bits that are needed to set things up
       double _loglikelihood, p0,p1;
       // Number of bins of the binned likelihood
       int _nbins;
       // Number of expected events of the templates
       VecD _NexpTemplates;

       // TCanvas and plot file
       TCanvas *_can; Str _ps;
       
              
};

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg);

#endif
