// -------------------------------------------------------------------------------------------
// MatrixElementFitter: a simple c++ interface to Minuit
// Stephan and Florian
// -------------------------------------------------------------------------------------------
// Include fitstep graphs with TGraph, putting in the theoretical formula with errorbands
#include <iostream>
#include <fstream>

#include "TF1.h"
#include "TMath.h"
#include "TMatrixTSym.h"

#include "FitClass.h"

FitClass *FitClass::_gblInstance;

using namespace std;

FitClass::FitClass(TEnv* set) {

	_set=set;

   Printf("\n\n\n FitClass -- initializing \n\n");

    //_plot=new PlotClass(set);

	  _bins = VectorizeD(_set->GetValue("Bins","")," ");  
            if(_bins.size() == 0) 
                Fatal("No binvalues could be found.","");

    _low = VectorizeD(_set->GetValue("LowerBinEdge","")," ");  
            if(_bins.size() == 0) 
                Fatal("No lower bin edge could be found.","");
    _high = VectorizeD(_set->GetValue("UpperBinEdge","")," ");  
            if(_bins.size() == 0) 
                Fatal("No upper bin edge could be found.","");          

    _values = VectorizeD(_set->GetValue("Data","")," ");  
            if(_values.size() == 0) 
                Fatal("No values could be found.","");
 
    _theocovvalues = VectorizeD(_set->GetValue("TheoMatrix","")," "); 

    //for(int k=0;k<_theocovvalues.size();k++) cout << _theocovvalues[k]<<endl;
   Printf("\n\n\n FitClass -- Preparing Covariance Matrix \n\n");
   
    Int_t MatDim=atoi(_set->GetValue("MatrixDimension",""));
    Int_t TheoMatDim=atoi(_set->GetValue("TheoMatDim",""));

 		Cov=new TMatrixDSym(MatDim);


    for(int i=0;i<MatDim;i++){
        _covvalues.push_back(VectorizeD(_set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
    	for(int j=0;j<MatDim;j++){
/*        if(_covvalues.size()!=MatDim) 
        Fatal("Matrix Dimension and number of entries does not match","");*/
		//cout << "Covariance Matrix Element (" <<i<<","<<j<<") is:" <<_covvalues[i+j]<<endl;
    		//Comment this in for the real matrix
        (*Cov)(i,j)=_covvalues[i][j];
        //Comment this in for the diagonal matrix
        //if(i==j) (*Cov)(i,j)=_covvalues[i+j]; else (*Cov)(i,j)=0;
    	}
    }

    Covinvert=new TMatrixDSym((*Cov)*1e12);
    //Covinvert=new TMatrixDSym((*Cov));
    Covinvert->Invert();
    (*Covinvert)*=1e12;

    Cov->Print();
    Covinvert->Print();

    TheoCor=new TMatrixDSym(TheoMatDim);
   
    if(sqrt(_theocovvalues.size())!=TheoMatDim) 
          Fatal("Theory Matrix Dimension and number of entries does not match","");

    for(int i=0;i<TheoMatDim;i++){
      for(int j=0;j<TheoMatDim;j++){
    //cout << "Covariance Matrix Element (" <<i<<","<<j<<") is:" <<_covvalues[i+j]<<endl;
        (*TheoCor)(i,j)=_theocovvalues[i+j];
      }
    }
	
   // Set global instance
   setAsGlobalInstance(); 
	//_gblInstance=this; 
       

   // Some constants we later will/might need
   p0 = 0; p1 = 1;


    //Get all Fitparameters
    FitPars=Vectorize(_set->GetValue("FitPars","")," ");
   cout <<"read in parameters: "<<FitPars.size()<<endl;
    // Create the fitter objects and bind it to our minuit function
   _tfitter = new TFitter( FitPars.size() );    
   _tfitter->ExecuteCommand("SET PRINTOUT",&p1,1);
   _tfitter->SetFitMethod("chisquare");
   _tfitter->SetFCN(minuitFunction);  

   int ivar(0);
   for(Str FitPar : FitPars ) {
    cout << "run param. "<<FitPar<<endl;
   //Get the initial fitparameter values and their constraints if possible
   FitParVals.push_back(VectorizeD(_set->GetValue(Form("%s.FitParVals",FitPar.Data()),"")," "));
   cout <<"read in value, size is: "<<FitParVals.size()<<endl;
   //Get the initial errors
   FitParErrs.push_back(VectorizeD(_set->GetValue(Form("%s.FitParErrs",FitPar.Data()),"")," "));
   // Sort out options in case we have not specified any we assume there are no limits
    VecD Opts; Opts.push_back(FitParVals[ivar][0]); 
    //if(FitPar=="Vub"){Opts.push_back( 0. ); Opts.push_back( 0. );}
    Opts.push_back( 0. ); Opts.push_back( 0. );
    if ( FitParErrs[ivar].size() == 2) {
      Opts.push_back( FitParVals[ivar][0]-FitParErrs[ivar][0] ); Opts.push_back( FitParVals[ivar][0]+FitParErrs[ivar][1] ); 
    } else if( FitParErrs[ivar].size() == 1) {
      Opts.push_back( FitParVals[ivar][0]-FitParErrs[ivar][0] ); Opts.push_back( FitParVals[ivar][0]+FitParErrs[ivar][0] ); 
    } else {
      Opts.push_back( 0. ); Opts.push_back( 0. ); 
    }
    _tfitter->SetParameter(ivar++, FitPar,  Opts[0], 0.01, Opts[1],  Opts[2] );       
   }
/*  TMatrixDSym* Diag=new TMatrixDSym(TheoMatDim);
  for(int k=0;k<TheoMatDim;k++){
    (*Diag)(k,k)=FitParErrs[k][0];
  }
    //Diag->Print();
  //TheoCov=Diag*TheoCor*Diag;
  TheoCov=new TMatrixDSym((*TheoCor));
  TheoCov->SimilarityT((*Diag));
  TheoCovinvert=new TMatrixDSym(*TheoCov);
     
  TheoCovinvert->Invert();

  Printf("\n\n\n FitClass -- Theory Covariance Matrix built \n\n");
  
//  TheoCov->Print();
//  TheoCovinvert->Print();
*/
  _tfitter->ExecuteCommand("CALL FCN", &p1 ,1);



}


void FitClass::Fit() {
Printf("\n\n\n FitClass -- Fit initialized \n\n");
      // Fix the maximal number of calls
      _tfitter->SetMaxIterations(10000);

      // Execute a simple minimization scheme
    _tfitter->ExecuteCommand("SIMPLEX",&p0,0);
    _tfitter->ExecuteCommand("MIGRAD",&p0,0); 
    _tfitter->ExecuteCommand("IMPROVE",&p0,0);
      _tfitter->ExecuteCommand("MINOS",&p0,0);    
      _tfitter->ExecuteCommand("HESSE",&p0,0);    
      //Just for testing purposes;
      //_tfitter->ExecuteCommand("MINOS",&p0,0);  
      // Show Covariance  
      _tfitter->ExecuteCommand("SHOW COVARIANCE",&p0,0);
      // Execute the function one last time
      _tfitter->ExecuteCommand("CALL FCN", &p1 ,1);

Printf("\n\n\n FitClass -- Fit done \n\n");

for(int i=0;i<FitPars.size();i++)
{
  _pars.push_back(_tfitter->GetParameter(i));
  _pars_err.push_back(_tfitter->GetParError(i));
}
double* covvec=_tfitter->GetCovarianceMatrix();
_Cov=new TMatrixDSym(FitPars.size(), covvec);
_tfitter->PrintResults(1,0);
cout << "# d.o.f. is: "<<_tfitter->GetNumberFreeParameters()<<endl;
fitstatus=_tfitter->GetMinuit()->GetStatus();
cout << "TMinuit status: "<<_tfitter->GetMinuit()->GetStatus()<<endl;
cout <<"The Fitstatus is: "<< fitstatus<<endl;
cout <<"p values are: "<<p0<<" "<<p1<<endl;
if(_tfitter->GetMinuit()->GetStatus()!=0) Fatal("here","");

}

void FitClass::Fcn() {

   _tfitter->ExecuteCommand("CALL FCN", &p1 ,1);

}

// Minuit Function
// -------------------------------------------------------------------------------------------

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {

  //if(FitClass::getGlobalInstance()->getToyBool()==true) result = FitClass::getGlobalInstance()->toychisq(par);
  result = FitClass::getGlobalInstance()->chisq(par);
	//result = FitClass::_gblInstance->chisq(par);

}

double FitClass::toychisq(double *x) { 
 double value = 0;
 TheoryClass* _theory=new TheoryClass(_set);
 if(FitPars.size()==1) _theory->setMatrixElement(x);
 else if(FitPars.size()==4) _theory->setAllPars(x);
 //_theory->setMatrixElement(x[0]);


//Covariance matrix only from FLAG
const double elements[9]={0.000169, 0.0003978, 0.00069888, 0.0003978, 0.01,0.054784, 0.00069888,0.054784,0.4096};

//TMatrixDSym constraintsMatrix(8,elements);
TMatrixDSym constraintsMatrix(3,elements);
//constraintsMatrix.SetTol(1.e-23);
//constraintsMatrix.Print();

  TVectorD binmid(_bins.size());
  TVectorD errvec(_bins.size());
  TVectorD erryvec(_bins.size());
  TVectorD binvec(_bins.size());
  TVectorD theovec(_bins.size());

  TMatrixD mDiff( _bins.size(), 1 );
  TMatrixD mCovTotal( _bins.size(), _bins.size() );

  for(int ibin = 0; ibin < _bins.size();ibin++){
    binvec(ibin)=_values[ibin];
    //theovec(ibin)=_theory->dBF(_bins[ibin]);
    theovec(ibin)=_theory->BF(0.,_low[ibin],_high[ibin])/((_high[ibin]-_low[ibin]));
    binmid(ibin)=_bins[ibin];
    mDiff(ibin,0)=binvec[ibin]-theovec[ibin];
    errvec(ibin)=sqrt((*Cov)(ibin,ibin));
    erryvec(ibin)=(_high[ibin]-_low[ibin])/2.;
 }

    double ch=0.;
    // chisq = ( f(x) - data )^T * Cov^-1 * ( f(x) - data )
    TMatrixD m1 = TMatrixD( mDiff,  TMatrixD::kTransposeMult,   *Covinvert  );
    TMatrixD m2 = TMatrixD( m1,     TMatrixD::kMult,            mDiff       );
    ch = m2(0,0);


  double testval=0;
for(int k=0;k<_bins.size();k++){
  testval+=pow(binvec(k)-theovec(k),2.)/(*Cov)(k,k);
}

 //TVectorD a=binvec-theovec;
 //a*=(*Covinvert);
// double c=a*(binvec-theovec);


 value+=ch;
 //Covinvert->Print();

// cout << "calculating Theoretical chi2"<<endl;

  TMatrixD mPars( FitPars.size()-1, 1 );

  TVectorD parvec(FitPars.size()-1);
  TVectorD theoparvec(FitPars.size()-1);

  for(int ivar = 1; ivar < FitPars.size();ivar++){
    theoparvec(ivar-1)=FitParVals[ivar][0];
    //theovec(ibin)=_theory->dBF(_bins[ibin]);
    parvec(ivar-1)=x[ivar];
    mPars(ivar-1,0)=FitParVals[ivar][0]-x[ivar];
 }
// theoparvec.Print();
// parvec.Print();
// mPars.Print();
 
// cout << "inverting constraints Matrix"<<endl;

TMatrixDSym InverseMatrix(constraintsMatrix);
InverseMatrix.Invert();
//InverseMatrix.Print();

 //TVectorD b=theoparvec-parvec;
 //b*=(InverseMatrix);
 //double d=b*(theoparvec-parvec);
 //value+=d;


      double ch2=0.;

if(FitPars.size()>1){
    // chisq = ( f(x) - data )^T * Cov^-1 * ( f(x) - data )
    TMatrixD p1 = TMatrixD( mPars,  TMatrixD::kTransposeMult,   InverseMatrix  );
    TMatrixD p2 = TMatrixD( p1,     TMatrixD::kMult,            mPars       );
    ch2 = p2(0,0);
//    cout << "second chi2 is "<<ch2<<endl;
    //value+=ch2;

//Adding the point at q2=0 from Bharucha
    double fp0=0.261; 
    double erp=0.020;
    double erm=0.023;
    double era=0.0215;
    double fpdif=fp0-_theory->fp(0.);

    double symval=0.;
    symval=pow(fpdif,2.)/pow(era,2.);
    double asymval=0.;
    if(fpdif<0) asymval=pow(fpdif,2.)/pow(erp,2.);
    else asymval=pow(fpdif,2.)/pow(erm,2.);

    //value+=asymval;
    //value+=symval;
}
    //cout << "symmetric error: "<<symval<<" asymmetric error: "<<asymval<<endl;

 //cout << "chi2 is "<<value<<" test w/o params: "<<ch<<" w/o correlations: "<<testval<<endl;

//cout << "Chi^2 is: "<<value<<" and the p-value is "<<TMath::Prob(value,12)<<endl;
//_plot->PlotMatrix(Cov);
//_plot->PlotMatrix(&constraintsMatrix);
//_plot->PlotFitSteps(binmid,binvec,theovec,errvec, erryvec ,x,value);
//_plot->PlotTheoryPred(binmid,x);
//_plot->PlotRatePred(binmid,x);
//_plot->Finalize();
//Fatal("","here");

/*  for(int ivar = 0; ivar < FitPars.size(); ivar++){
     //cout << "FitParVals are: "<<FitParVals[ivar][0]<<" with the errors "<<FitParErrs[ivar][0]<<" and the fit value "<<x[ivar]<<endl;
    value+=pow(FitParVals[ivar][0]-x[ivar],2.)/pow(FitParErrs[ivar][0],2.);  
  }*/
 chi2=value;
 return value;
}


 double FitClass::chisq(double *x) { 
 double value = 0;
 TheoryClass* _theory=new TheoryClass(_set);
 if(FitPars.size()==1) _theory->setMatrixElement(x);
 else if(FitPars.size()==4) _theory->setAllPars(x);
 //_theory->setMatrixElement(x[0]);


//Covariance matrix only from FLAG
const double elements[9]={0.000169, 0.0003978, 0.00069888, 0.0003978, 0.01,0.054784, 0.00069888,0.054784,0.4096};

//TMatrixDSym constraintsMatrix(8,elements);
TMatrixDSym constraintsMatrix(3,elements);
//constraintsMatrix.SetTol(1.e-23);
//constraintsMatrix.Print();

  TVectorD binmid(_bins.size());
  TVectorD errvec(_bins.size());
  TVectorD erryvec(_bins.size());
  TVectorD binvec(_bins.size());
  TVectorD theovec(_bins.size());

  TMatrixD mDiff( _bins.size(), 1 );
  TMatrixD mCovTotal( _bins.size(), _bins.size() );

  for(int ibin = 0; ibin < _bins.size();ibin++){
    binvec(ibin)=_values[ibin];
    //theovec(ibin)=_theory->dBF(_bins[ibin]);
    theovec(ibin)=_theory->BF(0.,_low[ibin],_high[ibin])/((_high[ibin]-_low[ibin]));
    binmid(ibin)=_bins[ibin];
    mDiff(ibin,0)=binvec[ibin]-theovec[ibin];
    errvec(ibin)=sqrt((*Cov)(ibin,ibin));
    erryvec(ibin)=(_high[ibin]-_low[ibin])/2.;
 }

    double ch=0.;
    // chisq = ( f(x) - data )^T * Cov^-1 * ( f(x) - data )
    TMatrixD m1 = TMatrixD( mDiff,  TMatrixD::kTransposeMult,   *Covinvert  );
    TMatrixD m2 = TMatrixD( m1,     TMatrixD::kMult,            mDiff       );
    ch = m2(0,0);


  double testval=0;
for(int k=0;k<_bins.size();k++){
  testval+=pow(binvec(k)-theovec(k),2.)/(*Cov)(k,k);
}

 //TVectorD a=binvec-theovec;
 //a*=(*Covinvert);
// double c=a*(binvec-theovec);


 value+=ch;
 //Covinvert->Print();

// cout << "calculating Theoretical chi2"<<endl;

  TMatrixD mPars( FitPars.size()-1, 1 );

  TVectorD parvec(FitPars.size()-1);
  TVectorD theoparvec(FitPars.size()-1);

  for(int ivar = 1; ivar < FitPars.size();ivar++){
    theoparvec(ivar-1)=FitParVals[ivar][0];
    //theovec(ibin)=_theory->dBF(_bins[ibin]);
    parvec(ivar-1)=x[ivar];
    mPars(ivar-1,0)=FitParVals[ivar][0]-x[ivar];
 }
// theoparvec.Print();
// parvec.Print();
// mPars.Print();
 
// cout << "inverting constraints Matrix"<<endl;

TMatrixDSym InverseMatrix(constraintsMatrix);
InverseMatrix.Invert();
//InverseMatrix.Print();

 //TVectorD b=theoparvec-parvec;
 //b*=(InverseMatrix);
 //double d=b*(theoparvec-parvec);
 //value+=d;


      double ch2=0.;

if(FitPars.size()>1){
    // chisq = ( f(x) - data )^T * Cov^-1 * ( f(x) - data )
    TMatrixD p1 = TMatrixD( mPars,  TMatrixD::kTransposeMult,   InverseMatrix  );
    TMatrixD p2 = TMatrixD( p1,     TMatrixD::kMult,            mPars       );
    ch2 = p2(0,0);
//    cout << "second chi2 is "<<ch2<<endl;
    value+=ch2;

//Adding the point at q2=0 from Bharucha
    //double fp0=0.261; 
    //double erp=0.020;
    //double erm=0.023;
    //double era=0.0215;
    double fp0=atof(_set->GetValue("LCSRval","")); 
    double erp=atof(_set->GetValue("LCSRerrup",""));
    double erm=atof(_set->GetValue("LCSRerrdown",""));
    double era=(erp+erm)/2.;
    double fpdif=fp0-_theory->fp(0.);

    double symval=0.;
    symval=pow(fpdif,2.)/pow(era,2.);
    double asymval=0.;
    if(fpdif<0) asymval=pow(fpdif,2.)/pow(erp,2.);
    else asymval=pow(fpdif,2.)/pow(erm,2.);

    //value+=asymval;
    value+=symval;
}
    //cout << "symmetric error: "<<symval<<" asymmetric error: "<<asymval<<endl;

 //cout << "chi2 is "<<value<<" test w/o params: "<<ch<<" w/o correlations: "<<testval<<endl;

//cout << "Chi^2 is: "<<value<<" and the p-value is "<<TMath::Prob(value,12)<<endl;
//_plot->PlotMatrix(Cov);
//_plot->PlotMatrix(&constraintsMatrix);
//_plot->PlotFitSteps(binmid,binvec,theovec,errvec, erryvec ,x,value);
//_plot->PlotTheoryPred(binmid,x);
//_plot->PlotRatePred(binmid,x);
//_plot->Finalize();
//Fatal("","here");

/*  for(int ivar = 0; ivar < FitPars.size(); ivar++){
     //cout << "FitParVals are: "<<FitParVals[ivar][0]<<" with the errors "<<FitParErrs[ivar][0]<<" and the fit value "<<x[ivar]<<endl;
    value+=pow(FitParVals[ivar][0]-x[ivar],2.)/pow(FitParErrs[ivar][0],2.);  
  }*/
fplus_0=_theory->fp(0.);
 chi2=value;
 return value;
}
