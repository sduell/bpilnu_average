/*
 *  Fitter
 */

#include "../utils/Utils.h"
#include <fstream>

#include "FitClass.h"
#include "PlotClass.h"
//#include "TheoryClass.h"

// ---------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    
    Printf("\n\n> Fitter -- doing stuff \n\n");
    
    SetAtlasStyle(); gStyle->SetPalette(1); gStyle->SetHistMinimumZero();
    
    // Reads in the settings, locations of files, etc
    Str SettingsFileName = argc == 1 ? "project.config" : argv[1];
    TEnv *set = OpenSettingsFile(SettingsFileName);
//    PlotClass *plot=new PlotClass(set);

    PlotClass *plot=new PlotClass(set);

    string toyopt=set->GetValue("GenerateToys"," ");

    VecVecD toyfitpars;
    VecVecD toyfiterrs;
    Str dbkpv=set->GetValue("Data","");
    vector<TVectorD> toyvectors;
    VecD toychi2;
    VecVecD covvalues;
    VecVecD rannum[atoi(set->GetValue("Toynum",""))];
    VecD rannum2;
    VecD toylcsrpoints;

//This is the normal fit procedure
    FitClass *fit = new FitClass(set);
    // Let's do the fit
    fit->Fit();
//    plot->Finalize();
    VecD pars=fit->getPars();
    VecD parerrs=fit->getParErrs();
    TMatrixDSym* FitCov=fit->getCovMat();
    double chi2=fit->getChi2();
    int numpar=fit->getNdf();
    double fitlcsr=fit->getfp0();
    cout << "fitlcsr is "<<fitlcsr<<endl;

    for(int i=0;i<pars.size();i++){
        cout << "parameter "<<i<<": "<<pars[i]<<" +- "<<parerrs[i]<<endl;
    }

    //Str newfitvals="Vub";
    //set->SetValue("FitPars",newfitvals);

    Str fitparvals=Form("%f %f %f",pars[1],pars[2],pars[3]);
    Str fitparerrs=Form("%f %f %f",parerrs[1],parerrs[2],parerrs[3]);
    //set->SetValue("apPars",fitparvals.Data());
    //set->SetValue("apErrs",fitparerrs.Data());


    //set->Print();
    //this is for marginalisation
    //FitClass *fit2 = new FitClass(set);
    //fit2->Fit();

    //VecD pars2=fit2->getPars();
    //VecD parerrs2=fit2->getParErrs();
    //TMatrixDSym* Cov2=fit2->getCovMat();
    //double secchi2=fit2->getChi2();
    //int numpar2=fit2->getNdf();

    //cout << "experimental error is: "<<parerrs2[0]<<" and theoretical error is: "<<sqrt(pow(parerrs[0],2.)-pow(parerrs2[0],2.))<<endl;

    //plot->PlotFitResult(pars,*Cov, chi2, numpar, parerrs2[0], sqrt(pow(parerrs[0],2.)-pow(parerrs2[0],2.)));
    //plot->PlotFitResult(pars,*FitCov, chi2, numpar, parerrs2[0], sqrt(pow(parerrs[0],2.)-pow(parerrs2[0],2.)));
/*    TMatrixDSym *parcov=new TMatrixDSym(pars.size());
    for(int i=0; i<pars.size();i++){
        for(int j=0; j<pars.size();j++){
            (*parcov)(i,j)=fit->GetCovarianceMatrixElement(i,j);
        }
    }
*/
    //plot->PlotFitResultMarg(pars, *FitCov, chi2, numpar);
    plot->PlotFitResult(pars,*FitCov, chi2, numpar);
cout << "plot finished"<<endl;

//set everything for the toys
    if(toyopt=="exp" || toyopt=="theo" || toyopt=="lcsr" || toyopt=="all"|| toyopt=="nolcsr"){
        cout << "Begin creating toys. Chosen option is "<<toyopt<<endl;
        set->SetValue("ap0.FitParVals", pars[1]);
        set->SetValue("ap1.FitParVals", pars[2]);
        set->SetValue("ap2.FitParVals", pars[3]);
        set->SetValue("LCSRval",fitlcsr);

        TRandom3 *rand = new TRandom3(2790);

        double x[pars.size()];
        for(int p=0;p<pars.size();p++){x[p]=pars[p];};
        TheoryClass* theo=new TheoryClass(set);
        theo->setAllPars(x);

        VecD low = VectorizeD(set->GetValue("LowerBinEdge","")," ");  
            if(low.size() == 0) 
                Fatal("No lower bin edge could be found.","");
        VecD high = VectorizeD(set->GetValue("UpperBinEdge","")," ");  
            if(high.size() == 0) 
                Fatal("No upper bin edge could be found.",""); 

        //get the data in order to generate toys
        VecD values = VectorizeD(set->GetValue("Data","")," ");  
            if(values.size() == 0) 
                Fatal("No values could be found.","");

        Int_t MatDim=atoi(set->GetValue("MatrixDimension",""));
        TMatrixDSym* Cov=new TMatrixDSym(MatDim);

        for(int i=0;i<MatDim;i++){
            covvalues.push_back(VectorizeD(set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
            for(int j=0;j<MatDim;j++){
                (*Cov)(i,j)=covvalues[i][j];
            }
        }

        const double elements[9]={0.000169, 0.0003978, 0.00069888, 0.0003978, 0.01,0.054784, 0.00069888,0.054784,0.4096};
        TMatrixDSym* constraintsMatrix=new TMatrixDSym(3,elements);

        TVectorD val(values.size());
        TVectorD valfit(values.size());
        TVectorD lambda;
        //prep TVector to rotate values into eigenbasis
        for(int ibin=0; ibin<values.size();ibin++){
            val(ibin)=theo->BF(0.,low[ibin],high[ibin])/((high[ibin]-low[ibin]));
            valfit(ibin)=val(ibin);
        }

        Str bkpvals="";
        for(int ibin=0; ibin<values.size();ibin++){
            bkpvals+=Form("%e ", val(ibin));
        }

        set->SetValue("Data",bkpvals);

        TVectorD parval(pars.size()-1);
        TVectorD parlambda;
        Str bkppars="";

        for(int ipar=0; ipar < pars.size()-1; ipar++){
            parval(ipar)=pars[ipar+1];
            bkppars+=Form("%e ", parval(ipar));
        }

        TVectorD lcsrval(1);
        lcsrval(0)=fitlcsr;
        Str bkplcsr=Form("%f",fitlcsr);

        // matrix with eigenvector: transformation matrix between original covariance
        // and diagonal form: C' = P^-1 C P
        TMatrixD P = Cov->EigenVectors(lambda);
        TMatrixD Pinvert = P; Pinvert.Invert();
        //Covariance matrix of FLAG average
        TMatrixD Ppar = constraintsMatrix->EigenVectors(parlambda);
        TMatrixD Pparinvert = Ppar; Pparinvert.Invert();
        //Get the error from the LCSR point
        double erp=atof(set->GetValue("LCSRerrup",""));
        double erm=atof(set->GetValue("LCSRerrdown",""));
        double err=(erp+erm)/2.;
        //create toys to validate the fit
        int toynum=atoi(set->GetValue("Toynum",""));
        if(toyopt=="all"){
            for(int toyindex=0; toyindex<toynum; toyindex++){
                //first vary the FF pars and extract the new toy values
                parval *= Pparinvert;
                for(int jpar=0; jpar<pars.size()-1;jpar++){
                    double rpnum=rand->Gaus(0.,1.);
                    parval(jpar)+=rpnum*sqrt(parlambda(jpar));
                }
                parval *= Ppar;
                Str newparvals="";
                for(int jpar=0; jpar<pars.size()-1;jpar++){
                    newparvals+=Form("%e ", parval(jpar));
                }
                set->SetValue("apPars", newparvals.Data());
                set->SetValue("ap0.FitParVals", parval(0));
                set->SetValue("ap1.FitParVals", parval(1));
                set->SetValue("ap2.FitParVals", parval(2));
                //This is another method which doesn't give the desired results here
                /*
                //now update the theory with the varied parameters
                TheoryClass* toytheo=new TheoryClass(set);
                toytheo->setFFpars(parval);
                for(int jbin=0; jbin<values.size();jbin++){
                    val(jbin)=toytheo->BF(0.,low[jbin],high[jbin])/((high[jbin]-low[jbin]));
                }
                //now reset the FF parameter values for the constraint in the fit
                delete toytheo;
                set->SetValue("apPars",bkppars.Data());
                set->SetValue("ap0.FitParVals", pars[1]);
                set->SetValue("ap1.FitParVals", pars[2]);
                set->SetValue("ap2.FitParVals", pars[3]);*/
                //now smear the new toy data with the exp. covariance matrix
                val *= Pinvert;
                for(int jbin=0; jbin<values.size();jbin++){
                    double rvnum=rand->Gaus(0.,1.);
                    val(jbin)+=rvnum*sqrt(lambda(jbin));
                }
                val *= P;
                toyvectors.push_back(val);
                Str newvals="";
                for(int jbin=0; jbin<values.size();jbin++){
                    newvals+=Form("%e ", val(jbin));
                }
                set->SetValue("Data", newvals.Data());
                //now vary the LCSR point which is independent of the other two parts
                double rlnum=rand->Gaus(0.,1.);
                lcsrval(0)+=rlnum*err;
                Str newlcsr=Form("%f ", lcsrval(0));
                set->SetValue("LCSRval", newlcsr.Data());
                //now start the fit and save the result
                FitClass *toyfit = new FitClass(set);
                toyfit->Fit();
                toyfitpars.push_back(toyfit->getPars());
                toyfiterrs.push_back(toyfit->getParErrs());
                toychi2.push_back(toyfit->getChi2());
                plot->setSettings(set);
                //delete the fitclass and reset the values to the nominal fit values
                delete toyfit;
                set->SetValue("Data",bkpvals.Data());
                set->SetValue("apPars",bkppars.Data());
                set->SetValue("ap0.FitParVals", pars[1]);
                set->SetValue("ap1.FitParVals", pars[2]);
                set->SetValue("ap2.FitParVals", pars[3]);
                set->SetValue("LCSRval",fitlcsr);
                for(int ipar=0; ipar < pars.size()-1; ipar++){
                    parval(ipar)=pars[ipar+1];
                }
                for(int ibin=0; ibin<values.size();ibin++){
                    val(ibin)=valfit(ibin);
                }
                lcsrval(0)=fitlcsr;                          
            }      
        }

        else if(toyopt=="exp"){
            for(int toyindex=0; toyindex<toynum; toyindex++){
                //smear the new toy data with the exp. covariance matrix
                val *= Pinvert;
                for(int jbin=0; jbin<values.size();jbin++){
                    double rvnum=rand->Gaus(0.,1.);
                    val(jbin)+=rvnum*sqrt(lambda(jbin));
                }
                val *= P;
                toyvectors.push_back(val);
                Str newvals="";
                for(int jbin=0; jbin<values.size();jbin++){
                    newvals+=Form("%e ", val(jbin));
                }
                set->SetValue("Data", newvals.Data());
                //now start the fit and save the result
                FitClass *toyfit = new FitClass(set);
                toyfit->Fit();
                toyfitpars.push_back(toyfit->getPars());
                toyfiterrs.push_back(toyfit->getParErrs());
                toychi2.push_back(toyfit->getChi2());
                plot->setSettings(set);
                //delete the fitclass and reset the values to the nominal fit values
                delete toyfit;
                set->SetValue("Data",bkpvals.Data());
                for(int ibin=0; ibin<values.size();ibin++){
                    val(ibin)=valfit(ibin);
                }        
            }      
        }

        else if(toyopt=="theo"){
            for(int toyindex=0; toyindex<toynum; toyindex++){
                //first vary the FF pars and extract the new toy values
                parval *= Pparinvert;
                for(int jpar=0; jpar<pars.size()-1;jpar++){
                    double rpnum=rand->Gaus(0.,1.);
                    parval(jpar)+=rpnum*sqrt(parlambda(jpar));
                }
                parval *= Ppar;
                Str newparvals="";
                for(int jpar=0; jpar<pars.size()-1;jpar++){
                    newparvals+=Form("%e ", parval(jpar));
                }
                set->SetValue("apPars", newparvals.Data());
                set->SetValue("ap0.FitParVals", parval(0));
                set->SetValue("ap1.FitParVals", parval(1));
                set->SetValue("ap2.FitParVals", parval(2));
                //now reset the FF parameter values for the constraint in the fit
                toyvectors.push_back(val);
                //now start the fit and save the result
                FitClass *toyfit = new FitClass(set);
                toyfit->Fit();
                toyfitpars.push_back(toyfit->getPars());
                toyfiterrs.push_back(toyfit->getParErrs());
                toychi2.push_back(toyfit->getChi2());
                plot->setSettings(set);
                //plot->PlotToyResult(val,pars,*FitCov, chi2, numpar);
                //if(toyindex==8){plot->Finalize();
                //Fatal("here","");}
                //delete the fitclass and reset the values to the nominal fit values
                delete toyfit;
                set->SetValue("apPars",bkppars.Data());
                set->SetValue("ap0.FitParVals", pars[1]);
                set->SetValue("ap1.FitParVals", pars[2]);
                set->SetValue("ap2.FitParVals", pars[3]);
                set->SetValue("Data",bkpvals.Data());
                for(int ipar=0; ipar < pars.size()-1; ipar++){
                    parval(ipar)=pars[ipar+1];
                }
                for(int ibin=0; ibin<values.size();ibin++){
                    val(ibin)=valfit(ibin);
                }
            }      
        }

        if(toyopt=="lcsr"){
            for(int toyindex=0; toyindex<toynum; toyindex++){
                //vary the LCSR point which is independent of the other two parts
                double rlnum=rand->Gaus(0.,1.);
                lcsrval(0)+=rlnum*err;
                Str newlcsr=Form("%f ", lcsrval(0));
                //set->SetValue("LCSRval", newlcsr.Data());
                //now start the fit and save the result
                FitClass *toyfit = new FitClass(set);
                toyfit->Fit();
                toyfitpars.push_back(toyfit->getPars());
                toyfiterrs.push_back(toyfit->getParErrs());
                toychi2.push_back(toyfit->getChi2());
                plot->setSettings(set);
                //delete the fitclass and reset the values to the nominal fit values
                delete toyfit;
                set->SetValue("LCSRval",fitlcsr);
                lcsrval(0)=fitlcsr;                          
            }      
        }
        if(toyopt=="nolcsr"){
            for(int toyindex=0; toyindex<toynum; toyindex++){
                //first vary the FF pars and extract the new toy values
                parval *= Pparinvert;
                for(int jpar=0; jpar<pars.size()-1;jpar++){
                    double rpnum=rand->Gaus(0.,1.);
                    parval(jpar)+=rpnum*sqrt(parlambda(jpar));
                }
                parval *= Ppar;
                Str newparvals="";
                for(int jpar=0; jpar<pars.size()-1;jpar++){
                    newparvals+=Form("%e ", parval(jpar));
                }
                set->SetValue("apPars", newparvals.Data());
                set->SetValue("ap0.FitParVals", parval(0));
                set->SetValue("ap1.FitParVals", parval(1));
                set->SetValue("ap2.FitParVals", parval(2));
                //This is another method which doesn't give the desired results here
                /*
                //now update the theory with the varied parameters
                TheoryClass* toytheo=new TheoryClass(set);
                toytheo->setFFpars(parval);
                for(int jbin=0; jbin<values.size();jbin++){
                    val(jbin)=toytheo->BF(0.,low[jbin],high[jbin])/((high[jbin]-low[jbin]));
                }
                //now reset the FF parameter values for the constraint in the fit
                delete toytheo;
                set->SetValue("apPars",bkppars.Data());
                set->SetValue("ap0.FitParVals", pars[1]);
                set->SetValue("ap1.FitParVals", pars[2]);
                set->SetValue("ap2.FitParVals", pars[3]);*/
                //now smear the new toy data with the exp. covariance matrix
                val *= Pinvert;
                for(int jbin=0; jbin<values.size();jbin++){
                    double rvnum=rand->Gaus(0.,1.);
                    val(jbin)+=rvnum*sqrt(lambda(jbin));
                }
                val *= P;
                toyvectors.push_back(val);
                Str newvals="";
                for(int jbin=0; jbin<values.size();jbin++){
                    newvals+=Form("%e ", val(jbin));
                }
                set->SetValue("Data", newvals.Data());
                //now start the fit and save the result
                FitClass *toyfit = new FitClass(set);
                toyfit->Fit();
                toyfitpars.push_back(toyfit->getPars());
                toyfiterrs.push_back(toyfit->getParErrs());
                toychi2.push_back(toyfit->getChi2());
                plot->setSettings(set);
                //delete the fitclass and reset the values to the nominal fit values
                delete toyfit;
                set->SetValue("Data",bkpvals.Data());
                set->SetValue("apPars",bkppars.Data());
                set->SetValue("ap0.FitParVals", pars[1]);
                set->SetValue("ap1.FitParVals", pars[2]);
                set->SetValue("ap2.FitParVals", pars[3]);
                for(int ipar=0; ipar < pars.size()-1; ipar++){
                    parval(ipar)=pars[ipar+1];
                }
                for(int ibin=0; ibin<values.size();ibin++){
                    val(ibin)=valfit(ibin);
                }                         
            }      
        }
    }    

cout << "Toys created "<<endl;



set->SetValue("Data",dbkpv.Data());











    // //Add LCSR point for variation

    //     if(toyopt=="exp"){

    //     double x[pars.size()];
    //     for(int p=0;p<pars.size();p++){x[p]=pars[p];};
    //     TheoryClass* theo=new TheoryClass(set);
    //     theo->setAllPars(x);

    //     VecD low = VectorizeD(set->GetValue("LowerBinEdge","")," ");  
    //         if(low.size() == 0) 
    //             Fatal("No lower bin edge could be found.","");
    //     VecD high = VectorizeD(set->GetValue("UpperBinEdge","")," ");  
    //         if(high.size() == 0) 
    //             Fatal("No upper bin edge could be found.",""); 

    //     //get the data in order to generate toys
    //     VecD values = VectorizeD(set->GetValue("Data","")," ");  
    //         if(values.size() == 0) 
    //             Fatal("No values could be found.","");

    //     Int_t MatDim=atoi(set->GetValue("MatrixDimension",""));
    //     TMatrixDSym* Cov=new TMatrixDSym(MatDim);

    //     for(int i=0;i<MatDim;i++){
    //         covvalues.push_back(VectorizeD(set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
    //         for(int j=0;j<MatDim;j++){
    //             (*Cov)(i,j)=covvalues[i][j];
    //         }
    //     }        

    //     TRandom3 *rand = new TRandom3(2790);

    //     TVectorD val(values.size());
    //     TVectorD lambda;
    //     //prep TVector to rotate values into eigenbasis
    //     for(int ibin=0; ibin<values.size();ibin++){
    //         values[ibin]=theo->BF(0.,low[ibin],high[ibin])/((high[ibin]-low[ibin]));
    //         val(ibin)=values[ibin];
    //     }

    //     Str bkpvals="";
    //     for(int ibin=0; ibin<values.size();ibin++){
    //         bkpvals+=Form("%e ", val(ibin));
    //     }

    //     // matrix with eigenvector: transformation matrix between original covariance
    //     // and diagonal form: C' = P^-1 C P
    //     TMatrixD P = Cov->EigenVectors(lambda);
    //     TMatrixD Pinvert = P; Pinvert.Invert();
    //     //create toys to validate the fit
    //     int toynum=atoi(set->GetValue("Toynum",""));
    //     //lambda.Print();
    //     for(int toyindex=0; toyindex<toynum; toyindex++){
    //         val *= Pinvert;
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             double rnum=rand->Gaus(0.,1.);
    //             val(ibin)+=rnum*sqrt(lambda(ibin));
    //             cout << "lambda "<< ibin << " is "<<lambda(ibin)<<endl;
    //             //rannum[toyindex].push_back(rnum);
    //             //if(ibin==0) rannum2.push_back(rnum);
    //         }
    //         val *= P;
    //         toyvectors.push_back(val);
    //         Str newvals="";
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             newvals+=Form("%e ", val(ibin));
    //             //if(ibin==0) rannum2.push_back(val(ibin)-values[ibin]);
    //         }
    //         set->SetValue("Data", newvals.Data());
    //         //set->Print();
    //         FitClass *fit = new FitClass(set);
    //         //fit->setToyBool(true);
    //         fit->Fit();
    //         toyfitpars.push_back(fit->getPars());
    //         toyfiterrs.push_back(fit->getParErrs());
    //         toychi2.push_back(fit->getChi2());
    //         plot->setSettings(set);
    //         //plot->PlotFitResult(fit->getPars(),*(fit->getCovMat()), fit->getChi2(), fit->getNdf());
    //         delete fit;
    //         set->SetValue("Data", bkpvals.Data());
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             val(ibin)=values[ibin];
    //         }
    //         /*if(toyindex>10){ 
    //             plot->Finalize();
    //             Fatal("here","");
    //         }*/
    //     }
    // }

    //         if(toyopt=="theo"){

    //     double x[pars.size()];
    //     for(int p=0;p<pars.size();p++){x[p]=pars[p];};
    //     TheoryClass* theo=new TheoryClass(set);
    //     theo->setAllPars(x);

    //     const double elements[9]={0.000169, 0.0003978, 0.00069888, 0.0003978, 0.01,0.054784, 0.00069888,0.054784,0.4096};

    //     TMatrixDSym* constraintsMatrix=new TMatrixDSym(3,elements);

    //     //get the data in order to generate toys
    //     VecD values = VectorizeD(set->GetValue("apPars","")," ");  
    //         if(values.size() == 0) 
    //             Fatal("No values could be found.","");

    //     TRandom3 *rand = new TRandom3(2790);

    //     TVectorD val(values.size());
    //     TVectorD lambda;
    //     //prep TVector to rotate values into eigenbasis
    //     for(int ibin=0; ibin<values.size();ibin++){
    //         val(ibin)=values[ibin];
    //     }

    //     Str bkpvals="";
    //     for(int ibin=0; ibin<values.size();ibin++){
    //         bkpvals+=Form("%e ", val(ibin));
    //     }

    //     // matrix with eigenvector: transformation matrix between original covariance
    //     // and diagonal form: C' = P^-1 C P
    //     TMatrixD P = constraintsMatrix->EigenVectors(lambda);
    //     TMatrixD Pinvert = P; Pinvert.Invert();
    //     //create toys to validate the fit
    //     int toynum=atoi(set->GetValue("Toynum",""));
    //     //lambda.Print();
    //     for(int toyindex=0; toyindex<toynum; toyindex++){
    //         val *= Pinvert;
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             double rnum=rand->Gaus(0.,1.);
    //             val(ibin)+=rnum*sqrt(lambda(ibin));
    //             //cout << "lambda "<< ibin << " is "<<lambda(ibin)<<endl;
    //             //rannum[toyindex].push_back(rnum);
    //             //if(ibin==0) rannum2.push_back(rnum);
    //         }
    //         val *= P;
    //         toyvectors.push_back(val);
    //         Str newvals="";
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             newvals+=Form("%e ", val(ibin));
    //             //if(ibin==0) rannum2.push_back(val(ibin)-values[ibin]);
    //         }
    //         set->SetValue("apPars", newvals.Data());
    //         set->SetValue("ap0.FitParVals", val(0));
    //         set->SetValue("ap1.FitParVals", val(1));
    //         set->SetValue("ap2.FitParVals", val(2));
    //         //set->Print();
    //         FitClass *fit = new FitClass(set);
    //         //fit->setToyBool(true);
    //         fit->Fit();
    //         toyfitpars.push_back(fit->getPars());
    //         toyfiterrs.push_back(fit->getParErrs());
    //         toychi2.push_back(fit->getChi2());
    //         plot->setSettings(set);
    //         //plot->PlotFitResult(fit->getPars(),*(fit->getCovMat()), fit->getChi2(), fit->getNdf());
    //         delete fit;
    //         set->SetValue("apPars", bkpvals.Data());
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             val(ibin)=values[ibin];
    //         }
    //         set->SetValue("ap0.FitParVals", val(0));
    //         set->SetValue("ap1.FitParVals", val(1));
    //         set->SetValue("ap2.FitParVals", val(2));
    //         /*if(toyindex>10){ 
    //             plot->Finalize();
    //             Fatal("here","");
    //         }*/
    //     }
    // }

    // if(toyopt=="lcsr"){

    //     double value = fitlcsr;
    //     cout << "original value is "<<value <<endl;

    //     TRandom3 *rand = new TRandom3(2790);

    //     Str bkpval=set->GetValue("LCSRval","");
    //     double erp=atof(set->GetValue("LCSRerrup",""));
    //     double erm=atof(set->GetValue("LCSRerrdown",""));
    //     double err=(erp+erm)/2.;

    //     //create toys to validate the fit
    //     int toynum=atoi(set->GetValue("Toynum",""));
    //     //lambda.Print();
    //     for(int toyindex=0; toyindex<toynum; toyindex++){

    //         double rnum=rand->Gaus(0.,1.);
    //         value+=rnum*err;
    //         TVectorD val(1);
    //         val(0)=value;
    //         toyvectors.push_back(val);
    //         toylcsrpoints.push_back(value);
    //         cout << "lcsr value is "<<value<<" "<<fitlcsr<<endl;
    //         Str newvals=Form("%e ", value);
 
    //         set->SetValue("LCSRval", newvals.Data());
    //         //set->Print();
    //         FitClass *fit = new FitClass(set);
    //         //fit->setToyBool(true);
    //         fit->Fit();
    //         toyfitpars.push_back(fit->getPars());
    //         toyfiterrs.push_back(fit->getParErrs());
    //         toychi2.push_back(fit->getChi2());
    //         plot->setSettings(set);
    //         //plot->PlotFitResult(fit->getPars(),*(fit->getCovMat()), fit->getChi2(), fit->getNdf());
    //         delete fit;
    //         set->SetValue("LCSRval", bkpval.Data());
    //         value=fitlcsr;
    //     }
    // }

    // if(toyopt=="all"){

    //     double x[pars.size()];
    //     for(int p=0;p<pars.size();p++){x[p]=pars[p];};
    //     TheoryClass* theo=new TheoryClass(set);
    //     theo->setAllPars(x);

    //     VecD low = VectorizeD(set->GetValue("LowerBinEdge","")," ");  
    //         if(low.size() == 0) 
    //             Fatal("No lower bin edge could be found.","");
    //     VecD high = VectorizeD(set->GetValue("UpperBinEdge","")," ");  
    //         if(high.size() == 0) 
    //             Fatal("No upper bin edge could be found.",""); 

    //     //get the data in order to generate toys
    //     VecD values = VectorizeD(set->GetValue("Data","")," ");  
    //         if(values.size() == 0) 
    //             Fatal("No values could be found.","");

    //     Int_t MatDim=atoi(set->GetValue("MatrixDimension",""));
    //     TMatrixDSym* Cov=new TMatrixDSym(MatDim);

    //     for(int i=0;i<MatDim;i++){
    //         covvalues.push_back(VectorizeD(set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
    //         for(int j=0;j<MatDim;j++){
    //             (*Cov)(i,j)=covvalues[i][j];
    //         }
    //     }

    //     const double elements[9]={0.000169, 0.0003978, 0.00069888, 0.0003978, 0.01,0.054784, 0.00069888,0.054784,0.4096};

    //     TMatrixDSym* constraintsMatrix=new TMatrixDSym(3,elements);

    //     //get the data in order to generate toys
    //     VecD parvalues = VectorizeD(set->GetValue("apPars","")," ");  
    //         if(parvalues.size() == 0) 
    //             Fatal("No parvalues could be found.","");

    //     double lcsrvalue = fitlcsr;

    //     TRandom3 *rand = new TRandom3(2790);

    //     TVectorD val(values.size());
    //     TVectorD lambda;
    //     //prep TVector to rotate values into eigenbasis
    //     for(int ibin=0; ibin<values.size();ibin++){
    //         values[ibin]=theo->BF(0.,low[ibin],high[ibin])/((high[ibin]-low[ibin]));
    //         val(ibin)=values[ibin];
    //     }

    //     Str bkpvals="";
    //     for(int ibin=0; ibin<values.size();ibin++){
    //         bkpvals+=Form("%e ", val(ibin));
    //     }


    //     TVectorD parval(parvalues.size());
    //     TVectorD parlambda;
    //     //prep TVector to rotate values into eigenbasis
    //     for(int ibin=0; ibin<parvalues.size();ibin++){
    //         parval(ibin)=parvalues[ibin];
    //     }

    //     Str bkpparvals="";
    //     for(int ibin=0; ibin<parvalues.size();ibin++){
    //         bkpparvals+=Form("%e ", parval(ibin));
    //     }

    //     Str bkplcsrval=set->GetValue("LCSRval","");
    //     double erp=atof(set->GetValue("LCSRerrup",""));
    //     double erm=atof(set->GetValue("LCSRerrdown",""));
    //     double err=(erp+erm)/2.;

    //     // matrix with eigenvector: transformation matrix between original covariance
    //     // and diagonal form: C' = P^-1 C P
    //     TMatrixD P = Cov->EigenVectors(lambda);
    //     TMatrixD Pinvert = P; Pinvert.Invert();

    //     TMatrixD Ppar = constraintsMatrix->EigenVectors(parlambda);
    //     TMatrixD Pparinvert = Ppar; Pparinvert.Invert();
    //     //create toys to validate the fit
    //     int toynum=atoi(set->GetValue("Toynum",""));
    //     //lambda.Print();
    //     for(int toyindex=0; toyindex<toynum; toyindex++){
    //         val *= Pinvert;
    //         parval *= Pparinvert;
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             double rnum=rand->Gaus(0.,1.);
    //             val(ibin)+=rnum*sqrt(lambda(ibin));
    //             //rannum[toyindex].push_back(rnum);
    //             //if(ibin==0) rannum2.push_back(rnum);
    //         }
    //         for(int ibin=0; ibin<parvalues.size();ibin++){
    //             double rnum=rand->Gaus(0.,1.);
    //             parval(ibin)+=rnum*sqrt(parlambda(ibin));
    //             //cout << "lambda "<< ibin << " is "<<lambda(ibin)<<endl;
    //             //rannum[toyindex].push_back(rnum);
    //             //if(ibin==0) rannum2.push_back(rnum);
    //         }
    //         val *= P;
    //         parval *= Ppar;

    //         double rlcsrnum=rand->Gaus(0.,1.);
    //         lcsrvalue+=rlcsrnum*err;
    //         TVectorD lcsrval(1);
    //         lcsrval(0)=lcsrvalue;

    //         toyvectors.push_back(val);
    //         toyvectors.push_back(parval);
    //         toyvectors.push_back(lcsrval);
    //         Str newvals="";
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             newvals+=Form("%e ", val(ibin));
    //             //if(ibin==0) rannum2.push_back(val(ibin)-values[ibin]);
    //         }
    //         set->SetValue("Data", newvals.Data());

    //         Str newparvals="";
    //         for(int ibin=0; ibin<parvalues.size();ibin++){
    //             newparvals+=Form("%e ", parval(ibin));
    //             //if(ibin==0) rannum2.push_back(val(ibin)-values[ibin]);
    //         }
    //         set->SetValue("apPars", newparvals.Data());
    //         set->SetValue("ap0.FitParVals", parval(0));
    //         set->SetValue("ap1.FitParVals", parval(1));
    //         set->SetValue("ap2.FitParVals", parval(2));

    //         Str newlcsrvals=Form("%e ", lcsrvalue);
    //         toylcsrpoints.push_back(lcsrvalue);
    //         set->SetValue("LCSRval", newlcsrvals.Data());

    //         //set->Print();
    //         FitClass *fit = new FitClass(set);
    //         //fit->setToyBool(true);
    //         fit->Fit();
    //         toyfitpars.push_back(fit->getPars());
    //         toyfiterrs.push_back(fit->getParErrs());
    //         toychi2.push_back(fit->getChi2());
    //         plot->setSettings(set);
    //         //plot->PlotFitResult(fit->getPars(),*(fit->getCovMat()), fit->getChi2(), fit->getNdf());
    //         delete fit;
    //         set->SetValue("Data", bkpvals.Data());
    //         for(int ibin=0; ibin<values.size();ibin++){
    //             val(ibin)=values[ibin];
    //         }

    //         set->SetValue("apPars", bkpparvals.Data());
    //         for(int ibin=0; ibin<parvalues.size();ibin++){
    //             parval(ibin)=parvalues[ibin];
    //         }
    //         set->SetValue("ap0.FitParVals", parval(0));
    //         set->SetValue("ap1.FitParVals", parval(1));
    //         set->SetValue("ap2.FitParVals", parval(2));

    //         set->SetValue("LCSRval", bkplcsrval.Data());
    //         lcsrvalue=fitlcsr;
    //         /*if(toyindex>10){ 
    //             plot->Finalize();
    //             Fatal("here","");
    //         }*/
    //     }
    // }

    //plot->PlotHistogram(rannum2);
    //plot->PlotToyPulls(pars, parerrs, toyfitpars);
    //plot->PlotToyPulls(pars, parcov, toyfitpars);
    //plot->Plotfpvalues(fitlcsr, toylcsrpoints);
    if(toyopt=="exp" || toyopt=="theo" || toyopt=="lcsr" || toyopt=="all" || toyopt=="nolcsr"){
    plot->Drawpval(toychi2);
    plot->PlotToyPulls(pars, toyfiterrs, toyfitpars);
    if(toyopt=="exp") plot->PlotFitResultMarg(pars,*FitCov, chi2, numpar);
    }
    plot->Finalize();

    //fit->Close();
    return 0;
}

