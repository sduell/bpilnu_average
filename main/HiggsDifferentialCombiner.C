/*
 *  HiggsCombiner: Florian Bernlochner
 */

#include "HiggsDifferentialCombiner.h"

HiggsCombiner::HiggsCombiner() { }

HiggsCombiner::HiggsCombiner(settings set) {
    
    Printf("\n\n> Initializing HiggsCombiner Class ");
    
    _set = set;
    // Read in the measurements into an input object
    _input = new HiggsInput(set);
    // Now setup the workspace given this input
    _workspace = new HiggsWorkspace(set,_input);
    // Create output file for saving results to
    Str ofn = Form("results/%s_%s%s.root",_set.Variable.Data(),_set.IncludeNP ? "floatNPs" : "fixNPs",_set.ManualScan=="true" ? "_manualScan" : "");
    _fout = new TFile(ofn,"recreate");
    // Initialize the plotting class
    _plot = new HiggsPlot(set,_input,_workspace,_fout);
    // Initialise vector of manual scans of nll
    _scans = new vector<TGraph*>;
    // Set result pointers to zero
    _result = 0; _stat_result = 0; _stat_result_yy = 0; _stat_result_zz = 0;

}

// ---------------------------------------------------------------------------------------
// This is our fitter function where we minimize the combination likelihood function
void HiggsCombiner::Combine() {
    
    Printf("\n\n> HiggsCombiner: Combine input \n ");
    
    // Get combination pdf from workspace
    RooAbsPdf *pdf = _workspace->GetWorkspace()->pdf("pdf_tot");
    Printf("\n\n> HiggsCombiner: Took pdf from workspace \n ");
    //cout << "pdf address is: "<<pdf<<endl;
    //pdf->printMultiline();
    // Construct -2 log Likelihood function we want to minimize
    _nll = new RooFormulaVar("nll", "nll", "-2*log(@0)", RooArgSet(*pdf));
    Printf("\n\n> HiggsCombiner: nll setup \n ");
    RooMsgService::instance().setGlobalKillBelow(ERROR);
    Printf("\n\n> HiggsCombiner: Survived Global Kill \n ");
    // Create a minuit instance
    RooMinuit m(*_nll);
    
    // Specify Error level, 1 = 1 sigma, 4 = 2 sigma coverage
    m.setErrorLevel(1.0);
    // Strategy 2 is slower, but more reliable
    m.setStrategy(2);
    // Allow offsetting of likelihood funcations
    //m.setOffsetting(true);
    // Profile 1 enables migrad timer
    // m.setProfile(0);
    // Do the fit
    m.simplex();
    m.migrad();
    m.improve();
    m.hesse();
    m.minos();
    
    RooMsgService::instance().setGlobalKillBelow(INFO);
    
    // Save and show the fit result
    _result = m.save();
    _result->Print();
    
    cout << "Nuisance parameters activation: "<<_set.IncludeNP<<endl;
    // Fix nuisance parameters and refit
    if(_set.IncludeNP) {
        cout<< "NP included"<<endl;
     SetStatusNPs();
     cout<< "SetStatusNPs included"<<endl;
     // Redo the fit
     m.simplex();
     cout<< "Simplex done"<<endl;
     m.migrad();
     cout<< "Migrad done"<<endl;
     m.improve();
     cout<< "Improve done"<<endl;
     m.hesse();
     cout<< "Hesse done"<<endl;
     m.minos();    
     _stat_result = m.save();
     _stat_result->Print();
     //This is the normal option
     //SetStatusNPs(false); // release NPs
     //use this option to fix the NP and get only the statistical covariance
     SetStatusNPs(true);
     cout << "Checking if workspace really fixed the NP"<<endl;
    _stat_result = m.save();
    _stat_result->Print();
    m.simplex();
    m.migrad();
    m.improve();
    m.hesse();
    m.minos();
    _stat_result = m.save();
    _stat_result->Print();


    TMatrixDSym ErrorMatrix(13);
    VecD comb_sigma, comb_sigmaError; double tot_sigma, tot_sigmaError;
    for(int i = 0; i < _workspace->GetNMaxBins()-1; i++) {
        RooRealVar *sigma = (RooRealVar*)(_stat_result->floatParsFinal().find(Form("sigma_%d",i+1)));
        comb_sigma.push_back( sigma->getVal() ); comb_sigmaError.push_back( sigma->getError() );
        tot_sigma += sigma->getVal(); tot_sigmaError = sqrt(pow(tot_sigmaError,2.) + pow(sigma->getError(),2.));
        //(ErrorMatrix)(i,i) = sigma->getError()/20.;
    }
    TG *comb=DrawTG(comb_sigma,comb_sigmaError,_workspace->GetMaxBins(),kBlack, true ) ;
    for (int o=0;o<comb->GetN();o++) {comb->GetY()[o] *= 0.1; comb->GetEY()[o] *=0.1;}
    vector <int> irow {0,5,6,7,8,9,10,11,12,1,2,3,4};
    vector <int> jrow {0,5,6,7,8,9,10,11,12,1,2,3,4};

    cout << "The reconstructed statistical covariance Matrix is "<<endl;
    TMatrixDSym FCorwrongorder = _stat_result->correlationMatrix();
    TMatrixDSym FCor(13);
    for(int i=0; i<13/*FitCorMatrix.GetNrows()*/;i++){
        for(int j=0; j<13/*FitCorMatrix.GetNrows()*/; j++){
            (FCor)(i,j)=(FCorwrongorder)(irow[i],jrow[j]);
            ErrorMatrix(i,i)=comb->GetErrorY(i);
            cout << comb->GetY()[i]<< " pm "<< comb->GetErrorY(i)<<endl;
        }
    }    
    ofstream corout;
    corout.open("/home/sduell/bpilnu_new_average/CombinationCode/results/stat_cormat.dat", ios::out);

    for(int i=0; i<13/*FitCovMatrix.GetNrows()*/;i++){
        for(int j=0; j<13/*FitCovMatrix.GetNrows()*/; j++){
            corout <<(double)(FCor)(i,j)<<" ";
            //cout<< "Matrixelement should be "<<(double)(FCor)(i,j)<<endl;
        }
    }
    corout.close();

    FCor.SimilarityT(ErrorMatrix);
    //TMatrixD FCortmp=ErrorMatrix*FCor;
    //FCortmp = FCortmp*FCor;
    FCor.Print();
    ofstream covout;

    covout.open("/home/sduell/bpilnu_new_average/CombinationCode/results/stat_covmat.dat", ios::out);
    
    for(int i=0; i<13/*FitCovMatrix.GetNrows()*/;i++){
        for(int j=0; j<13/*FitCovMatrix.GetNrows()*/; j++){
            covout <<(double)(FCor)(i,j)<<" ";
            cout<< "Matrixelement should be "<<(double)(FCor)(i,j)<<endl;
        }
    }
    covout.close();
    TMatrixDSym Cov = _stat_result->covarianceMatrix();
    cout << "The statistical covariance Matrix is "<<endl;
    Cov.Print();
    }
        
    // Manual scan
    if ( _set.ManualScan=="true" )
        ManualScan();
    
    // Plot the combination
    _plot->Plot(_result, _nll, _scans );
    
    
    // Perform fits with individual measurement pdfs and write to file
    if(_set.DoPerformIndividualScan) {
        NLL_yy_ZZ();
        
        // Write results to file
        _fout->cd();
        _result->SetName("combined_result");
        _result->Write();
        if(_stat_result) {
         _stat_result->SetName("combined_stat_result");
         _stat_result->Write();
        } 
        _nll->SetName("combined_nll");
        _nll->Write();
        for (int i=0; i<_scans->size(); ++i) {
            _scans->at(i)->SetName( Form("scan_sigma%d",i+1) );
            _scans->at(i)->Write();
        }
        gROOT->cd();
    }

    _plot->Finalize(); // Close the PDF
    cout << "Fit finished successfully"<<endl;
    
    
}


// ---------------------------------------------------------------------------------------
// Manually scan minimum log-likelihood for each cross section
void HiggsCombiner::ManualScan(  ) {
    
    // Get pointer to workspace
    RooWorkspace *ws = _workspace->GetWorkspace();
    
    // Store scan results in vector of pointers to TGraphs
    
    _scans->clear();
    
    // Minimum nll from 'nominal' fit (for shifting likelihood curves to zero)
    double min_nll_fit = _result->minNll();
    
    // Loop over floated cross sections
    for(int i = 0; i < _workspace->GetNMaxBins()-1; i++) {
        VecD min_nll, scan_pnts;
        
        // Get fitted sigma
        RooRealVar *sigma = (RooRealVar*)(_result->floatParsFinal().find(Form("sigma_%d",i+1)));
        
        // Define range to scan from sigma +/- N x error (specified in config)
        double range_low = sigma->getVal() - _set.ScanNSigma * sigma->getError() > 0 ? sigma->getVal() - _set.ScanNSigma * sigma->getError() : 0.;
        double range_high = sigma->getVal() + _set.ScanNSigma * sigma->getError();
        
        printf("fitted sigma = %f, scanning range %f - %f \n",sigma->getVal(),range_low,range_high);
        
        // Loop over scan points
        for(int j = 0; j < _set.ScanNPoints; ++j) {
            double sig = range_low + ( j * (range_high - range_low) / (_set.ScanNPoints - 1) );
            
            // reset sigmas
            for(int k = 0; k < _workspace->GetNMaxBins()-1; ++k)
                ws->var(Form("sigma_%d",k+1))->setVal(_set.CrossSectionOptions[0]);
            
            ws->var(Form("sigma_%d",i+1))->setVal(sig);
            ws->var(Form("sigma_%d",i+1))->setConstant(kTRUE);
            
            printf("%f \n",sig);
            
            RooAbsPdf *pdf = ws->pdf("pdf_tot");
            RooFormulaVar nll("nll_tmp", "nll_tmp", "-2*log(@0)", RooArgSet(*pdf));
            
            // Create instance of minuit and do fits
            RooMinuit m(nll);
            m.setErrorLevel(1.0);
            m.setStrategy(2);
            // m.setOffsetting(true);
            m.simplex();
            m.migrad();
            m.improve();
            m.hesse();
            m.minos();
            
            m.save()->status();
            
            // Print results
            m.save()->Print();
            
            // Save results
            scan_pnts.push_back(sig);
            min_nll.push_back(m.save()->minNll() - min_nll_fit); // Shift to zero
            //      min_nll.push_back(m.save()->minNll());
            
        }
        _scans->push_back(new TGraph(_set.ScanNPoints,&scan_pnts[0],&min_nll[0]));
        
        // Let cross section float again
        ws->var(Form("sigma_%d",i+1))->setConstant(kFALSE);
    }
    if(true){
        VecD min_nll, scan_pnts;
        
        // Get fitted sigma
        RooRealVar *sigma = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_Dummy")));
        
        // Define range to scan from sigma +/- N x error (specified in config)
        double range_low = sigma->getVal() - _set.ScanNSigma * sigma->getError() < 0 ? sigma->getVal() - _set.ScanNSigma * sigma->getError() : 0.;
        double range_high = sigma->getVal() + _set.ScanNSigma * sigma->getError();
        
        printf("fitted Dummy = %f, scanning range %f - %f \n",sigma->getVal(),range_low,range_high);
        
        // Loop over scan points
        for(int j = 0; j < _set.ScanNPoints; ++j) {
            double sig = range_low + ( j * (range_high - range_low) / (_set.ScanNPoints - 1) );
            
            // reset NPs
            ws->var(Form("Uncert_Dummy"))->setVal(_set.CrossSectionOptions[0]);
            
            ws->var(Form("Uncert_Dummy"))->setVal(sig);
            ws->var(Form("Uncert_Dummy"))->setConstant(kTRUE);
            
            printf("%f \n",sig);
            
            RooAbsPdf *pdf = ws->pdf("pdf_tot");
            RooFormulaVar nll("nll_tmp", "nll_tmp", "-2*log(@0)", RooArgSet(*pdf));
            
            // Create instance of minuit and do fits
            RooMinuit m(nll);
            m.setErrorLevel(1.0);
            m.setStrategy(2);
            // m.setOffsetting(true);
            m.simplex();
            m.migrad();
            m.improve();
            m.hesse();
            m.minos();
            
            m.save()->status();
            
            // Print results
            m.save()->Print();
            
            // Save results
            scan_pnts.push_back(sig);
            min_nll.push_back(m.save()->minNll() - min_nll_fit); // Shift to zero
            //      min_nll.push_back(m.save()->minNll());
            
        }
        _scans->push_back(new TGraph(_set.ScanNPoints,&scan_pnts[0],&min_nll[0]));
        
        // Let cross section float again
        ws->var(Form("Uncert_Dummy"))->setConstant(kFALSE);

    }

    // Reset sigmas in workspace to minimum values
    for(int m = 0; m < _workspace->GetNMaxBins()-1; m++) {
        RooRealVar *sig = (RooRealVar*)(_result->floatParsFinal().find(Form("sigma_%d",m+1)));
        ws->var( Form("sigma_%d",m+1) )->setVal( sig->getVal() );
    }
    
    // Reset NPs in workspace to minimum values
    if(_set.CombinationMethod == "full" && _set.IncludeNP) {
        for(int k = 0; k < _input->GetInput().size(); k++) {
            input in = _input->GetInput()[k];
            


            for(int m = 0; m < in.BF_Vub_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.BF_Vcb_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.BF_Vcb_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.BF_Vcb_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.FF_Vub_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.FF_Vub_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.FF_Vub_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.FF_Vcb_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.FF_Vcb_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.FF_Vcb_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.general_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.general_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.general_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }






            for(int m = 0; m < in.Bkg_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.Bkg_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.Bkg_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
        }
    }
}

// ---------------------------------------------------------------------------------------
// Perform fits with individual measurement pdfs
void HiggsCombiner::NLL_yy_ZZ () {
    
    // TO DO: change to loop over inputs! (would be much nicer...)
    
    // Get pointer to workspace
    RooWorkspace *ws = _workspace->GetWorkspace();
    
    // Get ZZ pdf from workspace
    RooAbsPdf *pdf_zz = ws->pdf("pdf_ZZ_NPs");
    
    // Construct -2 log Likelihood function we want to minimize
    RooFormulaVar* _nll_zz = new RooFormulaVar("zz_nll", "zz_nll", "-2*log(@0)", RooArgSet(*pdf_zz));
    
    // Create instance of minuit and do fits
    RooMinuit m_zz(*_nll_zz);
    m_zz.setErrorLevel(1.0);
    m_zz.setStrategy(2);
    // m_zz.setOffsetting(true);
    m_zz.simplex();
    m_zz.migrad();
    m_zz.improve();
    m_zz.hesse();
    m_zz.minos();
    
    RooFitResult* _result_zz = m_zz.save();
    _result_zz->Print();
    
    // Fix nuisance parameters and refit
    if(_set.IncludeNP) {
     SetStatusNPs();
     // Redo the fit
     m_zz.simplex();
     m_zz.migrad();
     m_zz.improve();
     m_zz.hesse();
     m_zz.minos();    
     _stat_result_zz = m_zz.save();
     _stat_result_zz->Print();
     SetStatusNPs(false); // release NPs
    }    
    
    // Get yy pdf from workspace
    RooAbsPdf *pdf_yy = ws->pdf("pdf_yy_NPs");
    
    // Construct -2 log Likelihood function we want to minimize
    RooFormulaVar* _nll_yy = new RooFormulaVar("yy_nll", "yy_nll", "-2*log(@0)", RooArgSet(*pdf_yy));
    
    // Create instance of minuit and do fits
    RooMinuit m_yy(*_nll_yy);
    m_yy.setErrorLevel(1.0);
    m_yy.setStrategy(2);
    // m_yy.setOffsetting(true);
    m_yy.simplex();
    m_yy.migrad();
    m_yy.improve();
    m_yy.hesse();
    m_yy.minos();
    
    RooFitResult* _result_yy = m_yy.save();
    _result_yy->Print();

    // Fix nuisance parameters and refit
    if(_set.IncludeNP) {
     SetStatusNPs();
     // Redo the fit
     m_yy.simplex();
     m_yy.migrad();
     m_yy.improve();
     m_yy.hesse();
     m_yy.minos();    
     _stat_result_yy = m_yy.save();
     _stat_result_yy->Print();
     SetStatusNPs(false); // release NPs
    } 
    
    // Write results to file
    _fout->cd();
    _result_yy->SetName("yy_result");
    _result_yy->Write();    
    if(_stat_result_yy) {
    _stat_result_yy->SetName("yy_stat_result");
    _stat_result_yy->Write();    
    }
    _nll_yy->Write();
    _result_zz->SetName("zz_result");
    _result_zz->Write();
    if(_stat_result_zz) {
    _stat_result_zz->SetName("zz_stat_result");
    _stat_result_zz->Write();    
    }
    _nll_zz->Write();
    gROOT->cd();
    
    // Reset sigmas in workspace to minimum values
    for(int m = 0; m < _workspace->GetNMaxBins()-1; m++) {
        RooRealVar *sig = (RooRealVar*)(_result->floatParsFinal().find(Form("sigma_%d",m+1)));
        ws->var( Form("sigma_%d",m+1) )->setVal( sig->getVal() );
    }
    
    // Reset NPs in workspace to minimum values
    if(_set.CombinationMethod == "full" && _set.IncludeNP) {
        for(int k = 0; k < _input->GetInput().size(); k++) {
            input in = _input->GetInput()[k];


            for(int m = 0; m < in.BF_Vub_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.BF_Vcb_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.BF_Vcb_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.BF_Vcb_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.FF_Vub_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.FF_Vub_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.FF_Vub_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.FF_Vcb_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.FF_Vcb_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.FF_Vcb_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
            for(int m = 0; m < in.general_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.general_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.general_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }





            for(int m = 0; m < in.Bkg_Uncertainties.size(); m++) {
                RooRealVar *NP = (RooRealVar*)(_result->floatParsFinal().find(Form("Uncert_%s",in.Bkg_Uncertainties[m].Data())));
                ws->var( Form("Uncert_%s",in.Bkg_Uncertainties[m].Data()) )->setVal( NP->getVal() );
            }
        }
    }
    
}

// ---------------------------------------------------------------------------------------
// Fix or release NPs 

void HiggsCombiner::SetStatusNPs(bool status) {

  // Get pointer to workspace
  RooWorkspace *ws = _workspace->GetWorkspace();
  cout << "Got Workspace"<<endl;
  // Reset NPs status
  for(int k = 0; k < _input->GetInput().size(); k++) {
    cout << "Begin loop"<<endl;
   input in = _input->GetInput()[k];
   cout << "Got Input"<<endl;

   for(int m = 0; m < in.Bkg_Uncertainties.size(); m++) 
     ws->var( Form("Uncert_%s",in.Bkg_Uncertainties[m].Data()) )->setConstant(status);
     cout << "Bkg Status set"<<endl;
     cout << "next up is Vub BF: "<<in.BF_Vub_Uncertainties.size() <<" with status " <<status<<endl;
     //ws->Print();
    // cout << "Content is: "<< Form("Uncert_%s",in.BF_Vub_Uncertainties[0].Data()) <<endl;
    for(int m = 0; m < in.BF_Vub_Uncertainties.size(); m++) {
     cout << "Content is: "<< Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data()) <<endl;
     cout << "Pointer is: "<< ws->var( Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data()) )<<endl;
     ws->var( Form("Uncert_%s",in.BF_Vub_Uncertainties[m].Data()) )->setConstant(status);
    }
     cout << "Vub BF Status set"<<endl;
    for(int m = 0; m < in.BF_Vcb_Uncertainties.size(); m++) 
     ws->var( Form("Uncert_%s",in.BF_Vcb_Uncertainties[m].Data()) )->setConstant(status);
     cout << "Vcb BF Status set"<<endl;
    for(int m = 0; m < in.FF_Vub_Uncertainties.size(); m++) 
     ws->var( Form("Uncert_%s",in.FF_Vub_Uncertainties[m].Data()) )->setConstant(status);
      cout << "Vub FF Status set"<<endl;
    for(int m = 0; m < in.FF_Vcb_Uncertainties.size(); m++) 
    ws->var( Form("Uncert_%s",in.FF_Vcb_Uncertainties[m].Data()) )->setConstant(status);
     cout << "Vcb FF Status set"<<endl;
    for(int m = 0; m < in.general_Uncertainties.size(); m++) 
    ws->var( Form("Uncert_%s",in.general_Uncertainties[m].Data()) )->setConstant(status);
     cout << "Loop done"<<endl;
   }
   
}


