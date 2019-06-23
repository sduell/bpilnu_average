/*
 *  HiggsCombiner: Florian Bernlochner
 */

#include "HiggsDifferentialPlot.h"
#include "TPaveText.h"

HiggsPlot::HiggsPlot() { }

HiggsPlot::HiggsPlot(settings set, HiggsInput *input, HiggsWorkspace *workspace, TFile* fout) {
    
    Printf("> Initializing HiggsPlot Class\n");
    
    _set = set; _input = input; _workspace = workspace;
    _fout = fout;
    
    // We need to replace these
    SetAtlasStyle(); gStyle->SetPalette(1); gStyle->SetHistMinimumZero();
    
    _can = new TCanvas(); _ps= _set.PlotFileName;
    _can->Print(_ps+"[");
    
}

// ---------------------------------------------------------------------------------------

void HiggsPlot::Plot(RooFitResult *result, RooFormulaVar *nll, vector<TGraph*>* scans) {
    
    // Plot differential cross section and fit
    PlotDifferentialCrossSections(result);
    
    PlotComparison(result);
    
    if((_set.CombinationMethod == "full" || _set.CombinationMethod == "shape" || _set.CombinationMethod== "Multi") && _set.IncludeNP){
        PlotNuisanceParameters(result);
        PlotNuisanceParameterProjection(result,nll);
        //PlotNPScans(result,nll);
    }
    
    if ( _set.ManualScan=="true" ){
        PlotScans(result,nll,scans);
        //This was originally commented in!!!
        //PlotNuisanceScans(result,nll,scans);
    }
    else
        Projections(result,nll);
    
    // More control plots
    
    /*
     PlotMeasurements();
     PlotDifferentialMeasurements();
     PlotCorrelations();
     */
        
}

// ---------------------------------------------------------------------------------------
// Plot the differential cross sections and the combination result

void HiggsPlot::PlotDifferentialCrossSections(RooFitResult *result) {
    //VecI col; col.push_back( kAzure-3 ); col.push_back( kOrange+7 ); col.push_back( kGreen+2 ); col.push_back(kRed); col.push_back(kMagenta);
   VecI col; col.push_back( kAzure-3 ); col.push_back( kViolet-7 ); col.push_back( kGreen+2 ); col.push_back(kRed); col.push_back(kMagenta);
    VecI markerstyle; markerstyle.push_back(27); markerstyle.push_back(32);markerstyle.push_back(26); markerstyle.push_back(28); markerstyle.push_back(25); markerstyle.push_back(20);
    //VecI col; col.push_back( kAzure-3 ); col.push_back( kOrange+7 ); col.push_back(kRed); col.push_back(kMagenta);
    // Draw Input Data
    VecD dq2, dbhawon, dbsibidanov, dbpsibidanov, dblees, dbsanchez, dbfit;
    VecD dbhawonerr, dbsibidanoverr, dbpsibidanoverr, dbleeserr, dbsanchezerr, dbfiterr;
    dq2.push_back(0.);
    TH1F *axis;
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        Str xtit = _set.Unit == "NONE" ? _set.Label
        : _set.Label + "  [" + _set.Unit + "]";
        //Str ytit = Form("BF(B #rightarrow #pi l #nu)(%s)", _set.Label.Data());
        //Str ytit = Form("d#it{BF}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
        Str ytit = Form("d#it{B}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
        //Str ytit = _set.Unit == "NONE" ? Form("d#it{#sigma^{incl}} / d(#it{%s})  [pb]", _set.Label.Data() )
        //: Form("d#it{#sigma^{incl}} / d(#it{%s})  [pb/%s]", _set.Label.Data(), _set.Unit.Data() );
        
        if(_set.CombinationMethod == "shape") { ytit = " #frac{1}{#it{#sigma^{incl}}}  "+ytit; ytit.ReplaceAll("pb","1"); }
        
        axis = MakeAxis(in.Bins,xtit,ytit,0);
        TG *data = _set.CombinationMethod == "shape" ?
        DrawTG( UnfoldNormalized(in), UnfoldNormalizedError(in,i), in.Bins, col.size() > i ? col[i] : kBlack, true ) :
        DrawTG( Unfold(in), UnfoldError(in,i), in.Bins, col.size() > i ? col[i] : kBlack, true,0.5,1, markerstyle.size() > i ? markerstyle[i] : 1 ) ;

        TG *data_stat = DrawTG( Unfold(in), UnfoldStatError(in,i), in.Bins, col.size() > i ? col[i] : kBlack, true,0.5,1, markerstyle.size() > i ? markerstyle[i] : 1 ) ;
        TG *data_syst = DrawTG( Unfold(in), UnfoldSysError(in,i), in.Bins, col.size() > i ? col[i] : kBlack, true,0.5,1, markerstyle.size() > i ? markerstyle[i] : 1 ) ;


        for (int o=0;o<data->GetN();o++) {data->GetY()[o] *= 0.1; data->GetEY()[o] *=0.1;}
        TGA* asymdata=new TGA(data->GetN());
        for (int ox=0;ox<data->GetN();ox++) {asymdata->SetPoint(ox,data->GetX()[ox]+(-0.2+i*0.1), data->GetY()[ox]); asymdata->SetPointError(ox, (data->GetErrorX(ox)+(-0.2+i*0.1)),(data->GetErrorX(ox)-(-0.2+i*0.1)),(data->GetErrorY(ox)),(data->GetErrorY(ox)));}
        //asymdata->SetPointError(ox, (data->GetX()[ox]-data->GetErrorX(ox)),(data->GetX()[ox]+data->GetErrorX(ox)),(data->GetY()[ox]-data->GetErrorY(ox)),(data->GetY()[ox]+data->GetErrorY(ox)));}

        asymdata->SetMaximum(data->GetMaximum());
        asymdata->GetXaxis()->SetTitle(""); asymdata->GetYaxis()->SetTitle("");
        asymdata->SetLineColor(col[i]); asymdata->SetLineWidth(2);
        asymdata->SetMarkerStyle(markerstyle[i]); asymdata->SetMarkerSize(0.8); asymdata->SetMarkerColor(col[i]); asymdata->SetMinimum(0);
  

        //data->SetMarkerStyle(7); data->SetMarkerSize(1.5);
        axis->SetMaximum( 2.5 *data->GetMaximum()*0.1 );
        if(i == 0) axis->Draw();
        //data->Draw("PE+");
        //data->Draw("PE+"); data->Draw("[]");
        asymdata->Draw("PE+"); asymdata->Draw("[]");

        cout << "Get input Points from TGraph!"<<endl;
        VecD x, y, err, err_stat, err_syst;

        double sigma(0), error_tot_stat(0), error_tot_syst(0); 
        
        for(int i=0; i < data->GetN(); i++){
        double xpos, ypos;
        data->GetPoint(i,xpos,ypos);
        x.push_back(xpos);
        y.push_back(ypos);
        err.push_back(data->GetErrorY(i));
        err_stat.push_back(data_stat->GetErrorY(i));
        err_syst.push_back(data_syst->GetErrorY(i));
        cout << y[i] <<"\t"<<err_stat[i] << "\t" << err_syst[i] << endl;
        cout << y[i]*( in.Bins[i+1] - in.Bins[i] ) <<"\t +- "<<err_stat[i]*( in.Bins[i+1] - in.Bins[i] ) << "\t +- " << err_syst[i]*( in.Bins[i+1] - in.Bins[i] ) <<" ( "<<err_syst[i]/y[i]*10<<" % ) " << endl;
        sigma += ypos * ( in.Bins[i+1] - in.Bins[i] );
        //cout <<   ( in.Bins[i+1] - in.Bins[i] ) << endl;
        cout<<"Name is "<<in.Name.Data()<<endl;
        if(in.Name=="HaWon") {
            cout <<"Went into HaWon! "<<endl;
            dq2.push_back(in.Bins[i+1]);
            dbhawon.push_back(ypos *10.);
            dbhawonerr.push_back(err[i] * 10.);
        }
        if(in.Name=="Sibidanov") {
            dbsibidanov.push_back(ypos * 10.);
            dbsibidanoverr.push_back(err[i] * 10.);
        }
        if(in.Name=="SibidanovBp") {
            dbpsibidanov.push_back(ypos * 10.);
            dbpsibidanoverr.push_back(err[i] * 10.);
        }
        if(in.Name=="Lees") {
            dblees.push_back(ypos * 10.);
            dbleeserr.push_back(err[i] * 10.);
        }
        if(in.Name=="Sanchez") {
            dbsanchez.push_back(ypos * 10.);
            dbsanchezerr.push_back(err[i] * 10.);
        }


        }

                //Now read in the statistical correlation matrices
        ifstream infile;
        string namae;
        Str line;
        cout << "opening file "<< Form("/home/sduell/sources/bpilnu_new_average/CombinationCode/Correlations/%s_stat.dat", in.Name.Data())<<endl;
        infile.open(Form("/home/sduell/sources/bpilnu_new_average/CombinationCode/Correlations/%s_stat.dat", in.Name.Data()));
        vector<VecD> a;
        int n=0;
        cout << "begin reading file"<<endl;
        while(!infile.eof()){
            getline(infile,namae);
            line=Form("%s", namae.c_str());
            a.push_back(VecD());
            a[n]=VectorizeD(line," ");
            //cout << "line: "<<line<<endl;
            //cout << "element: "<<a[n][1]<<endl;
            n++;
        }
//Comment this line in to modify the cov matrix with a factor to account for the change in tauB0/tauB+
        //if(in.Name.Data()=="SibidanovBp") covfactor=1.0;
        cout << "file successfully read"<<endl;

        TMatrixDSym *CorMatrix = new TMatrixDSym(a[0].size());
        //normal Correlation matrix
        for(int j=0; j<a.size(); j++){
            for(int k=0; k<a[j].size(); k++){
                    (*CorMatrix)(j,k)=a[j][k];
                    //cout << "Matrix is: "<<a[j][k]<<endl;
            }   
        }


        cout << "Read in Correlation Matrix: "<<endl;
        (CorMatrix)->Print();

        //Now build the (statistical) covariance matrix from the correlation matrix
        VecD errs;
        
        for(int j = 0; j < in.NsigError.size(); j++)
            errs.push_back(in.NsigError[j]);
        
        TMatrixDSym *ErrMatrix = new TMatrixDSym(errs.size());
        
        for(int j=0; j<a.size(); j++){
            (*ErrMatrix)(j,j)=errs[j];
        }

        TMatrixDSym *CovMatrix = new TMatrixDSym(*CorMatrix);
        CovMatrix->SimilarityT(*ErrMatrix);
        cout << "Statistical covariance matrix built"<<endl;

        //CorMatrix->Print();
        //CovMatrix->Print();

        //Now calculate the statistical error of the BF!
        for(int o=0;o<a.size();o++){
            for(int p=0;p<a.size();p++){
                error_tot_stat+=(*CovMatrix)(o,p)*pow(0.1,2.);//* ( in.Bins[o+1] - in.Bins[o] )*( in.Bins[p+1] - in.Bins[p] );
            }
        }
        

        //Now get the systematic covariance matrix

        TMatrixDSym SystCovMatrix=(*_workspace->GetSyscovMatrix());
        SystCovMatrix.Print();
        int pos = _workspace->GetCovPos(i);
        int npoints = in.Nsig.size();
        for(int o=0;o<npoints;o++){
            //error_tot_syst+=SystCovMatrix(o+pos,o+pos)*pow(0.1,2.);
            for(int p=0;p<npoints;p++){
                //cout << "pos "<< pos<<endl;
                //cout << (( in.Bins[o+1] - in.Bins[o] )*( in.Bins[p+1] - in.Bins[p] )) << endl;
                error_tot_syst+=SystCovMatrix(o+pos,p+pos)*pow(0.1,2.);//(( in.Bins[o+1] - in.Bins[o] )*( in.Bins[p+1] - in.Bins[p] ));
            }
        }

        cout << " sigma_tot = " << sigma <<" +- "<<sqrt(error_tot_stat)<< " +- " << sqrt(error_tot_syst)<< endl;

        /*
        for(int i=0; i < data->GetN(); i++){
         for(int j=0; j < data->GetN(); j++){
        
         }
        }
        */


        TMarker *m;
        //DrawText("+  #scale[1.0]{#it{ "+_set.MeasLabel[i]+" }}",col.size() > i ? col[i] : kBlack,0.935-(i+2)*0.045,0.6); //was 0.875-(i+2)
        //if(i<3) DrawText("   #scale[1.0]{#it{ "+_set.MeasLabel[i]+" }}",col.size() > i ? col[i] : kBlack,0.935-(i+2)*0.045,0.2); //was 0.875-(i+2)
        //else DrawText("   #scale[1.0]{#it{ "+_set.MeasLabel[i]+" }}",col.size() > i ? col[i] : kBlack,0.935-(i+2-3)*0.045,0.565); //was 0.875-(i+2)
        
        
       // if(i<3) DrawText("   #scale[0.7]{#it{ "+_set.MeasLabel[i]+" }}",col.size() > i ? col[i] : kBlack,0.935-(i+2)*0.045,0.2); //was 0.875-(i+2)
        //else DrawText("   #scale[0.7]{#it{ "+_set.MeasLabel[i]+" }}",col.size() > i ? col[i] : kBlack,0.935-(i+2-3)*0.045,0.565); //was 0.875-(i+2)
        DrawText("   #scale[1.2]{#it{ "+_set.MeasLabel[i]+" }}",col.size() > i ? col[i] : kBlack,0.935-(i+2)*0.045,0.2); //was 0.875-(i+2)
        
       // if(i<3) m=new TMarker(2.675,12.175-i*0.74,markerstyle[i]);
       // else m=new TMarker(14.2,12.175-(i-3)*0.74,markerstyle[i]);
        m=new TMarker(2.675,17.-i*1.1,markerstyle[i]);
        m->SetMarkerColor(col[i]);
        m->Draw();
        
        _fout->cd();
        data->SetName(Form("%s_%s",_set.Variable.Data(),in.Name.Data()));
        data->Write();
        gROOT->cd();
    }


    //Draw the data point marker
    TMarker *m;
    //m=new TMarker(14.2,12.175-(5-3)*0.74,markerstyle[5]);
    m=new TMarker(2.69,17.-(5)*1.1,markerstyle[5]);
    m->SetMarkerColor(kBlack);
    m->Draw();


    // Now draw combination result
    cout << "Take values from here!"<<endl;
    VecD comb_sigma, comb_sigmaError; double tot_sigma, tot_sigmaError;
    for(int i = 0; i < _workspace->GetNMaxBins()-1; i++) {
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("sigma_%d",i+1)));
        comb_sigma.push_back( sigma->getVal() ); comb_sigmaError.push_back( sigma->getError() );
        tot_sigma += sigma->getVal(); tot_sigmaError = sqrt(pow(tot_sigmaError,2.) + pow(sigma->getError(),2.));
        cout << sigma->getVal()<< " pm "<< sigma->getError()<<endl;
    }
     
    
    ShapeErrors(comb_sigma,comb_sigmaError,result->correlationMatrix());
    
        
    TG *comb = _set.CombinationMethod == "shape" ?    
    DrawTG(Normalize(comb_sigma), ShapeErrors(comb_sigma,comb_sigmaError,result->correlationMatrix()),_workspace->GetMaxBins(),kBlack, true ) :
    DrawTG(comb_sigma,comb_sigmaError,_workspace->GetMaxBins(),kBlack, true, 0.5,1, 20 ) ;
    //DrawTG(comb_sigma,comb_sigmaError,_workspace->GetMaxBins(),kBlack, true ) ;
    
    for (int o=0;o<comb->GetN();o++) {comb->GetY()[o] *= 0.1; comb->GetEY()[o] *=0.1;}

    cout << "Get Points from TGraph!"<<endl;
    //comb->Draw("*PE"); comb->Draw("[]");
    comb->Draw("PE"); comb->Draw("[]");

    VecD x, y, err;
    for(int i=0; i<13; i++){
    double xpos, ypos;
    comb->GetPoint(i,xpos,ypos);
    x.push_back(xpos);
    y.push_back(ypos);
    err.push_back(comb->GetErrorY(i));
    cout << y[i] <<"\t"<<err[i]<<endl;
    dbfit.push_back(ypos*10);
    dbfiterr.push_back(err[i]*10.);
    }

    cout<<"The Fit Result is: "<<endl;
    result->Print();

    //Print out a table for Latex
    //for(int i=0; i<13; i++){
    //    if(i==0 || i==2 || i==4 || i==6 || i==8 || i==10) Printf("$ %f - %f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\", dq2[i], dq2[i+1], dbhawon[i], dbhawonerr[i], dbsibidanov[i], dbsibidanoverr[i], dbpsibidanov[i], dbpsibidanoverr[i], dblees[i], dbleeserr[i], dbsanchez[i], dbsanchezerr[i], dbfit[i], dbfiterr[i]);
    //}
    cout <<"plotting latex table"<<endl;
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[0], dq2[1], dbhawon[0], dbhawonerr[0], dbsibidanov[0], dbsibidanoverr[0], dbpsibidanov[0], dbpsibidanoverr[0], dblees[0], dbleeserr[0], dbsanchez[0], dbsanchezerr[0], dbfit[0], dbfiterr[0]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[1], dq2[2], dbhawon[1], dbhawonerr[1], dbsibidanov[1], dbsibidanoverr[1], dbpsibidanov[0], dbpsibidanoverr[0], dblees[1], dbleeserr[1], dbsanchez[0], dbsanchezerr[0], dbfit[1], dbfiterr[1]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[2], dq2[3], dbhawon[2], dbhawonerr[2], dbsibidanov[2], dbsibidanoverr[2], dbpsibidanov[1], dbpsibidanoverr[1], dblees[2], dbleeserr[2], dbsanchez[1], dbsanchezerr[1], dbfit[2], dbfiterr[2]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[3], dq2[4], dbhawon[3], dbhawonerr[3], dbsibidanov[3], dbsibidanoverr[3], dbpsibidanov[1], dbpsibidanoverr[1], dblees[3], dbleeserr[3], dbsanchez[1], dbsanchezerr[1], dbfit[3], dbfiterr[3]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[4], dq2[5], dbhawon[4], dbhawonerr[4], dbsibidanov[4], dbsibidanoverr[4], dbpsibidanov[2], dbpsibidanoverr[2], dblees[4], dbleeserr[4], dbsanchez[2], dbsanchezerr[2], dbfit[4], dbfiterr[4]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[5], dq2[6], dbhawon[5], dbhawonerr[5], dbsibidanov[5], dbsibidanoverr[5], dbpsibidanov[2], dbpsibidanoverr[2], dblees[5], dbleeserr[5], dbsanchez[2], dbsanchezerr[2], dbfit[5], dbfiterr[5]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[6], dq2[7], dbhawon[6], dbhawonerr[6], dbsibidanov[6], dbsibidanoverr[6], dbpsibidanov[3], dbpsibidanoverr[3], dblees[6], dbleeserr[6], dbsanchez[3], dbsanchezerr[3], dbfit[6], dbfiterr[6]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[7], dq2[8], dbhawon[7], dbhawonerr[7], dbsibidanov[7], dbsibidanoverr[7], dbpsibidanov[3], dbpsibidanoverr[3], dblees[7], dbleeserr[7], dbsanchez[3], dbsanchezerr[3], dbfit[7], dbfiterr[7]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[8], dq2[9], dbhawon[8], dbhawonerr[8], dbsibidanov[8], dbsibidanoverr[8], dbpsibidanov[4], dbpsibidanoverr[4], dblees[8], dbleeserr[8], dbsanchez[4], dbsanchezerr[4], dbfit[8], dbfiterr[8]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[9], dq2[10], dbhawon[9], dbhawonerr[9], dbsibidanov[9], dbsibidanoverr[9], dbpsibidanov[4], dbpsibidanoverr[4], dblees[9], dbleeserr[9], dbsanchez[4], dbsanchezerr[4], dbfit[9], dbfiterr[9]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[10], dq2[11], dbhawon[10], dbhawonerr[10], dbsibidanov[10], dbsibidanoverr[10], dbpsibidanov[5], dbpsibidanoverr[5], dblees[10], dbleeserr[10], dbsanchez[5], dbsanchezerr[5], dbfit[10], dbfiterr[10]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[11], dq2[12], dbhawon[11], dbhawonerr[11], dbsibidanov[11], dbsibidanoverr[11], dbpsibidanov[5], dbpsibidanoverr[5], dblees[11], dbleeserr[11], dbsanchez[5], dbsanchezerr[5], dbfit[11], dbfiterr[11]);
    Printf("$ %1.0f - %1.0f $ & $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $ &  $ %2.1f \\pm %2.1f $\\\\", dq2[12], dq2[13], dbhawon[12], dbhawonerr[12], dbsibidanov[12], dbsibidanoverr[12], dbpsibidanov[6], dbpsibidanoverr[6], dblees[11], dbleeserr[11], dbsanchez[5], dbsanchezerr[5], dbfit[12], dbfiterr[12]);
 

    //const TMatrixDSym &FitCovMatrix = _result->covarianceMatrix();

    //FitCovMatrix.Print();

    //ofstream covout;
    //covout.open("/home/sduell/bpilnu_new_average/CombinationCode/results/fit_covmat.dat", ios::out);

    //this is wrong, only if the total cov matrix is 47x47
    //vector <int> irow {35,40,41,42,43,44,45,46,47,36,37,38,39};
    //vector <int> jrow {35,40,41,42,43,44,45,46,47,36,37,38,39};

    //this is for 52x52
    //vector <int> irow {38,43,44,45,46,47,48,49,50,39,40,41,42};
    //vector <int> jrow {38,43,44,45,46,47,48,49,50,39,40,41,42};

    //this is for a 57x57 matrix (B0 B+ combination)
    //vector <int> irow {43,48,49,50,51,52,53,54,55,44,45,46,47};
    //vector <int> jrow {43,48,49,50,51,52,53,54,55,44,45,46,47};

    //this is for a 55x55 matrix (B0 B+ combination)
    vector <int> irow {41,46,47,48,49,50,51,52,53,42,43,44,45};
    vector <int> jrow {41,46,47,48,49,50,51,52,53,42,43,44,45};

    //this is for when using totcormat
    //vector <int> irow {9,14,15,16,17,18,19,20,21,10,11,12,13};
    //vector <int> jrow {9,14,15,16,17,18,19,20,21,10,11,12,13};

    const TMatrixDSym &FitCorMatrix = result->correlationMatrix();
    TMatrixDSym ReducedFitCorMatrix(13);

    cout << "FitCorMatrix is: "<<endl;
    FitCorMatrix.Print();

    ofstream corout;
    corout.open("/home/sduell/bpilnu_new_average/CombinationCode/results/fit_cormat.dat", ios::out);

    for(int i=0; i<13/*FitCorMatrix.GetNrows()*/;i++){
        for(int j=0; j<13/*FitCorMatrix.GetNrows()*/; j++){
            //corout <<(double)(FitCorMatrix)(i,j)<<" ";
            corout <<(double)(FitCorMatrix)(irow[i]+1,jrow[j]+1)<<" ";
            (ReducedFitCorMatrix)(i,j)=(FitCorMatrix)(irow[i]+1,jrow[j]+1);
        }
        corout <<""<<endl;
    }

    corout.close();

    //Calculate Covariance matrix
    TMatrixDSym CMatrix = ReducedFitCorMatrix; //Watch out, this will overwrite the correlation matrix
    //create error matrix
    TMatrixDSym ErrorMatrix(13);
    for(int i=0;i<13;i++){
        (ErrorMatrix)(i,i) = err[i];
    }
    cout << "ErrorMatrix "<<endl;
    ErrorMatrix.Print();
    CMatrix.SimilarityT(ErrorMatrix);
    cout << "Covariance Matrix is "<<endl;
    CMatrix.Print();

    //ofstream covout;

    corout.open("/home/sduell/bpilnu_new_average/CombinationCode/results/fit_covmat.dat", ios::out);

    for(int i=0; i<13/*FitCovMatrix.GetNrows()*/;i++){
        for(int j=0; j<13/*FitCovMatrix.GetNrows()*/; j++){
            corout <<(double)(CMatrix)(i,j)<<" ";
            cout<< "Matrixelement should be "<<(double)(CMatrix)(i,j)<<endl;
        }
    }
    corout.close();

    //DrawText("Input Measurements: ",kBlack,0.935-(1)*0.045,0.69);
    //DrawText("#scale[1.0]{*  Likelihood fit average}",kBlack,0.935-(_input->GetInput().size()+2)*0.045,0.69);
    //DrawText("Input Measurements: ",kBlack,0.935-(1)*0.045,0.6);
    DrawText("Input Measurements: ",kBlack,0.935-(1)*0.045,0.2);
    
    //DrawText("#scale[1.0]{*  Likelihood fit average}",kBlack,0.935-(_input->GetInput().size()+2)*0.045,0.6);
    //DrawText("#scale[1.0]{*   Likelihood fit average}",kBlack,0.935-(_input->GetInput().size()+2-3)*0.045,0.564);
    //DrawText("#scale[0.7]{     Likelihood fit average}",kBlack,0.935-(_input->GetInput().size()+2-3)*0.045,0.564);
    
    
    //DrawText("#scale[0.7]{      Likelihood fit average}",kBlack,0.934-(_input->GetInput().size()+2-3)*0.045,0.565);
    DrawText("#scale[1.2]{      Likelihood fit average}",kBlack,0.93-(_input->GetInput().size()+2)*0.045,0.185);
    
    axis->Draw("same");

    //this is the HFLAV box 1.25,1.5
    Double_t xpos=0.175;
    Double_t ypos=0.175;
    Double_t scale=1;
    const char * labelChar = "2018";
    TString label(labelChar);
    TVirtualPad* thePad;
  
    if ((thePad = TVirtualPad::Pad()) == 0) return;

    UInt_t pad_width(thePad->XtoPixel(thePad->GetX2()));
    UInt_t pad_height(thePad->YtoPixel(thePad->GetY1()));
  
    Double_t ysiz_pixel(25);
    Double_t ysiz(Double_t(ysiz_pixel)/Double_t(pad_height));
    Double_t xsiz(4.8*ysiz*Double_t(pad_height)/Double_t(pad_width));

    Double_t x1, x2, y1, y2;
    xsiz = scale*xsiz;
    ysiz = scale*ysiz;

    if (xpos >= 0) {
        x1 = xpos;
        x2 = xpos + xsiz;
    } else {
        x1 = 1 + xpos - xsiz;
        x2 = 1 + xpos;
    }

    if (ypos >= 0) {
        y1 = ypos+0.9*ysiz;
        y2 = ypos+0.9*ysiz + ysiz;
    } else {
        y1 = 1 + ypos - ysiz;
        y2 = 1 + ypos;
    }

    TPaveText *tbox1 = new TPaveText(x1, y1, x2, y2, "BRNDC");
    // tbox1->SetLineColor(1);
    // tbox1->SetLineStyle(1);
    // tbox1->SetLineWidth(2);
    tbox1->SetFillColor(kBlack);
    tbox1->SetFillStyle(1001);
    // tbox1->SetBorderSize(1);
    tbox1->SetShadowColor(kWhite);
    tbox1->SetTextColor(kWhite);
    tbox1->SetTextFont(76);
    tbox1->SetTextSize(24*scale);
    tbox1->SetTextAlign(22); //center-adjusted and vertically centered
    tbox1->AddText(TString("HFLAV"));
    tbox1->Draw();
    //
    TPaveText *tbox2 = new TPaveText(x1, y1-0.9*ysiz, x2, y2-ysiz, "BRNDC");
    // tbox2->SetLineColor(1);
    // tbox2->SetLineStyle(1);
    // tbox2->SetLineWidth(2);
    tbox2->SetFillColor(kWhite);
    tbox2->SetFillStyle(1001);
    // tbox2->SetBorderSize(1);
    tbox2->SetShadowColor(kWhite);
    tbox2->SetTextColor(kBlack);
    tbox2->SetTextFont(76);
    tbox2->SetTextSize(18*scale);
    tbox2->SetTextAlign(22); //center-adjusted and vertically centered
    tbox2->AddText(label);
    tbox2->Draw();


    //this is the HFAG box old code
/*   //TBox *box = new TBox(0.1585227,-8.476744,0.2011364,-7.534884);
   TBox *box = new TBox(1.25,1,5.75,2.5);
   box->SetFillColor(19);
   box->SetFillStyle(0);
   box->SetLineWidth(2);
   box->Draw();
   //box = new TBox(0.1585227,-8.162791,0.2011364,-7.534884);
   box = new TBox(1.25,1.5,5.75,2.5);
   box->SetFillColor(1);
   box->SetFillStyle(1000);
   box->SetLineWidth(2);
   box->Draw();
   //TText *text = new TText(0.1798295,-7.848837,"HFAG");
   //TText *text = new TText(3.5,2.05,"HFAG");
   TText *text = new TText(3.5,2.05,"HFLAV");
   text->SetTextAlign(22);
   text->SetTextColor(0);
   text->SetTextFont(53);
   text->SetTextSize(30);
   text->Draw();
   text = new TText(3.5,1.3,"Summer 2019");
   text->SetTextAlign(22);
   text->SetTextFont(53);
   text->SetTextSize(15);
   text->Draw();
*/
    _fout->cd();
    comb->SetName(Form("%s_comb",_set.Variable.Data()));
    comb->Write();
    gROOT->cd();
    
    _can->Print("average.eps");
    _can->Print("average.pdf");
    _can->Print("average.gif");
    _can->Print(_ps); _can->Print(Form("plots/%s_diffxsec.pdf",_set.Variable.Data()));
}

// ---------------------------------------------------------------------------------------
// Compare with Dag's results
void HiggsPlot::PlotComparison(RooFitResult *result) {
    
    VecI col; col.push_back( kAzure-3 ); col.push_back( kOrange+7 ); col.push_back(kRed); col.push_back(kMagenta);
    // Draw Input Data
    TH1F *axis;
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        Str xtit = _set.Unit == "NONE" ? _set.Label
        : _set.Label + "  [" + _set.Unit + "]";
        //Str ytit = Form("BF(B #rightarrow #pi l #nu)(%s)", _set.Label.Data());
        Str ytit = Form("d#it{B}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
        //Str ytit = _set.Unit == "NONE" ? Form("d#it{#sigma^{incl}} / d(#it{%s})  [pb]", _set.Label.Data() )
        //: Form("d#it{#sigma^{incl}} / d(#it{%s})  [pb/%s]", _set.Label.Data(), _set.Unit.Data() );
        
        if(_set.CombinationMethod == "shape") { ytit = " #frac{1}{#it{#sigma^{incl}}}  "+ytit; ytit.ReplaceAll("pb","1"); }
        
        TG *data = _set.CombinationMethod == "shape" ?
        DrawTG( UnfoldNormalized(in), UnfoldNormalizedError(in,i), in.Bins, col.size() > i ? col[i] : kBlack, true ) :
        DrawTG( Unfold(in), UnfoldError(in,i), in.Bins, col.size() > i ? col[i] : kBlack, true ) ;
        
        axis = MakeAxis(in.Bins,xtit,ytit,0);
        axis->SetMaximum( 1.8 *data->GetMaximum() );
        if(i == 0) axis->Draw();
    }
    // Now draw combination result
    VecD comb_sigma, comb_sigmaError; double tot_sigma, tot_sigmaError;
    for(int i = 0; i < _workspace->GetNMaxBins()-1; i++) {
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("sigma_%d",i+1)));
        comb_sigma.push_back( sigma->getVal() ); comb_sigmaError.push_back( sigma->getError() );
        tot_sigma += sigma->getVal(); tot_sigmaError = sqrt(pow(tot_sigmaError,2.) + pow(sigma->getError(),2.));
    }
    TG *comb = _set.CombinationMethod == "shape" ?
    DrawTG(Normalize(comb_sigma), ShapeErrors(comb_sigma,comb_sigmaError,result->correlationMatrix()),_workspace->GetMaxBins(),kBlack, true ) :
    DrawTG(comb_sigma,comb_sigmaError,_workspace->GetMaxBins(),kBlack, true ) ;
    comb->Draw("*PE"); comb->Draw("[]");
    DrawText("#scale[1.0]{*  Likelihood fit average}",kBlack,0.875-(_input->GetInput().size()+2)*0.045,0.69);
    //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
    //DrawText("#scale[0.8]{#int} #it{L} d#it{t} = 20.3/fb, #sqrt{s} = 8 TeV",kBlack,(0.93-1.*0.035),0.69);
    //DrawText(Form("#sigma_{incl} = %4.2f #pm %4.2f pb",tot_sigma,tot_sigmaError),kBlack,0.875-(_input->GetInput().size()+3.3)*0.045,0.69);
    // Now draw chi2 combination result if specified

    axis->Draw("same");
    _can->Print(_ps); _can->Print(Form("plots/%s_compxsec.pdf",_set.Variable.Data()));
    
}

// ---------------------------------------------------------------------------------------
// Plot Nuisance Parameters
void HiggsPlot::PlotNuisanceParameters(RooFitResult *result) {
    
    // Collect the nuisance parameters from the workspace
    VecD NuisPar, NuisParError; StrV NuisParName; std::map<Str,bool> pull_plotted;
    // Fetch Nuisance Parameters from workspace
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];


        for(int j = 0; j < in.BF_Vub_Uncertainties.size(); j++) {
            RooRealVar *NP = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",in.BF_Vub_Uncertainties[j].Data())));
            if(!NP) Fatal( Form("Could not find Uncert_%s!",in.BF_Vub_Uncertainties[j].Data()), " ");
            if(pull_plotted[Form("Uncert_%s",in.BF_Vub_Uncertainties[j].Data())]) continue;
            NuisPar.push_back( NP->getVal() ); NuisParError.push_back( NP->getError() ); NuisParName.push_back( Form("Uncert_%s",in.BF_Vub_Uncertainties[j].Data()) );
            pull_plotted[Form("Uncert_%s",in.BF_Vub_Uncertainties[j].Data())] = true;
        }
        for(int j = 0; j < in.BF_Vcb_Uncertainties.size(); j++) {
            RooRealVar *NP = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",in.BF_Vcb_Uncertainties[j].Data())));
            if(!NP) Fatal( Form("Could not find Uncert_%s!",in.BF_Vcb_Uncertainties[j].Data()), " ");
            if(pull_plotted[Form("Uncert_%s",in.BF_Vcb_Uncertainties[j].Data())]) continue;
            NuisPar.push_back( NP->getVal() ); NuisParError.push_back( NP->getError() ); NuisParName.push_back( Form("Uncert_%s",in.BF_Vcb_Uncertainties[j].Data()) );
            pull_plotted[Form("Uncert_%s",in.BF_Vcb_Uncertainties[j].Data())] = true;
        }
        for(int j = 0; j < in.FF_Vub_Uncertainties.size(); j++) {
            RooRealVar *NP = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",in.FF_Vub_Uncertainties[j].Data())));
            if(!NP) Fatal( Form("Could not find Uncert_%s!",in.FF_Vub_Uncertainties[j].Data()), " ");
            if(pull_plotted[Form("Uncert_%s",in.FF_Vub_Uncertainties[j].Data())]) continue;
            NuisPar.push_back( NP->getVal() ); NuisParError.push_back( NP->getError() ); NuisParName.push_back( Form("Uncert_%s",in.FF_Vub_Uncertainties[j].Data()) );
            pull_plotted[Form("Uncert_%s",in.FF_Vub_Uncertainties[j].Data())] = true;
        }
        for(int j = 0; j < in.FF_Vcb_Uncertainties.size(); j++) {
            RooRealVar *NP = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",in.FF_Vcb_Uncertainties[j].Data())));
            if(!NP) Fatal( Form("Could not find Uncert_%s!",in.FF_Vcb_Uncertainties[j].Data()), " ");
            if(pull_plotted[Form("Uncert_%s",in.FF_Vcb_Uncertainties[j].Data())]) continue;
            NuisPar.push_back( NP->getVal() ); NuisParError.push_back( NP->getError() ); NuisParName.push_back( Form("Uncert_%s",in.FF_Vcb_Uncertainties[j].Data()) );
            pull_plotted[Form("Uncert_%s",in.FF_Vcb_Uncertainties[j].Data())] = true;
        }
        for(int j = 0; j < in.general_Uncertainties.size(); j++) {
            RooRealVar *NP = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",in.general_Uncertainties[j].Data())));
            if(!NP) Fatal( Form("Could not find Uncert_%s!",in.general_Uncertainties[j].Data()), " ");
            if(pull_plotted[Form("Uncert_%s",in.general_Uncertainties[j].Data())]) continue;
            NuisPar.push_back( NP->getVal() ); NuisParError.push_back( NP->getError() ); NuisParName.push_back( Form("Uncert_%s",in.general_Uncertainties[j].Data()) );
            pull_plotted[Form("Uncert_%s",in.general_Uncertainties[j].Data())] = true;
        }



    }
    VecD bins, binlims; for(int i = 0; i < NuisPar.size(); i++) bins.push_back(i+1.5); binlims.push_back(1.5); binlims.push_back(NuisPar.size()+1.5);
    TH2D *axis = new TH2D("axis","axis;Nuisance Parameter Value / Error",7,-3,3,NuisParName.size()+2,0,NuisParName.size()+2);
    axis->SetMinimum(-1); axis->SetMaximum(NuisParName.size()-1);
    for(int i = NuisParName.size()-1; i >= 0; i--)
        axis->GetYaxis()->SetBinLabel(i+2,"#scale[0.5]{"+NuisParName[i]+"}");
    axis->Draw();
    TGraphErrors* tg_nuis = Draw1D(NuisPar,NuisParError,"",NuisPar.size(),bins,binlims, kAzure-3,1,"Nuisance Parameter Value / Error","Nuisance Parameter",false);
    tg_nuis->SetMarkerStyle(25); tg_nuis->SetMarkerColor(kAzure-3); tg_nuis->SetLineColor(12);
    tg_nuis->Draw("P"); tg_nuis->Print("all");
    //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
    TArrow sig_min(-1,1,-1,NuisPar.size()+1.,0.05,""); sig_min.SetLineColor(kAzure-3); sig_min.SetLineStyle(7); sig_min.Draw("same");
    TArrow sig_max(1,1,1,NuisPar.size()+1.,0.05,""); sig_max.SetLineColor(kAzure-3); sig_max.SetLineStyle(7); sig_max.Draw("same");
    _can->Print(_ps); _can->Print(Form("plots/%s_NP.pdf",_set.Variable.Data()));
    
    // Now collect the background nuisance parameters
    VecD BkgNuisPar, BkgNuisParError; StrV BkgNuisParName;
    // Fetch Nuisance Parameters from workspace
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        for(int j = 0; j < in.Bkg_Uncertainties.size(); j++) {
            RooRealVar *NP = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",in.Bkg_Uncertainties[j].Data())));
            if(pull_plotted[Form("Uncert_%s",in.Bkg_Uncertainties[j].Data())]) continue;
            if(!NP) Fatal( Form("Could not find Uncert_%s!",in.Bkg_Uncertainties[j].Data()), " ");
            BkgNuisPar.push_back( NP->getVal() ); BkgNuisParError.push_back( NP->getError() ); BkgNuisParName.push_back( Form("Uncert_%s",in.Bkg_Uncertainties[j].Data()) );
            pull_plotted[Form("Uncert_%s",in.Bkg_Uncertainties[j].Data())] = true;
        }
    }
    bins.clear(); binlims.clear(); for(int i = 0; i < BkgNuisPar.size(); i++) bins.push_back(i+1.5); binlims.push_back(1.5); binlims.push_back(BkgNuisPar.size()+1.5);
    TH2D *bkgaxis = new TH2D("axis","axis; Bkg Nuisance Parameter Value / Error",7,-3,3,BkgNuisParName.size()+2,0,BkgNuisParName.size()+2);
    bkgaxis->SetMinimum(-1); bkgaxis->SetMaximum(BkgNuisParName.size()-1);
    for(int i = BkgNuisParName.size()-1; i >= 0; i--)
        bkgaxis->GetYaxis()->SetBinLabel(i+2,"#scale[0.5]{"+BkgNuisParName[i]+"}");
    bkgaxis->Draw();
    TGraphErrors* tg_bkgnuis = Draw1D(BkgNuisPar,BkgNuisParError,"",BkgNuisPar.size(),bins,binlims, kAzure-3,1,"Bkg Nuisance Parameter Value / Error","Nuisance Parameter",false);
    tg_bkgnuis->SetMarkerStyle(25); tg_bkgnuis->SetMarkerColor(kAzure-3); tg_bkgnuis->SetLineColor(12);
    tg_bkgnuis->Draw("P"); tg_bkgnuis->Print("all");
    //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
    TArrow bkgsig_min(-1,1,-1,BkgNuisPar.size()+1.,0.05,""); bkgsig_min.SetLineColor(kAzure-3); bkgsig_min.SetLineStyle(7); bkgsig_min.Draw("same");
    TArrow bkgsig_max(1,1,1,BkgNuisPar.size()+1.,0.05,""); bkgsig_max.SetLineColor(kAzure-3); bkgsig_max.SetLineStyle(7); bkgsig_max.Draw("same");
    _can->Print(_ps); _can->Print(Form("plots/%s_BkgNP.pdf",_set.Variable.Data()));
   
}

// ---------------------------------------------------------------------------------------
// Control plot functions

void HiggsPlot::Projections(RooFitResult *result, RooFormulaVar *nll) {
    
    for(int i = 0; i < _workspace->GetNMaxBins()-1; i++) {
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("sigma_%d",i+1)));
        double range_low = sigma->getVal() - 3. * sigma->getError() > 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
        double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
        RooPlot* frame = sigma->frame(Name(Form("nll_frame_%d",i+1)),Range(range_low,range_high),Title(Form("-log(L) scan vs #sigma_{%d}",i+1)));
        frame->SetXTitle(Form("#sigma_{%d}  [pb]",i+1)); frame->SetYTitle("#it{Projection of '-2 log(L)'}");
        nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        frame->SetMinimum(0.);
        frame->Draw();
        TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kBlack); error_p68.SetLineWidth(2); error_p68.SetLineStyle(7); error_p68.Draw("same");
        TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(7); error_p95.Draw("same");
        DrawText(Form("#scale[1.5]{#sigma_{%d} = %4.2f #pm %4.2f pb}",i+1,sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
        //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
        _can->Print(_ps); _can->Print(Form("plots/%s_sigma_%d.pdf",_set.Variable.Data(),i+1));
    }
    
}

// ---------------------------------------------------------------------------------------
// Control plot functions


void HiggsPlot::PlotNuisanceParameterProjection_Dummy(RooFitResult *result, RooFormulaVar *nll) {
        
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_Dummy")));
        //RooAbsReal* pll=nll->createProfile(*sigma);
        double range_low = sigma->getVal() - 3. * sigma->getError() < 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
        //double range_low = sigma->getVal() - 5. * sigma->getError() > 0 ? sigma->getVal() - 5. * sigma->getError() : 0. , range_high = sigma->getVal() + 5. * sigma->getError();
        double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
        RooPlot* frame = sigma->frame(Name(Form("nll_frame")),Range(range_low,range_high),Title(Form("-log(L) scan vs Dummy")));
        frame->SetXTitle(Form("Dummy")); frame->SetYTitle("#it{Projection of '-2 log(L)'}");
        //nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        frame->SetMinimum(0.);
        frame->Draw();
        TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kBlack); error_p68.SetLineWidth(2); error_p68.SetLineStyle(7); error_p68.Draw("same");
        TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(7); error_p95.Draw("same");
        DrawText(Form("#scale[1.5]{Dummy = %4.4f #pm %4.4f}",sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
        //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
        _can->Print(_ps); _can->Print(Form("plots/%s_Dummy.pdf",_set.Variable.Data()));   

}

void HiggsPlot::PlotNuisanceParameterProjection(RooFitResult *result, RooFormulaVar *nll) {
        cout << "Begin plotting Nuisance Parameter Likelihood Profiles"<<endl;
        vector<Str> npname;
        for(int j = 0; j < _input->GetInput().size(); j++) {
            input in = _input->GetInput()[j];
            //cout << "got input"<<endl;

            //cout << "got input 3"<<endl;
            for(int i=0; i<in.Bkg_Uncertainties.size(); i++){
                npname.push_back(in.Bkg_Uncertainties[i].Data());
            }
            //cout << "got input 4"<<endl;
            for(int i=0; i<in.BF_Vub_Uncertainties.size(); i++){
                //cout << in.BF_Vub_Uncertainties[i].Data()<<endl;
                npname.push_back(in.BF_Vub_Uncertainties[i]);
            }
            //cout << "got input 5"<<endl;
            for(int i=0; i<in.BF_Vcb_Uncertainties.size(); i++){
                npname.push_back(in.BF_Vcb_Uncertainties[i].Data());
            }
            //cout << "got input 6"<<endl;
            for(int i=0; i<in.FF_Vub_Uncertainties.size(); i++){
                npname.push_back(in.FF_Vub_Uncertainties[i].Data());
            }
            //cout << "got input 7"<<endl;
            for(int i=0; i<in.FF_Vcb_Uncertainties.size(); i++){
                npname.push_back(in.FF_Vcb_Uncertainties[i].Data());
            }
            for(int i=0; i<in.general_Uncertainties.size(); i++){
                npname.push_back(in.general_Uncertainties[i].Data());
            }
        }
        cout << "Nuisance Parameters read in: "<<endl;
        for(int j=0; j<npname.size(); j++){
            cout << npname[j].Data()<<endl;
        }
        for(int i = 0; i < npname.size(); i++) {
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",npname[i].Data())));
        //RooAbsReal* pll=nll->createProfile(*sigma);
        double range_low = sigma->getVal() - 3. * sigma->getError() < 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
        //double range_low = sigma->getVal() - 5. * sigma->getError() > 0 ? sigma->getVal() - 5. * sigma->getError() : 0. , range_high = sigma->getVal() + 5. * sigma->getError();
        double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
        RooPlot* frame = sigma->frame(Name(Form("nll_frame")),Range(range_low,range_high),Title(Form("-log(L) scan vs %s", npname[i].Data())));
        frame->SetXTitle(Form("%s", npname[i].Data())); frame->SetYTitle("#it{Projection of '-2 log(L)'}");
        //nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        frame->SetMinimum(0.);
        frame->Draw();
        TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kBlack); error_p68.SetLineWidth(2); error_p68.SetLineStyle(7); error_p68.Draw("same");
        TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(7); error_p95.Draw("same");
        DrawText(Form("#scale[1.5]{%s = %4.4f #pm %4.4f}",npname[i].Data(),sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
        //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
        _can->Print(_ps); _can->Print(Form("plots/%s_%s.pdf",_set.Variable.Data(), npname[i].Data()));   
    }
}

// ---------------------------------------------------------------------------------------
// Control plot functions


void HiggsPlot::PlotScans(RooFitResult *result, RooFormulaVar *nll, vector<TGraph*>* scans) {
    
    for(int i = 0; i < _workspace->GetNMaxBins()-1; i++) {
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("sigma_%d",i+1)));
        double range_low = sigma->getVal() - 3. * sigma->getError() > 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
        double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
        RooPlot* frame = sigma->frame(Name(Form("nll_frame_%d",i+1)),Range(range_low,range_high),Title(Form("-log(L) scan vs #sigma_{%d}",i+1)));
        frame->SetXTitle(Form("#sigma_{%d}  [pb]",i+1)); frame->SetYTitle("#it{Scan of '-2 log(L)'}");
        nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        frame->SetMinimum(0.);
        frame->Draw();
        TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kBlack); error_p68.SetLineWidth(2); error_p68.SetLineStyle(7); error_p68.Draw("same");
        TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(7); error_p95.Draw("same");
        DrawText(Form("#scale[1.5]{#sigma_{%d} = %4.2f #pm %4.2f pb}",i+1,sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
        //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
        
        scans->at(i)->SetLineWidth(3);
        scans->at(i)->SetLineColor(kOrange+1);
        scans->at(i)->Draw("L");
        DrawText("#scale[1.2]{Projection}",kRed,0.9,0.75);
        DrawText("#scale[1.2]{Manual scan}",kOrange+1,0.85,0.75);
        
        double pos_err = DetermineQuantile(scans->at(i), 1, sigma->getVal(), sigma->getVal() + 3.*sigma->getError()) - sigma->getVal();
        double neg_err = -1*(DetermineQuantile(scans->at(i), 1, sigma->getVal(), sigma->getVal() - 3.*sigma->getError()) - sigma->getVal());
        
        DrawText(Form("#scale[1]{Error from manual scan = + %4.2f - %4.2f pb}", pos_err, neg_err),kBlack,0.675-(2)*0.045,0.45);
        
        
        _can->Print(_ps); _can->Print(Form("plots/%s_sigma_%d.pdf",_set.Variable.Data(),i+1));
    }
    
}

// ---------------------------------------------------------------------------------------
// Control plot functions


void HiggsPlot::PlotNuisanceScans(RooFitResult *result, RooFormulaVar *nll, vector<TGraph*>* scans) {
        vector<Str> npname;
        for(int j = 0; j < _input->GetInput().size(); j++) {
            input in = _input->GetInput()[j];
            //cout << "got input"<<endl;
            //cout << "got input 3"<<endl;
            for(int i=0; i<in.Bkg_Uncertainties.size(); i++){
                npname.push_back(in.Bkg_Uncertainties[i].Data());
            }
            //cout << "got input 4"<<endl;
            for(int i=0; i<in.BF_Vub_Uncertainties.size(); i++){
                //cout << in.BF_Vub_Uncertainties[i].Data()<<endl;
                npname.push_back(in.BF_Vub_Uncertainties[i]);
            }
            //cout << "got input 5"<<endl;
            for(int i=0; i<in.BF_Vcb_Uncertainties.size(); i++){
                npname.push_back(in.BF_Vcb_Uncertainties[i].Data());
            }
            //cout << "got input 6"<<endl;
            for(int i=0; i<in.FF_Vub_Uncertainties.size(); i++){
                npname.push_back(in.FF_Vub_Uncertainties[i].Data());
            }
            //cout << "got input 7"<<endl;
            for(int i=0; i<in.FF_Vcb_Uncertainties.size(); i++){
                npname.push_back(in.FF_Vcb_Uncertainties[i].Data());
            }
            for(int i=0; i<in.general_Uncertainties.size(); i++){
                npname.push_back(in.general_Uncertainties[i].Data());
            }
        }
        cout << "Nuisance Parameters read in: "<<endl;
        for(int j=0; j<npname.size(); j++){
            cout << npname[j].Data()<<endl;
        }
        for(int i = 0; i < npname.size(); i++) {
        RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_%s",npname[i].Data())));
        double range_low = sigma->getVal() - 3. * sigma->getError() < 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
        double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
        RooPlot* frame = sigma->frame(Name(Form("Uncert_%s",npname[i].Data())),Range(range_low,range_high),Title(Form("-log(L) scan vs %s", npname[i].Data())));
        frame->SetXTitle(Form("%s",npname[i].Data())); frame->SetYTitle("#it{Scan of '-2 log(L)'}");
        nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        frame->SetMinimum(0.);
        frame->Draw();
        TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kBlack); error_p68.SetLineWidth(2); error_p68.SetLineStyle(7); error_p68.Draw("same");
        TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(7); error_p95.Draw("same");
        DrawText(Form("#scale[1.5]{%s = %4.4f #pm %4.4f pb}",npname[i].Data(),sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
        //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
        
        scans->back()->SetLineWidth(3);
        scans->back()->SetLineColor(kOrange+1);
        scans->back()->Draw("L");
        DrawText("#scale[1.2]{Projection}",kRed,0.9,0.75);
        DrawText("#scale[1.2]{Manual scan}",kOrange+1,0.85,0.75);
        
        double pos_err = DetermineQuantile(scans->back(), 1, sigma->getVal(), sigma->getVal() + 3.*sigma->getError()) - sigma->getVal();
        double neg_err = -1*(DetermineQuantile(scans->back(), 1, sigma->getVal(), sigma->getVal() - 3.*sigma->getError()) - sigma->getVal());
        
        DrawText(Form("#scale[1]{Error from manual scan = + %4.4f - %4.4f pb}", pos_err, neg_err),kBlack,0.675-(2)*0.045,0.45);
        
        
        _can->Print(_ps); _can->Print(Form("plots/%s_%s.pdf",_set.Variable.Data(),npname[i].Data()));
    }
}

// ---------------------------------------------------------------------------------------
// Control plot functions


void HiggsPlot::PlotNPScans(RooFitResult *result, RooFormulaVar *nll, vector<TGraph*>* scans) {
    
            RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("Uncert_Dummy")));
        double range_low = sigma->getVal() - 3. * sigma->getError() < 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
        double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
        RooPlot* frame = sigma->frame(Name(Form("Uncert_Dummy")),Range(range_low,range_high),Title(Form("-log(L) scan vs Dummy")));
        frame->SetXTitle(Form("Dummy")); frame->SetYTitle("#it{Scan of '-2 log(L)'}");
        nll->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kRed ),LineStyle( kSolid ),Precision(1e-4));
        frame->SetMinimum(0.);
        frame->Draw();
        TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kBlack); error_p68.SetLineWidth(2); error_p68.SetLineStyle(7); error_p68.Draw("same");
        TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(7); error_p95.Draw("same");
        DrawText(Form("#scale[1.5]{Dummy = %4.4f #pm %4.4f pb}",sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
        //ATLAS_LABEL(0.21,0.88,1); DrawText("internal",kBlack,(0.93-0.038),0.325);
        
        scans->back()->SetLineWidth(3);
        scans->back()->SetLineColor(kOrange+1);
        scans->back()->Draw("L");
        DrawText("#scale[1.2]{Projection}",kRed,0.9,0.75);
        DrawText("#scale[1.2]{Manual scan}",kOrange+1,0.85,0.75);
        
        double pos_err = DetermineQuantile(scans->back(), 1, sigma->getVal(), sigma->getVal() + 3.*sigma->getError()) - sigma->getVal();
        double neg_err = -1*(DetermineQuantile(scans->back(), 1, sigma->getVal(), sigma->getVal() - 3.*sigma->getError()) - sigma->getVal());
        
        DrawText(Form("#scale[1]{Error from manual scan = + %4.4f - %4.4f pb}", pos_err, neg_err),kBlack,0.675-(2)*0.045,0.45);
        
        
        _can->Print(_ps); _can->Print(Form("plots/%s_Dummy.pdf",_set.Variable.Data()));
}


// ---------------------------------------------------------------------------------------
// Control plot functions

void HiggsPlot::PlotMeasurements() {
    
    VecI col; col.push_back( kAzure-3 ); col.push_back( kOrange+7 ); col.push_back(kRed); col.push_back(kMagenta);
    // Draw Input Data
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        Str xtit = _set.Label + "  " + _set.Unit;
        Str ytit = !in.hasBkg ? "#it{N^{sig}}" : "#it{N^{data}}" ;
        TH1F *axis = MakeAxis(in.Bins,xtit,ytit,0);
        TG *data = !in.hasBkg ? DrawTG(in.Nsig,in.NsigError,in.Bins, col.size() > i ? col[i] : kBlack )
        :  DrawTG(in.Ndata,SqrtVecD(in.Ndata),in.Bins, col.size() > i ? col[i] : kBlack) ;
        axis->SetMaximum( 1.5 *data->GetMaximum() );
        axis->Draw(); data->Draw("PE"); DrawText("#scale[2.0]{Input: #it{"+in.Name+"}}",kBlack,0.8,0.6);
        _can->Print(_ps); _can->Clear();
    }
    
}

void HiggsPlot::PlotDifferentialMeasurements() {
    VecI col; col.push_back( kAzure-3 ); col.push_back( kOrange+7 ); col.push_back(kRed); col.push_back(kMagenta);
    // Draw Input Data
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        Str xtit = _set.Label + "  " + _set.Unit;
        Str ytit = !in.hasBkg ? Form("d#it{N^{sig}} / d(#it{%s})", _set.Label.Data() ) : Form("d#it{N^{data}} / d(#it{%s})", _set.Label.Data() ) ;
        TH1F *axis = MakeAxis(in.Bins,xtit,ytit,0);
        TG *data = !in.hasBkg ? DrawTG(in.Nsig,in.NsigError,in.Bins, col.size() > i ? col[i] : kBlack, true )
        :  DrawTG(in.Ndata,SqrtVecD(in.Ndata),in.Bins, col.size() > i ? col[i] : kBlack, true) ;
        axis->SetMaximum( 1.5 *data->GetMaximum() );
        axis->Draw(); data->Draw("PE"); DrawText("#scale[2.0]{Input: #it{"+in.Name+"}}",kBlack,0.8,0.6);
        _can->Print(_ps); _can->Clear();
    }
}

void HiggsPlot::PlotCovariances() {
    TH2* h_sys  = DrawMatrix( (* _workspace->GetSyscovMatrix() ) );
    TH2* h_stat = DrawMatrix( (* _workspace->GetStatcovMatrix() ) );
    TH2* h_tot  = DrawMatrix( (* _workspace->GetTotcovMatrix() ) );
    h_sys->Draw("colztext");  _can->Print(_ps);
    h_stat->Draw("colztext"); _can->Print(_ps);
    h_tot->Draw("colztext");  _can->Print(_ps);
}

void HiggsPlot::PlotCorrelations() {
    TH2* h_sys  = DrawMatrix( (* _workspace->GetSyscorMatrix() ) );
    TH2* h_stat = DrawMatrix( (* _workspace->GetStatcorMatrix() ) );
    TH2* h_tot  = DrawMatrix( (* _workspace->GetTotcorMatrix() ) );
    h_sys->Draw("colztext");  _can->Print(_ps);
    h_stat->Draw("colztext"); _can->Print(_ps);
    h_tot->Draw("colztext");  _can->Print(_ps);
}

// ----------------------------------------------------------------------------------------

VecD HiggsPlot::Unfold(input in) {
    VecD xsec; int npoints = !in.hasBkg ? in.Nsig.size() : in.Ndata.size();
    for(int i = 0; i < npoints; i++) {
        if(in.hasBkg)
            xsec.push_back( (in.Ndata[i] - in.Nbkg[i]) / in.BR / in.Luminosity);
        else
            xsec.push_back( (in.Nsig[i]) / in.BR / in.Luminosity );
    }
    return xsec;
}

VecD HiggsPlot::UnfoldNormalized(input in) {
    VecD xsec = Unfold(in); int npoints = !in.hasBkg ? in.Nsig.size() : in.Ndata.size();
    double xsec_tot = Sum(xsec);
    for(int i = 0; i < npoints; i++)
        xsec[i] = xsec[i]/xsec_tot;
    return xsec;
}

VecD HiggsPlot::UnfoldStatError(input in, int index) {
    VecD xsecError; int pos = _workspace->GetCovPos(index), npoints = !in.hasBkg ? in.Nsig.size() : in.Ndata.size();
    for(int i = 0; i < npoints; i++)
        xsecError.push_back( sqrt( (*_workspace->GetStatcovMatrix())(i+pos,i+pos) ) );
    return xsecError;
}

VecD HiggsPlot::UnfoldSysError(input in, int index) {
    VecD xsecError; int pos = _workspace->GetCovPos(index), npoints = !in.hasBkg ? in.Nsig.size() : in.Ndata.size();
    for(int i = 0; i < npoints; i++)
        xsecError.push_back( sqrt( (*_workspace->GetSyscovMatrix())(i+pos,i+pos) ) );
    return xsecError;
}

VecD HiggsPlot::UnfoldError(input in, int index) {
    VecD xsecError; int pos = _workspace->GetCovPos(index), npoints = !in.hasBkg ? in.Nsig.size() : in.Ndata.size();
    for(int i = 0; i < npoints; i++)
        xsecError.push_back( sqrt( (*_workspace->GetTotcovMatrix())(i+pos,i+pos) ) );
    return xsecError;
}

VecD HiggsPlot::UnfoldNormalizedError(input in, int index) {
    VecD xsec = Unfold(in), xsecError = UnfoldError(in,index);
    int npoints = !in.hasBkg ? in.Nsig.size() : in.Ndata.size();
    double xsec_tot = Sum(xsec);
    for(int i = 0; i < npoints; i++)
        xsecError[i] = xsecError[i] / xsec_tot;
    return xsecError;
}

// ----------------------------------------------------------------------------------------
// Perform the proper error propagation from free parameters to

VecD HiggsPlot::ShapeErrors(VecD xsec, VecD errors, TMatrixD RooFitCor) {


    cout << "Input: " << endl;
    for(int i = 0; i < xsec.size(); i++) 
     cout << xsec[i] << "\t" <<  errors[i] << endl;

    int npoints = xsec.size(); TMatrixD A(npoints,npoints), Cov(npoints,npoints), Cor(npoints,npoints);
    double norm = Sum(xsec), norm2 = pow(Sum(xsec),2.);
    
    RooFitCor.Print("all");
    
    // Construct Jacobian
    for(int i = 0; i < npoints; i++)
        for(int j = 0; j < npoints; j++) {
            A(i,j) = i == j ? -xsec[i]/norm2 + 1/norm :
            -xsec[i]/norm2 ;
        }
    // Construct Covariance
    for(int i = 1; i < npoints; i++)
        for(int j = 1; j < npoints; j++) {
            Cov(i,j) = RooFitCor(RooFitCor.GetNrows()-xsec.size()+i-1, RooFitCor.GetNrows()-xsec.size()+j-1)*errors[i]*errors[j];
            Cor(i,j) = RooFitCor(RooFitCor.GetNrows()-xsec.size()+i-1,RooFitCor.GetNrows()-xsec.size()+j-1);
        }
    cout <<"Print constructed Correlation and Covariance in ShapeErrors"<<endl;
    Cor.Print("all");   
    Cov.Print("all");   
 
   
    for(int i = 0; i < errors.size(); i++)
        cout << xsec[i] << "\t" << errors[i] << endl;
    
    // Propagate the error
    //TMatrixD CovTr = A*Cov*A.T();
    cout <<"Print correlation and covariance with propagated error"<<endl;
    TMatrixD CovTr = A*Cov*A.T();
    VecD xsec_errors;
    for(int i = 0; i < npoints; i++)
        xsec_errors.push_back( sqrt(CovTr(i,i)) );
    CovTr.Print("all");
    TMatrixD CorTr = CovTr;
    for(int i = 0; i < npoints; i++)
        for(int j = 0; j < npoints; j++)
            CorTr(i,j) = CovTr(i,j)/sqrt(CovTr(i,i))/sqrt(CovTr(j,j));
    CorTr.Print("all");
    

    for(int i = 0; i < xsec_errors.size(); i++) 
      cout << Normalize(xsec)[i] << "\t" << xsec_errors[i]/ Normalize(xsec)[i] << endl;  
        
    return xsec_errors;
    
}

