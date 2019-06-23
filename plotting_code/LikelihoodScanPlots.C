/*
 *  LikelihoodScanPlots: Michaela Queitsch-Maitland
 */

#include "../utils/Utils.C"
#include "../utils/AtlasStyle.C"
#include "../utils/AtlasLabels.C"
#include "../utils/AtlasUtils.C"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooMultiVarGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
using namespace RooFit;

void fatal(Str msg) { printf("\nFATAL\n  %s\n\n",msg.Data()); abort(); }

int GetNBins(Str var) {
  if (var=="pTH") return 8;
  else if (var=="yAbsH") return 6;
  else if (var=="Njets") return 4;
  else if (var=="pTj1") return 5;
  else if (var=="Inclusive") return 1;
  else {
    fatal(Form("Could not understand variable %s!",var.Data()));
    return -1;
  }
}

Str VarLabel(Str var) {
  if (var=="pTH") return "p_{T}^{H}";
  else if (var=="yAbsH") return "|y^{H}|";
  else if (var=="Njets") return "N_{jets}";
  else if (var=="pTj1") return "p_{T}^{j1}";
  else if (var=="Inclusive") return "Inclusive";
  else {
    fatal(Form("Could not understand variable %s!",var.Data()));
    return -1;
  }
}

int main(int argc, char *argv[]) {

  if (argc<2) fatal(Form("Only provided %d args!",argc));

  SetAtlasStyle(); gStyle->SetPalette(1); gStyle->SetHistMinimumZero();

  for (int i=1; i<argc; ++i) {
    Str var = argv[i];
    int nbins = GetNBins(var);

    Str ifn_fixed = Form("results/%s_fixNPs_manualScan.root",var.Data());
    Str ifn_float = Form("results/%s_floatNPs_manualScan.root",var.Data());

    if ( (gSystem->AccessPathName(ifn_fixed)) )
	 fatal(Form("Could not read file: %s",ifn_fixed.Data()));
    if ( (gSystem->AccessPathName(ifn_float)) )
	 fatal(Form("Could not read file: %s",ifn_float.Data()));
    
    TFile f(ifn_fixed,"read");
    RooFitResult* result_fixNPs = (RooFitResult*)f.Get("combined_result");
    vector<TGraph*> scans_fixnps;
    for (int i=0; i<nbins; ++i)
      scans_fixnps.push_back( (TGraph*)f.Get( Form("scan_sigma%d",i+1) ) );
    
    TFile f_nps(ifn_float,"read");
    RooFitResult* result = (RooFitResult*)f_nps.Get("combined_result");
    vector<TGraph*> scans;
    for (int i=0; i<nbins; ++i)
      scans.push_back( (TGraph*)f_nps.Get( Form("scan_sigma%d",i+1) ) );
    RooFormulaVar* nll_yy = (RooFormulaVar*)f_nps.Get("yy_nll");
    RooFormulaVar* nll_zz = (RooFormulaVar*)f_nps.Get("zz_nll");
    
    TCanvas *can = new TCanvas();
    can->Draw();
    can->Print(Form("%s_likelihood_scans.pdf[",var.Data()));
    
    for (int i=0; i<nbins; ++i) {
      //    RooRealVar *sigma = (RooRealVar*)(result_fixNPs->floatParsFinal().find(Form("sigma_%d",i+1)));   
      RooRealVar *sigma = (RooRealVar*)(result->floatParsFinal().find(Form("sigma_%d",i+1)));   
      double range_low = sigma->getVal() - 3. * sigma->getError() > 0 ? sigma->getVal() - 3. * sigma->getError() : 0. , range_high = sigma->getVal() + 3. * sigma->getError();
      double range_low_p = sigma->getVal() - sigma->getError(), range_high_p = sigma->getVal() + sigma->getError();
      RooPlot* frame = sigma->frame(Name(Form("nll_frame_%d",i+1)),Range(range_low,range_high),Title(Form("-log(L) scan vs #sigma_{%d}",i+1)));
      frame->SetXTitle(Form("#sigma_{%d}  [pb]",i+1)); frame->SetYTitle("#it{-2 log(L)}");
      nll_zz->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kOrange+7 ),LineStyle( kSolid ),Precision(1e-4));
      nll_yy->plotOn(frame,PrintEvalErrors(0),ShiftToZero(),LineColor( kAzure-3 ),LineStyle( kSolid ),Precision(1e-4));
      frame->SetMinimum(0.);
      frame->SetMaximum(12.);
      frame->Draw();
      TArrow error_p68(range_low,1,range_high,1,0.15,""); error_p68.SetLineColor(kGray+1); error_p68.SetLineWidth(2); error_p68.SetLineStyle(1); error_p68.Draw("same");
      TArrow error_p95(range_low,4,range_high,4,0.15,""); error_p95.SetLineColor(kGray); error_p95.SetLineWidth(2); error_p95.SetLineStyle(1); error_p95.Draw("same");
      DrawText(Form("#scale[1.5]{%s}",VarLabel(var).Data()),kBlack,0.68,0.45);
      DrawText(Form("#scale[1.2]{#sigma_{%d} = %4.2f #pm %4.2f pb}",i+1,sigma->getVal(),sigma->getError()),kBlack,0.675-(1)*0.045,0.45);
      ATLAS_LABEL(0.21,0.88,1); DrawText("#scale[1.5]{internal}",kBlack,(0.93-0.038),0.325); 
      DrawText("#scale[1.2]{H #rightarrow #gamma#gamma}",kAzure-3,0.9,0.7);
      DrawText("#scale[1.2]{H #rightarrow ZZ}",kOrange+7,0.85,0.7);
      DrawText("#scale[1.2]{Combined (stat+syst)}",kGray+2,0.8,0.7);
      DrawText("#scale[1.2]{Combined (stat)}",kGray+1,0.75,0.7);
      
      scans[i]->SetLineWidth(3);
      scans[i]->SetLineColor(kGray+2);
      scans[i]->Draw("L");
      
      scans_fixnps[i]->SetLineWidth(3);
      scans_fixnps[i]->SetLineColor(kGray+1);
      scans_fixnps[i]->SetLineStyle(2);
      scans_fixnps[i]->Draw("L");
      
      double pos_err = DetermineQuantile(scans[i], 1, sigma->getVal(), sigma->getVal() + 3.*sigma->getError()) - sigma->getVal();
      double neg_err = -1*(DetermineQuantile(scans[i], 1, sigma->getVal(), sigma->getVal() - 3.*sigma->getError()) - sigma->getVal());
      DrawText(Form("manual scan #scale[1.2]{#Delta#sigma_{%d} = {}^{+%4.2f}_{-%4.2f} pb}", i+1, pos_err, neg_err),kBlack,0.675-(2)*0.045,0.45);
      
      can->SaveAs(Form("plots/%s_sigma_%d_likelihood_scans.pdf",var.Data(),i+1));
      can->Print(Form("%s_likelihood_scans.pdf",var.Data()));
    }
    
    can->Print(Form("%s_likelihood_scans.pdf]",var.Data()));
  }

}
