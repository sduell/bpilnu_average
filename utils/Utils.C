//
// Utils Class
//

#include <iostream>
#include "Utils.h"
#include "TROOT.h"

///////////////////////////////////////////////////////////////////////////////////////
// ---- Many from Dag's utils file - thanks Dag! -----/
///////////////////////////////////////////////////////////////////////////////////////

static int hi=0;

void error(Str msg) {
  printf("ERROR:\n\n  %s\n\n",msg.Data()); 
  abort();
}


StrV Vectorize(Str str, Str sep) {
  StrV result; TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries()==0) { delete strings; return result; }
  TIter istr(strings);
  while (TObjString* os=(TObjString*)istr()) { 
    if (os->GetString()[0]=='#') break; 
    add(result,os->GetString()); 
  }
  delete strings; return result;
}

VecD VectorizeD(Str str, Str sep) {
  VecD result; StrV vecS = Vectorize(str,sep);
  for (uint i=0;i<vecS.size();++i) 
    result.push_back(atof(vecS[i]));
  return result;
}

VecI VectorizeI(Str str, Str sep) {
  VecI result; StrV vecS = Vectorize(str,sep);
  for (uint i=0;i<vecS.size();++i) {
    result.push_back(atoi(vecS[i]));
    }
  return result;
}

TEnv *OpenSettingsFile(Str fileName) {
  if (fileName=="") error("No config file name specified. Cannot open file!");
  TEnv *settings = new TEnv();
  int status=settings->ReadFile(fileName.Data(),EEnvLevel(0));
  if (status!=0) error(Form("Cannot read file %s",fileName.Data()));
  return settings;
}

TFile *OpenFile(Str fn) {
  TFile *f = TFile::Open(fn); if (f==NULL) error("Cannot open "+fn); return f; 
}

TTree *GetTree(TFile *f, Str tn) {
  TTree *t = (TTree*)f->Get(tn); if (t==NULL) error("Cannot access tree "+tn+" in "+f->GetName());
  return t;
}

TTree *GetTree(Str fn, Str tn) {
//   std::cout << "Will open " << fn << std::endl;
  return GetTree(OpenFile(fn),tn);
}

StrV ReadFile(TString fileName) {
  StrV lines;
  ifstream file(fileName.Data());
  if (!file.good()) error("Cannot open file "+fileName);
  string line, lastline="weeee";
  while (getline(file,line)) {
    if (line==lastline) continue;
    //if (line[0]==' ') continue;
    StrV subLines=Vectorize(line,",");
    for (uint i=0;i<subLines.size();++i) 
      lines.push_back(subLines[i]);
  }
  return lines;
}

TGraphErrors* Draw1D(VecD content, VecD errors, Str var, int Nbins, VecD bins, VecD binlims, 
	    int col, int ms, Str xtit, Str ytit, bool orn) {    
  Str hname(Form("hist%d",++hi)); // unique name  
  int n = bins.size(); double x[n], y[n], ex[n], ey[n];   
  for(int i = 0; i < bins.size(); i++) { x[i] = bins[i]; y[i] = content[i]; ex[i] = 0.0; ey[i] = errors[i]; }
  TGraphErrors *h; orn == true ? h = new TGraphErrors(n,x,y,ex,ey) : h = new TGraphErrors(n,y,x,ey,ex);
  h->GetXaxis()->SetTitle(xtit); h->GetYaxis()->SetTitle(ytit);  
  h->GetXaxis()->SetLimits(binlims[0]-1.0,binlims[binlims.size()-1]+1.0);  
  if (h==NULL) error("Cannot access TGraphErrors: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  return h;
}

TGraphAsymmErrors* Draw1D(VecD content, double scale, int Nbins, VecD bins, VecD binlims, 
	    int col, int ms, Str xtit, Str ytit, bool orn) {    
  Str hname(Form("hist%d",++hi)); // unique name  
  int n = bins.size(); double x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];   
  for(int i = 0; i < bins.size(); i++) { x[i] = bins[i]; y[i] = content[i]*scale; exl[i] = 0.0; eyl[i] = 0.0; exh[i] = 0.0; eyh[i] = 0.0; } 
  TGraphAsymmErrors *h; orn == true ? h = new TGraphAsymmErrors(n,x,y,exl,eyl,exh,eyh) : h = new TGraphAsymmErrors(n,y,x,exl,eyl,exh,eyh);
  h->GetXaxis()->SetTitle(xtit); h->GetYaxis()->SetTitle(ytit);  
  h->GetXaxis()->SetLimits(binlims[0]-1.0,binlims[binlims.size()-1]+1.0);  
  if (h==NULL) error("Cannot access TGraphErrors: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  return h;
}


TGraphAsymmErrors* Draw1D(VecD content, VecD errors_low, VecD errors_up, Str var, int Nbins, VecD bins, VecD binlims, 
	    int col, int ms, Str xtit, Str ytit, bool orn) {    
  Str hname(Form("hist%d",++hi)); // unique name  
  int n = bins.size(); double x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];   
  for(int i = 0; i < bins.size(); i++) { x[i] = bins[i]; y[i] = content[i]; exl[i] = 0.0; eyl[i] = errors_low[i]; exh[i] = 0.0; eyh[i] = errors_up[i]; }
  TGraphAsymmErrors *h; orn == true ? h = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh) : h = new TGraphAsymmErrors(n,y,x,eyl,eyh,exl,exh);
  h->GetXaxis()->SetTitle(xtit); h->GetYaxis()->SetTitle(ytit);  
  h->GetXaxis()->SetLimits(binlims[0]-1.0,binlims[binlims.size()-1]+1.0);  
  if (h==NULL) error("Cannot access TGraphErrors: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  return h;
}

TGraphAsymmErrors* Draw1D(VecD content, Str var, int Nbins, VecD bins, VecD binlims, 
	    int col, int ms, Str xtit, Str ytit, bool orn) {    
  Str hname(Form("hist%d",++hi)); // unique name  
  int n = bins.size(); double x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];   
  for(int i = 0; i < bins.size(); i++) { x[i] = bins[i]; y[i] = content[i]; exl[i] = 0.0; eyl[i] = 0.0; exh[i] = 0.0; eyh[i] = 0.0; }
  TGraphAsymmErrors *h; orn == true ? h = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh) : h = new TGraphAsymmErrors(n,y,x,eyl,eyh,exl,exh);
  h->GetXaxis()->SetTitle(xtit); h->GetYaxis()->SetTitle(ytit);  
  h->GetXaxis()->SetLimits(binlims[0]-1.0,binlims[binlims.size()-1]+1.0);  
  if (h==NULL) error("Cannot access TGraphErrors: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  return h;
}


TH1* Draw1D( TTree *t, Str var, Str cut, 
	     int Nbins, double min, double max, 
	     int col, int sty, int ms, Str xtit, Str ytit);

// TH1* Draw2D( TTree *t, Str var, Str cut="1", 
// 	     int Nbins=100, double xmin=-1, double xmax=1, 
// 	     double ymin=-1, double ymax=1, int col=kBlack, int ms=1, Str xtit="", Str ytit="");

TH1* Draw1D(TFile *f, Str var, Str cut, int Nbins, double min, double max, 
	    int col, int ms, Str xtit, Str ytit) {
  TTree *t = (TTree*)f->Get("BDT_Tree"); if (t==NULL) error("can't find tree BDT_Tree");
  return Draw1D(t,var,cut,Nbins,min,max,col,1,ms,xtit,ytit);
}

TH1* Draw1D(Str fn, Str hist, Str var, int col, int sty, int ms, Str xtit, Str ytit) {
   TFile *f = new TFile(fn);
   TH1D *t = (TH1D*)f->Get(hist); if (t==NULL) error("can't find hist "+hist);
   return Draw1D(t,var,col,1,ms,xtit,ytit);
}

TH1* Draw1D(TTree *t, Str var, Str cut, int Nbins, double min, double max, 
	    int col, int sty, int ms, Str xtit, Str ytit) {
  Str hname(Form("hist%d",++hi)); // unique name 
  t->Draw(var+Form(">>%s(%d,%.4f,%.4f)",hname.Data(),Nbins,min,max),cut,"goffe");
  TH1* h = (TH1*)gROOT->FindObject(hname);  
  if (h==NULL) error("Cannot access histo: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(sty); h->SetStats(0); h->SetXTitle(var); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  if (xtit!="") h->SetXTitle(xtit); if (ytit!="") h->SetYTitle(ytit);
  return h;
}

TH1* Draw1D(TTree *t, Str var, Str cut, int Nbins, double *xbins,
	    int col, int sty, int ms, Str xtit, Str ytit) {
  Str hname(Form("hist%d",++hi)); // unique name 
  TH1* h = new TH1D(hname,hname,Nbins,xbins);
  t->Draw(var+Form(">>%s",hname.Data()),cut,"goffe");
  if (h==NULL) error("Cannot access histo: "+hname);
  for(int i = 0; i < h->GetXaxis()->GetNbins(); i++)
   h->SetBinContent(i+1,h->GetBinContent(i+1));
  h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(sty); h->SetStats(0); h->SetXTitle(var); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  if (xtit!="") h->SetXTitle(xtit); if (ytit!="") h->SetYTitle(ytit);
  return h;
}



TH1* Draw1D(TH1D *h, Str var, int col, int sty, int ms, Str xtit, Str ytit) {
  if (h==NULL) error("Cannot access histo");
  h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(sty); h->SetStats(0); h->SetXTitle(var); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  if (xtit!="") h->SetXTitle(xtit); if (ytit!="") h->SetYTitle(ytit);
  return h;
}

TH2* Draw2D(TTree *t, Str var, Str cut, int Nbins, double xmin, double xmax, 
            double ymin, double ymax, int col, int ms, Str xtit, Str ytit) {
  Str hname(Form("hist%d",++hi)); // unique name  
  t->Draw(var+Form(">>%s(%d,%.4f,%.4f,%d,%.4f,%.4f)",hname.Data(),Nbins,xmin,xmax,Nbins,ymin,ymax),cut,"goffe");
  TH2* h = (TH2*)gROOT->FindObject(hname);  
  if (h==NULL) error("Cannot access histo: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); h->SetStats(0); h->SetXTitle(var); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  if (xtit!="") h->SetXTitle(xtit); if (ytit!="") h->SetYTitle(ytit);
  return h;
}

TH2* Draw2D(TTree *t, Str var, Str cut, int col, int ms, Str xtit, Str ytit) {
  Str hname(Form("hist%d",++hi)); // unique name  
  t->Draw(var+Form(">>%s",hname.Data()),cut,"goffe");
  TH2* h = (TH2*)gROOT->FindObject(hname);  
  if (h==NULL) error("Cannot access histo: "+hname);
  h->SetLineColor(col); h->SetLineWidth(2); h->SetStats(0); h->SetXTitle(var); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  if (xtit!="") h->SetXTitle(xtit); if (ytit!="") h->SetYTitle(ytit);
  return h;
}


void DrawText(TString txt, int col, double y, double x, int align)
{ 
  static TLatex *_tex = new TLatex(); _tex->SetNDC(); _tex->SetTextSize(0.03);
  _tex->SetTextAlign(align); _tex->SetTextColor(col); _tex->DrawLatex(x,y,txt); 
}

void DrawText(TString txt, int col, int style, double y, double x, int align)
{ 
 static TLatex *_tex = new TLatex(); _tex->SetNDC(); _tex->SetTextSize(0.03);
  _tex->SetTextAlign(align); _tex->SetTextColor(col); _tex->DrawLatex(x,y,txt); 
  return;
  TMarker *marker = new TMarker(x-(0.4*1),y,8);
  marker->SetMarkerColor(col);  marker->SetNDC();
  marker->SetMarkerStyle(style);
  marker->SetMarkerSize(style);
  //marker->Draw("same");
}

void DrawSubPad() {
  double TEXTSIZE=0.12, FIGURE2_RATIO = 0.36, SUBFIGURE_MARGIN = 0.18; // white space between figures
  double _maxDev=0.25;
  gPad->SetBottomMargin(FIGURE2_RATIO); gPad->SetTopMargin(0.0);
  // create new pad, fullsize to have equal font-sizes in both plots
  TPad *p = new TPad( "p_test", "", 0, 0.0, 0.99, 1.0 - SUBFIGURE_MARGIN, 0, 0, 0);
  p->SetTopMargin(1.0 - FIGURE2_RATIO); p->SetFillStyle(0);
  p->SetMargin(0.16,0.04,0.18,1.0 - FIGURE2_RATIO);
  p->Draw(); p->cd(); p->SetGridy(kTRUE);
}

TH1* Draw1D(VecD Values, VecD Errors, VecD Bins,int col, int sty, int ms, Str xtit, Str ytit) {
  TH1D *hist = new TH1D(Form("hist_%d",++hi),"",Bins.size()-1,&Bins[0]);
  hist->SetXTitle(xtit); hist->SetYTitle(ytit);
  hist->SetLineColor(col); hist->SetLineWidth(2); hist->SetLineStyle(sty); hist->SetStats(0); 
  hist->SetMarkerStyle(ms); hist->SetMarkerSize(0.8); hist->SetMarkerColor(col); hist->SetMinimum(0);
  for(int i = 0; i < Values.size(); i++) {
   hist->SetBinContent(i+1,Values[i]);
   hist->SetBinError(i+1,Errors[i]);
  }
  return hist;
}

TG* DrawTG(VecD Values, VecD Errors, VecD Bins, int col, bool diff, double offset, int sty, int ms, Str xtit, Str ytit) {
  int n = Bins.size()-1; double x[n], y[n], ex[n], ey[n];  
  VecD BinCenter; for(int i = 0; i < n; i++) BinCenter.push_back( Bins[i] + (Bins[i+1]-Bins[i])*offset );
  VecD BinError; for(int i = 0; i < n; i++) BinError.push_back( (Bins[i+1]-Bins[i])/2 );
  double max(-99999);
  VecD UsedValues, UsedErrors; 
  if(!diff) { for(int i = 0; i < Values.size(); i++) { UsedValues.push_back( Values[i] ); UsedErrors.push_back( Errors[i] ); }  }
   else {  for(int i = 0; i < Values.size(); i++) { UsedValues.push_back( Values[i] / (2 * BinError[i]) ); UsedErrors.push_back( Errors[i] / ( 2 * BinError[i])  ); } }
  for(int i = 0; i < n; i++) { x[i] = BinCenter[i]; y[i] = UsedValues[i]; ex[i] = BinError[i]; ey[i] = UsedErrors[i]; UsedValues[i] > max ? max = UsedValues[i] : max = max; }
  TGraphErrors *h; h = new TGraphErrors(n,x,y,ex,ey);
  h->SetMaximum(max);
  h->GetXaxis()->SetTitle(xtit); h->GetYaxis()->SetTitle(ytit);  
  h->SetLineColor(col); h->SetLineWidth(2); 
  h->SetMarkerStyle(ms); h->SetMarkerSize(0.8); h->SetMarkerColor(col); h->SetMinimum(0);
  return h;
}

TH1D* residualHist(const RooHist* rhist, const RooCurve* curve) {

  double r = 0.2; double sr = 1. / r;
  // Grab info from the histogram.
  int n = rhist->GetN(); double *x = rhist->GetX(); double *y = rhist->GetY(); 
  // Create residual histogram.
  double xMin = x[0]; double xMax = x[n-1];
  TH1D* residuals_temp = new TH1D("r","",n,xMin,xMax);
  double datum = 0.; double pdf = 0.;

  // Fill the histogram.
   if ( curve ) 
    for ( int bin = 0; bin < n; bin++ ) {
      datum = y[bin]; pdf = curve->Eval(x[bin]);
      //residuals_temp->SetBinContent(bin+1,pdf-datum);
      //residuals_temp->SetBinError(bin+1,sqrt(datum));
      //residuals_temp->SetBinContent(bin+1,(pdf-datum)/sqrt(datum));
      //residuals_temp->SetBinError(bin+1,0.001);      
      residuals_temp->SetBinContent(bin+1,residual(datum,pdf));
      residuals_temp->SetBinError(bin+1,0.001);
     }

  //residuals_temp->SetMinimum(-2.*residuals_temp->GetMaximum()); residuals_temp->SetMaximum(2.*residuals_temp->GetMaximum());
  residuals_temp->SetMinimum(-4.); residuals_temp->SetMaximum(4.);
  residuals_temp->SetMarkerStyle(8); residuals_temp->SetMarkerSize(0.8);

  residuals_temp->GetYaxis()->SetTitle("pulls   ");
  residuals_temp->GetYaxis()->SetNdivisions(4);
  residuals_temp->GetXaxis()->SetLabelSize(0.06); residuals_temp->GetYaxis()->SetLabelSize(0.06); 
     
  return residuals_temp;
  
}

double residual( double datum, double pdf) {
    double chi2 = 0.;
    if ( pdf > 0 ) chi2 += 2. * ( pdf - datum );
    if ( datum > 0 && pdf > 0 ) chi2 += 2. * datum * log( datum / pdf );
    return ( ( datum >= pdf ) ? sqrt( chi2 ) : -sqrt( chi2 ) );
}

VecD MakeUniformVecD(int N, double min, double max) {
  VecD vec; double dx=(max-min)/N;
  for (int i=0;i<=N;++i) vec.push_back(min+i*dx);
  return vec;
}

TH1F *MakeAxis(VecD bins, Str xtit, Str ytit, double value, bool sub) {
  static int s_i=0; TH1F *axis = new TH1F(Form("axis_%d",++s_i),"",bins.size()-1,&bins[0]);
  axis->SetXTitle(xtit); axis->SetYTitle(ytit);
//   if (xtit.Contains("#it{N}")&&bins.size()==5) {
//     for (int i=0;i<3;++i) axis->GetXaxis()->SetBinLabel(i+1,Form("%d",i));
//     axis->GetXaxis()->SetBinLabel(4,"#geq3");
//     axis->GetXaxis()->SetLabelSize(0.07);
//     axis->GetXaxis()->SetLabelOffset(0.02);
//   } else if (bins.size()==4) {
//     for (int i=1;i<=3;++i) axis->GetXaxis()->SetBinLabel(i,Form("#sigma_{%d} / #sigma_{#geq %d}",i-1,i-1));
//     for (int i=1;i<=3;++i) axis->SetBinContent(i,value);
//     axis->GetXaxis()->SetLabelSize(0.07); axis->GetXaxis()->SetLabelOffset(0.02);
//   }
  for (int i=1;i<bins.size();++i) axis->SetBinContent(i,value);
  axis->GetYaxis()->SetTitleOffset(1.0); axis->GetXaxis()->SetTitleOffset(1.0);
  if (sub) { axis->GetXaxis()->SetTitleOffset(20); axis->GetXaxis()->SetLabelOffset(20); }
  return axis;
}

TH2* DrawMatrix(TMatrixDSym mat) {
  TH2D *h_mat = new TH2D(Form("hist2_%d",++hi),"",mat.GetNrows(),1,mat.GetNrows(),mat.GetNcols(),1,mat.GetNcols()); 
  for(int i = 0; i < mat.GetNrows(); i++)
   for(int j = 0; j < mat.GetNcols(); j++)
    h_mat->SetBinContent(i+1,j+1, mat(i,j));
  return h_mat;
}

TH1F *MakeAxis(int N, double xmin, double xmax, Str xtit, Str ytit, double value) 
{ return MakeAxis(MakeUniformVecD(N,xmin,xmax),xtit,ytit,value); }

double DetermineQuantile(TH1* h, double observation) {
  TH1* h_norm = (TH1*) h->Clone(); h_norm->Scale(1.0/h_norm->Integral("width")); 
  int binobs = h_norm->FindBin(observation); double d_sum(0);  
  for(int i = 0; i < binobs; i++)
    d_sum += h_norm->GetBinContent(i+1)*h_norm->GetXaxis()->GetBinWidth(i+1);    
  return d_sum / h_norm->Integral("width");
}

double GetQuantile(TH1* h, double median, double quantile) {
  TH1* h_norm = (TH1*) h->Clone(); h_norm->Scale(1.0/h_norm->Integral("width")); 
  int binobs = h_norm->FindBin(median); double d_sum = 0, error = 0; int j(binobs);
  d_sum += h_norm->GetBinContent(binobs)*h_norm->GetXaxis()->GetBinWidth(binobs);
  for(int i = binobs; i <= h_norm->GetXaxis()->GetNbins(); i++) {
    d_sum += h_norm->GetBinContent(i)*h_norm->GetXaxis()->GetBinWidth(i)
             + h_norm->GetBinContent(j)*h_norm->GetXaxis()->GetBinWidth(j); 
    j--;         
    //cout << median << "\t" << quantile << "\t" << d_sum << endl;         
    if(d_sum > quantile && error == 0) error = fabs(median-h_norm->GetBinCenter(i));
  }
  // cout << " error = " << error << endl;
  return error;
}

double GetOnesidedQuantile(TH1* h, double quantile) {
  TH1* h_norm = (TH1*) h->Clone(); h_norm->Scale(1.0/h_norm->Integral("width")); 
  double d_sum = 0, error = 0;
  for(int i = 0; i <= h_norm->GetXaxis()->GetNbins(); i++) {
    d_sum += h_norm->GetBinContent(i)*h_norm->GetXaxis()->GetBinWidth(i); 
    if(d_sum > quantile && error == 0) error = h_norm->GetBinCenter(i);
  }
  return error;
}


double DetermineMedian(TH1 *h) { 
   int nbins = h->GetXaxis()->GetNbins(); 
   double *x = new double[nbins], *y = new double[nbins]; 
   for (Int_t i=0;i<nbins;i++) {
      x[i] = h->GetXaxis()->GetBinCenter(i+1); 
      y[i] = h->GetBinContent(i+1); 
   } 
   Double_t median = TMath::Median(nbins,x,y); 
   delete [] x; 
   delete [] y; 
   return median; 
} 

// Normal chisq(x;ndf) Function
double myProbFunction(double *_x, double *_par) {
 double x = _x[0], ndf = _par[0], norm = _par[1]; 
 return norm * ( ROOT::Math::chisquared_pdf(x,ndf) );
}

// Function for finite parts of 0.5 * delta(x) + 0.5 * chisq(x;ndf)
double myProbFunction2(double *_x, double *_par) {
 double x = _x[0], ndf = _par[0], norm = _par[1], param = _par[2]; 
 return x == 0? param : norm * ( ROOT::Math::chisquared_pdf(x,ndf) );
}

// Function for finite parts of Eq. 8 in my write-up
double myProbFunction3(double *_x, double *_par) {
 double x = _x[0], ndf = _par[0], norm = _par[1], param = _par[2], mu = _par[3], sigma = _par[4];
 return x == 0 ? param :
           // mu < 0.3 ? 0.5 * norm * ( 1. / sqrt(2*TMath::Pi()) * 1. / sqrt(x) * exp( -0.5*pow( sqrt(x) - mu/sigma, 2. ) ) ) + 0.5* ( 1. / sqrt(2*TMath::Pi()) * 1. / (2 * mu / sigma) * exp( -0.5*pow( x - pow(mu/sigma,2.) ,2.) / ( pow( 2 * mu / sigma ,2.) ) ) ) :
           x <= pow(mu/sigma,2.) ?  
           norm * ( 1. / sqrt(2*TMath::Pi()) * 1. / sqrt(x) * exp( -0.5*pow( sqrt(x) - mu/sigma, 2. ) ) ) :
          ( 1. / sqrt(2*TMath::Pi()) * 1. / (2 * mu / sigma) * exp( -0.5*pow( x - pow(mu/sigma,2.) ,2.) / ( pow( 2 * mu / sigma ,2.) ) ) );
}

// Function for finite parts of Eq. 5 in my write-up
double myProbFunction4(double *_x, double *_par) {
 double x = _x[0], ndf = _par[0], norm = _par[1], param = _par[2], mu = _par[3], sigma = _par[4];
 return x == 0? param : x < pow(mu/sigma,2.) ?
        norm * ( ROOT::Math::chisquared_pdf(x,ndf) ) : 
        ( 1. / sqrt(2*TMath::Pi()) * 1. / (2 * mu / sigma) * exp( -0.5*pow( x + pow(mu/sigma,2.) ,2.) / ( pow( 2 * mu / sigma ,2.) ) ) );
}

VecD myVecDivide(VecD one, VecD two, double a, double b, double c) {
 VecD res; if(one.size() != two.size()) Warning("inconsistent vector sizes passed to myVecDivide call","");
 for(int i = 0; i < one.size(); i++) { 
   if( (a + b*two[i] ) > 0 ) res.push_back( (one[i] / (a + b*two[i] ) > c && c != -1) ? c : one[i] / (a + b*two[i] ) );
    else res.push_back(1.);
 }
 return res;
}

VecD myVecDivideError(VecD one, VecD two, VecD oneError, VecD twoError, double a, double b) {
 VecD res; if(one.size() != two.size()) Warning("inconsistent vector sizes passed to myVecDivide call","");
 for(int i = 0; i < one.size(); i++) {
  if( (a + b*two[i] ) > 0 ) res.push_back( sqrt( pow( 1./(a+b*two[i]) ,2.) * pow( oneError[i], 2.) 
                             + pow( (-b*one[i]) / ( pow( a + b*two[i] ,2.) ) , 2. ) * pow( twoError[i], 2.) ) );
  else res.push_back(0.);
 }                     
 return res;
}

TGraph* _graph; TF1* _func; double _quant;

double brent_graph_func(double x) { 
 return pow(_graph->Eval(x) - _quant,2.);
}

double DetermineQuantile(TGraph* graph, double quantile, double xmin, double xmax) {
   // Set global variables
   _graph = graph; _quant = quantile;
   // Initialize function
   ROOT::Math::Functor1D func(&brent_graph_func);
   ROOT::Math::BrentMinimizer1D bm;
   bm.SetFunction(func,xmin,xmax);
   bm.Minimize(1000,0,0);
   
   cout << bm.XMinimum() << endl;
   
   return bm.XMinimum();
}

double brent_func(double x) { 
 return pow(_func->Eval(x) - _quant,2.);
}

double DetermineQuantile(TF1* tf_func, double quantile, double xmin, double xmax) {
   // Set global variables
   _func = tf_func; _quant = quantile;
   // Initialize function
   ROOT::Math::Functor1D func(&brent_func);
   ROOT::Math::BrentMinimizer1D bm;
   bm.SetFunction(func,xmin,xmax);
   bm.Minimize(1000,0,0);
   
   cout << bm.XMinimum() << endl;
   
   return bm.XMinimum();
}

VecD SqrtVecD(VecD Values) {
 VecD Sqrt; for(int i = 0; i < Values.size(); i++) Sqrt.push_back( sqrt(Values[i]) );
 return Sqrt;
}

double Sum(VecD vec) {
  double sum(0);
  for(int i = 0; i < vec.size(); i++)
   sum += vec[i];
  return sum; 
}

VecD Normalize(VecD vec) {
  double sum = Sum(vec);
  for(int i = 0; i < vec.size(); i++)
   vec[i] = vec[i]/sum;   
  return vec; 
}

VecD NormalizeError(VecD vec, VecD errors) {
  double sum = Sum(vec);
  for(int i = 0; i < vec.size(); i++)
   errors[i] = errors[i]/sum;
  return errors; 
}

double KroneckerDelta(int i, int j) { return i == j ? 1.0 : 0.; };

///////////////////////////////////////////////////////////////////////////////////////
