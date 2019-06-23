#ifndef _UTILS
#define _UTILS

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <time.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TEnv.h"
#include "TTree.h"
#include "TChain.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TH2.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGraph2DErrors.h>
#include <TLatex.h>
#include <TGraphQQ.h>
#include "TArrow.h"
#include "TPad.h"
#include <TVectorD.h>
#include "TF1.h"
#include "TMarker.h"

#include "RooHist.h"
#include "RooCurve.h"

#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "AtlasUtils.h"
#include <vector>

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"

using namespace std;

typedef TString Str;
typedef unsigned int uint;
typedef vector<TString> StrV;
typedef vector< vector <TString> > StrVV;
typedef vector< vector < vector <TString> > > StrVVV;
typedef vector<float> VecF;
typedef vector<int> VecI;
typedef vector<double> VecD;
typedef vector <vector <double> > VecVecD;
typedef vector <vector<vector <double> > > VecVecVecD;
typedef vector<TH1*> VecTH1;
typedef vector<vector <TH1*> > VecVecTH1;
typedef vector<TGraph*> VecTG;
typedef vector<TGraphAsymmErrors*> VecTGA;
typedef TGraphAsymmErrors TGA;
typedef TGraphErrors TG;
typedef vector<TMatrixD*> VecTMD;
typedef vector<TF1*> VecTF1;
typedef vector<TTree*> VecTTree;
typedef vector<TChain*> VecTChain;
typedef vector <vector<TTree*> > VecVecTTree;
typedef vector <vector<TChain*> > VecVecTChain;
typedef TLorentzVector TLV;
typedef std::vector<TLorentzVector> TLVs;
typedef vector <TVectorD> VecTVectorD;

///////////////////////////////////////////////////////////////////////////////////////

void DrawText(TString txt, int col=kBlack, double y=0.88, double x=0.145, int align=12);
void DrawText(TString txt, int col, int style, double y, double x, int align);

void error(Str msg);

template <class T> void add(vector<T> &vec, T a) { vec.push_back(a); };
template <class T> void add(vector<T> &v, T a, T b) { add(v,a); add(v,b); };
template <class T> void add(vector<T> &v, T a, T b, T c) { add(v,a,b); add(v,c); };
template <class T> void add(vector<T> &v, T a, T b, T c, T d) { add(v,a,b); add(v,c,d); };
template <class T> void add(vector<T> &vec, T a[]) { 
  uint n=sizeof(a)/sizeof(a[0]);
  for (uint i=0;i<n;++i) vec.push_back(a[i]);
}
template <class T> vector<T> vec(T a) { vector<T> v; add(v,a); return v; };
template <class T> vector<T> vec(T a, T b) { vector<T> v; add(v,a,b); return v; };
template <class T> vector<T> vec(T a, T b, T c) { vector<T> v; add(v,a,b,c); return v; };

StrV Vectorize(Str str, Str sep=" ");
VecD VectorizeD(Str str, Str sep);
VecI VectorizeI(Str str, Str sep);
TEnv *OpenSettingsFile(Str fileName);
TFile *OpenFile(Str fn);
TTree *GetTree(TFile *f, Str tn="BDT_Tree");
TTree *GetTree(Str fn, Str tn="BDT_Tree");
StrV ReadFile(TString fileName);

TGraphErrors* Draw1D(VecD content, VecD errors, Str var, int Nbins, VecD bins, VecD binlims, int col, int ms, Str xtit, Str ytit, bool orn);
TGraphAsymmErrors* Draw1D(VecD content, VecD errors_low, VecD errors_up, Str var, int Nbins, VecD bins, VecD binlims, int col, int ms, Str xtit, Str ytit, bool orn);
TGraphAsymmErrors* Draw1D(VecD content, Str var, int Nbins, VecD bins, VecD binlims, int col, int ms, Str xtit, Str ytit, bool orn);

TH1* Draw1D(VecD Values, VecD Errors, VecD Bins, int col=kBlack, int sty=1, int ms=1, Str xtit="", Str ytit="");
TG* DrawTG(VecD Values, VecD Errors, VecD Bins, int col=kBlack, bool diff = false, double offset=0.5, int sty=1, int ms=1, Str xtit="", Str ytit="");
TH2* DrawMatrix(TMatrixDSym mat);

VecD SqrtVecD(VecD Values);

TH1* Draw1D(TTree *t, Str var, Str cut="1",int Nbins=100, double min=-1, double max=1,int col=kBlack, int sty=1, int ms=1, Str xtit="", Str ytit="");
TH1* Draw1D(TFile *f, Str var, Str cut="1", int Nbins=100, double min=-1, double max=1,int col=kBlack, int ms=1, Str xtit="", Str ytit="");
TH1* Draw1D(TTree *t, Str var, Str cut, int Nbins, double min, double max, int col, int sty, int ms, Str xtit, Str ytit);
TH2* Draw2D(TTree *t, Str var, Str cut, int Nbins, double xmin, double xmax, double ymin, double ymax, int col, int ms, Str xtit, Str ytit);
TH2* Draw2D(TTree *t, Str var, Str cut, int col, int ms, Str xtit, Str ytit);
TH1* Draw1D(Str fn, Str hist, Str var, int col, int sty, int ms, Str xtit, Str ytit);
TH1* Draw1D(TH1D *h, Str var, int col, int sty, int ms, Str xtit, Str ytit);
TH1* Draw1D(TTree *t, Str var, Str cut, int Nbins, double *xbins, int col, int sty, int ms, Str xtit, Str ytit);
TH1* Draw1D(TTree *t, Str var, Str cut, VecD xbins, int col, int sty, int ms, Str xtit, Str ytit);

VecD MakeUniformVecD(int N, double min, double max);
TH1F *MakeAxis(VecD bins, Str xtit, Str ytit, double value ,bool sub=false);
TH1F *MakeAxis(int N, double xmin, double xmax, Str xtit, Str ytit, double value = -9999999);

void DrawSubPad();
TH1D* residualHist(const RooHist* rhist, const RooCurve* curve);
double residual( double datum, double pdf);

double DetermineQuantile(TH1* h, double observation);
double GetQuantile(TH1* h, double median, double quantile);
double GetOnesidedQuantile(TH1* h, double quantile);
double DetermineMedian(TH1 *h);

double myProbFunction(double *_x, double *_par);
double myProbFunction2(double *_x, double *_par);
double myProbFunction3(double *_x, double *_par);
double myProbFunction4(double *_x, double *_par);

VecD myVecDivide(VecD one, VecD two, double a = 0., double b = 1., double c = -1.);
VecD myVecDivideError(VecD one, VecD two, VecD oneError, VecD twoError, double a = 0., double b = 1.);

// Simple minimizer to invert interpolated TGraphs

double brent_graph_func(double x);
double DetermineQuantile(TGraph* graph, double quantile, double xmin, double xmax);

double brent_func(double x);
double DetermineQuantile(TF1* tf_func, double quantile, double xmin, double xmax);

double Sum(VecD vec);
VecD Normalize(VecD vec);
VecD NormalizeError(VecD vec, VecD errors);

double KroneckerDelta(int i, int j);


#endif
