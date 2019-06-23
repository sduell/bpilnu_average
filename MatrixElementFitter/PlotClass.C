/*
 *  PlotClass: PlotClass
 */

#include "PlotClass.h"
#include "TPaveText.h"

PlotClass::PlotClass() { }

PlotClass::PlotClass(TEnv *set) {
   
    Printf("> Initializing PlotClass\n");
    
    _set = set;
    
    // We need to replace these
    SetAtlasStyle(); gStyle->SetPalette(1); gStyle->SetHistMinimumZero(); gStyle->SetPaintTextFormat("4.2e");
    
    _can = new TCanvas(); _ps= _set->GetValue("PlotFileName","");
    _can->Print(_ps+"[");

    _file=TFile::Open("output.root","RECREATE");

}


void PlotClass::PlotFitResultMarg(VecD pars, TMatrixDSym Cov, double chi2, int ndf){

	double globfac=1e6;	
	double sum=0., sumerr=0.;

 _theory=new TheoryClass(_set);
 double* x=new double(pars.size());
 TVectorD* param=new TVectorD(pars.size());
 for(int i=0;i<pars.size();i++){x[i]=pars[i]; (*param)(i)=pars[i];};
 _theory->setAllPars(x);

    StrV FitPars=Vectorize(_set->GetValue("FitPars","")," ");

	VecD bins = VectorizeD(_set->GetValue("Bins","")," ");  
            if(bins.size() == 0) 
                Fatal("No binvalues could be found.","");

    VecD low = VectorizeD(_set->GetValue("LowerBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No lower bin edge could be found.","");
    
    VecD high = VectorizeD(_set->GetValue("UpperBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No upper bin edge could be found.",""); 

	VecD values = VectorizeD(_set->GetValue("Data","")," ");  
            if(values.size() == 0) 
                Fatal("No values could be found.","");

Int_t MatDim=atoi(_set->GetValue("MatrixDimension",""));

TMatrixDSym datCov(MatDim);

VecVecD covvalues;

for(int i=0;i<MatDim;i++){
    covvalues.push_back(VectorizeD(_set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
  	for(int j=0;j<MatDim;j++){
        datCov(i,j)=covvalues[i][j];
    	}
    }

TVectorD data(values.size());
TVectorD errs(values.size());
TVectorD bin(values.size());
TVectorD xerr(values.size());

  // eigenvalues: elements of the diagonal form of covariance matrix
  TVectorD lambda;

  // matrix with eigenvector: transformation matrix between original covariance
  // and diagonal form: C' = P^-1 C P

//use this line for the estimated projection for Belle II 50 ab^-1
  //datCov*=0.281;

  TMatrixD P = datCov.EigenVectors(lambda);
  TMatrixD Pinvert = P; Pinvert.Invert();
  TMatrixD diag=Pinvert*datCov*P;



for(int i=0;i<values.size();i++){
	data(i)=values[i];
	errs(i)=sqrt(datCov(i,i));
	//errs(i)=sqrt(diag(i,i));
	bin(i)=bins[i];
	xerr(i)=(high[i]-low[i])/2.;
	sum+=values[i]*(high[i]-low[i]);
	for(int h=0; h<values.size(); h++){
		sumerr+=datCov(i,h);
	}
}
	sumerr=sqrt(sumerr);

	TGraphErrors *graph=new TGraphErrors(bin,data*globfac, xerr, errs*globfac);
/*	graph->SetMinimum(0.);
	graph->SetMaximum(11.);
	graph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
	graph->GetYaxis()->SetTitle("dBF(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
	graph->Draw("AP");*/

int sampling=100;
TVectorD theobin(sampling);
TVectorD theoval(sampling);
TVectorD theo1su(sampling);
TVectorD theo1sd(sampling);
TVectorD theo2su(sampling);
TVectorD theo2sd(sampling);

double* y=new double(*x);
double q2=low[0]+0.000000261; 
//double q2=low[0];

TGraph *grshade = new TGraph(2*sampling-1);
TGraph *gr2shade = new TGraph(2*sampling-1);

for(int i=0;i<sampling;i++){
	theobin(i)=q2;
	theoval(i)=_theory->dBF(q2);
	double temp1su=0;
	double temp1sd=0;
	double temp2su=0;
	double temp2sd=0;
	for(int j=0;j<pars.size();j++){
		this->calcVariatedFF(j,1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};
	}
	theo1su(i)=theoval(i)+sqrt(temp1su);
	theo1sd(i)=theoval(i)-sqrt(temp1sd);
	theo2su(i)=theoval(i)+sqrt(temp2su);
	theo2sd(i)=theoval(i)-sqrt(temp2sd);
	//cout << "testing theovals: "<<q2<<" "<<_theory->dBF(q2)<<endl;
	//cout << theoval(i)<<" "<<theo1su(i)<<" "<<theo1sd(i)<<" "<<theo2su(i)<<" "<<theo2sd(i)<<endl;
	_theory->setAllPars(x);
   	q2+=((high[bins.size()-1]-low[0])/(sampling-1));
}

for(int i=0;i<sampling;i++){
	grshade->SetPoint(i,theobin(i),theo1su(i)*globfac);
    grshade->SetPoint(sampling+i,theobin(sampling-i-1),theo1sd(sampling-i-1)*globfac);
    gr2shade->SetPoint(i,theobin(i),theo2su(i)*globfac);
    gr2shade->SetPoint(sampling+i,theobin(sampling-i-1),theo2sd(sampling-i-1)*globfac);
}

TColor* col=new TColor();

//TGraphAsymmErrors* theograph=new TGraphAsymmErrors(theobin, theoval, xerr, xerr, theo1sd, theo1su);
TGraph* theograph=new TGraph(theobin,theoval*globfac);
theograph->SetLineWidth(2);
theograph->SetLineColor(col->GetColor("#31a354"));

theograph->SetMinimum(0.);
theograph->SetMaximum(11.);
theograph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
//theograph->GetYaxis()->SetTitle("d#it{BF}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
theograph->GetYaxis()->SetTitle("d#it{B}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
//theograph->GetXaxis()->SetRangeUser(-0.5,26.42);
theograph->GetXaxis()->SetLimits(-0.9,26.42);
theograph->Draw("Ac");

double lcsr_fp=0.261;
double z=(sqrt(29.322225-0.)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-0.)+sqrt(29.322225-20.178728518));
double ppi=1./(2.*5.280)*sqrt(pow(pow(5.280,2.)+pow(0.135,2.)-0.,2.)-4.*pow(5.280,2.)*pow(0.135,2.));
double lcsrpoint=globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*lcsr_fp,2.);

double lcsr_errlow=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp-0.023),2.));
double lcsr_errhigh=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp+0.020),2.));

double nul=0.;

TGraphAsymmErrors* lcsr=new TGraphAsymmErrors(1,&nul,&lcsrpoint,&nul,&nul,&lcsr_errlow, &lcsr_errhigh);

lcsr->SetLineColor(kBlue);
lcsr->SetMarkerColor(kBlue);
lcsr->SetMarkerStyle(26);
//lcsr->Draw("P");

gr2shade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.5);
//gr2shade->SetFillColorAlpha(col->GetColor("#e5f5e0"),1.0);
gr2shade->Draw("f");
grshade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.75);
grshade->Draw("f");

//theograph->Draw("Acsame");
TGraph* theoline=new TGraph(theobin,theoval*globfac);
theoline->SetLineWidth(2);
theoline->SetLineColor(col->GetColor("#31a354"));
theoline->Draw("csame");
lcsr->Draw("P");
//graph->SetMarkerStyle(21);
graph->Draw("P");

   //TLegend* leg = new TLegend(0.1,0.7,0.48,0.9,"Fit","C");
   TLegend* leg = new TLegend(0.6,0.8,0.9,0.9);
   leg->SetBorderSize(0);
   leg->AddEntry(graph,"Average Belle + BaBar","leP");
   leg->AddEntry(lcsr,"LCSR (Bharucha)","leP");
   //leg->AddEntry("theograph","Average Belle + BaBar","l");
   leg->Draw();
   //TBox *b = new TBox(17.25,9.25,19.,8.625);
   TBox *b = new TBox(15.25,9.25,17.,8.625);
   //TBox *b = new TBox(15.25,9.125,17.,8.5); 
   //b->SetFillColorAlpha(col->GetColor("#a1d99b"),0.3); 
   b->SetFillColorAlpha(col->GetColor("#e5f5e0"),1.0);
   b->Draw();
   TBox *b2 = new TBox(15.25,9.125,17.,8.75);
   //TBox *b2 = new TBox(15.25,9.,17.,8.625); 
   //TBox *b2 = new TBox(17.25,9.125,19.,8.75);
   //b2->SetFillColorAlpha(col->GetColor("#a1d99b"),0.6); 
   b2->SetFillColorAlpha(col->GetColor("#a1d99b"),1.0); 
   b2->Draw();
   //TLine *lin = new TLine(17.25,8.9375,19.,8.9375);
   //TLine* lin=new TLine(15.25,8.8125,17.,8.8125);
   TLine* lin=new TLine(15.25,8.9375,17.,8.9375);
   //TLine* lin=new TLine(15.25,8.9375,16.975,8.9375);
   lin->SetLineColor(col->GetColor("#31a354"));
   lin->SetLineWidth(2);
   lin->Draw();
   double vubgran=1e3;
   double experror=sqrt(Cov(0,0))*_pull[0];
   //double theoerror=sqrt((Cov(0,0))*(1-pow(_pull[0],2.)));
   double theoerror=sqrt(Cov(0,0)-pow(experror,2.));
   DrawText(Form("    #scale[1.0]{#bf{ BCL fit (3 + 1 parameter)}}"), kBlack,0.779,0.6475);
   DrawText(Form("    #scale[1.0]{#bf{ Data & LQCD (FLAG) & LCSR}}"), kBlack,0.779-(0.75)*0.045,0.6475);
   //DrawText(Form("    #scale[1.0]{#bf{ with constraints from}}"), kBlack,0.779-(0.75)*0.045,0.6475);
   //DrawText(Form("    #scale[1.0]{#bf{ LQCD (FLAG average)}}"), kBlack,0.779-(1.5)*0.045,0.6475);
   //DrawText(Form("    #scale[1.0]{#bf{ and LCSR (Bharucha)}}"), kBlack,0.779-(2.25)*0.045,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ |V_{ub}|= [ %1.2f #pm %1.2f (exp) #pm %1.2f (theo) ] x 10^{-3}}}",pars[0]*vubgran, experror*vubgran,theoerror*vubgran,vubgran), kBlack,0.95-(0.275+2)*0.045,0.14);
   //DrawText(Form("   #scale[1.0]{#it{ #chi^{2} / ndf:  %2.2f / %i}}",chi2, bins.size()+3-ndf-1), kBlack,0.95-(4+2)*0.045,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ Fit prob.:  %3.0f%}}",TMath::Prob(chi2, bins.size()+4-ndf-1)*100.), kBlack,0.95-(1.725+2)*0.045,0.14);

   for(int i=0; i<pars.size(); i++){
   		cout << "parameter "<<i<<": "<<pars[i]<<" +- " << (double)(sqrt(Cov(i,i))*_pull[i]) << "(exp.) +- " << (double)(sqrt(Cov(i,i)*(1-pow(_pull[i],2.))))<< "(theo.)\nTotal Error: "<<(double)(sqrt(Cov(i,i)))<<endl;
   }

//	for(int l=0;l<pars.size();l++){
//		DrawText(Form("-  #scale[1.0]{#it{ %s: %1.2e #pm %1.1e}}",FitPars[l].Data(),pars[l], sqrt(Cov(l,l))), kBlack,0.95-(l+2)*0.045,0.675);
//	}

    //DrawText(Form("-  #scale[1.0]{#it{ average Belle + BaBar  }}"), kBlack,0.95-(0+2)*0.045,0.44);

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
   TBox *box = new TBox(1,1,6,2.5);
   box->SetFillColor(19);
   box->SetFillStyle(0);
   box->SetLineWidth(2);
   box->Draw();
   //box = new TBox(0.1585227,-8.162791,0.2011364,-7.534884);
   box = new TBox(1,1.5,6,2.5);
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
   	cout << "The total branching fraction is: "<<sum <<" +- "<<sumerr<<endl;

	_can->Print(_ps);
	_can->Print("fit.eps");
	_can->Print("fit.gif");
	_can->Print("fit.pdf");
	_can->Write();
	_file->Write();
	_can->Clear();
}

void PlotClass::PlotFitResult(VecD pars, TMatrixDSym Cov, double chi2, int ndf){

	double globfac=1e6;	
	double sum=0., sumerr=0.;

 _theory=new TheoryClass(_set);
 double* x=new double(pars.size());
 TVectorD* param=new TVectorD(pars.size());
 for(int i=0;i<pars.size();i++){x[i]=pars[i]; (*param)(i)=pars[i];};
 _theory->setAllPars(x);

    StrV FitPars=Vectorize(_set->GetValue("FitPars","")," ");

	VecD bins = VectorizeD(_set->GetValue("Bins","")," ");  
            if(bins.size() == 0) 
                Fatal("No binvalues could be found.","");

    VecD low = VectorizeD(_set->GetValue("LowerBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No lower bin edge could be found.","");
    
    VecD high = VectorizeD(_set->GetValue("UpperBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No upper bin edge could be found.",""); 

	VecD values = VectorizeD(_set->GetValue("Data","")," ");  
            if(values.size() == 0) 
                Fatal("No values could be found.","");

Int_t MatDim=atoi(_set->GetValue("MatrixDimension",""));

TMatrixDSym datCov(MatDim);

VecVecD covvalues;

for(int i=0;i<MatDim;i++){
    covvalues.push_back(VectorizeD(_set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
  	for(int j=0;j<MatDim;j++){
        datCov(i,j)=covvalues[i][j];
    	}
    }

TVectorD data(values.size());
TVectorD errs(values.size());
TVectorD bin(values.size());
TVectorD xerr(values.size());

  // eigenvalues: elements of the diagonal form of covariance matrix
  TVectorD lambda;

  // matrix with eigenvector: transformation matrix between original covariance
  // and diagonal form: C' = P^-1 C P

//use this line for the estimated projection for Belle II 50 ab^-1
  //datCov*=0.281;

  TMatrixD P = datCov.EigenVectors(lambda);
  TMatrixD Pinvert = P; Pinvert.Invert();
  TMatrixD diag=Pinvert*datCov*P;



for(int i=0;i<values.size();i++){
	data(i)=values[i];
	errs(i)=sqrt(datCov(i,i));
	//errs(i)=sqrt(diag(i,i));
	bin(i)=bins[i];
	xerr(i)=(high[i]-low[i])/2.;
	sum+=values[i]*(high[i]-low[i]);
	for(int h=0; h<values.size(); h++){
		sumerr+=datCov(i,h);
	}
}
	sumerr=sqrt(sumerr);

	TGraphErrors *graph=new TGraphErrors(bin,data*globfac, xerr, errs*globfac);
/*	graph->SetMinimum(0.);
	graph->SetMaximum(11.);
	graph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
	graph->GetYaxis()->SetTitle("dBF(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
	graph->Draw("AP");*/

int sampling=100;
TVectorD theobin(sampling);
TVectorD theoval(sampling);
TVectorD theo1su(sampling);
TVectorD theo1sd(sampling);
TVectorD theo2su(sampling);
TVectorD theo2sd(sampling);

double* y=new double(*x);
double q2=low[0]+0.000000261; 
//double q2=low[0];

TGraph *grshade = new TGraph(2*sampling-1);
TGraph *gr2shade = new TGraph(2*sampling-1);

for(int i=0;i<sampling;i++){
	theobin(i)=q2;
	theoval(i)=_theory->dBF(q2);
	double temp1su=0;
	double temp1sd=0;
	double temp2su=0;
	double temp2sd=0;
	for(int j=0;j<pars.size();j++){
		this->calcVariatedFF(j,1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};
	}
	theo1su(i)=theoval(i)+sqrt(temp1su);
	theo1sd(i)=theoval(i)-sqrt(temp1sd);
	theo2su(i)=theoval(i)+sqrt(temp2su);
	theo2sd(i)=theoval(i)-sqrt(temp2sd);
	//cout << "testing theovals: "<<q2<<" "<<_theory->dBF(q2)<<endl;
	//cout << theoval(i)<<" "<<theo1su(i)<<" "<<theo1sd(i)<<" "<<theo2su(i)<<" "<<theo2sd(i)<<endl;
	_theory->setAllPars(x);
   	q2+=((high[bins.size()-1]-low[0])/(sampling-1));
}

for(int i=0;i<sampling;i++){
	grshade->SetPoint(i,theobin(i),theo1su(i)*globfac);
    grshade->SetPoint(sampling+i,theobin(sampling-i-1),theo1sd(sampling-i-1)*globfac);
    gr2shade->SetPoint(i,theobin(i),theo2su(i)*globfac);
    gr2shade->SetPoint(sampling+i,theobin(sampling-i-1),theo2sd(sampling-i-1)*globfac);
}

TColor* col=new TColor();

//TGraphAsymmErrors* theograph=new TGraphAsymmErrors(theobin, theoval, xerr, xerr, theo1sd, theo1su);
TGraph* theograph=new TGraph(theobin,theoval*globfac);
theograph->SetLineWidth(2);
theograph->SetLineColor(col->GetColor("#31a354"));

theograph->SetMinimum(0.);
theograph->SetMaximum(11.);
theograph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
//theograph->GetYaxis()->SetTitle("d#it{BF}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
//theograph->GetYaxis()->SetTitle("d#it{B}(B#rightarrow #pi l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
theograph->GetYaxis()->SetTitle("d#it{B}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
//theograph->GetYaxis()->SetTitle("d#it{B}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
//theograph->GetXaxis()->SetRangeUser(-0.5,26.42);
theograph->GetXaxis()->SetLimits(-0.9,26.42);
theograph->Draw("Ac");

double lcsr_fp=0.261;
double z=(sqrt(29.322225-0.)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-0.)+sqrt(29.322225-20.178728518));
double ppi=1./(2.*5.280)*sqrt(pow(pow(5.280,2.)+pow(0.135,2.)-0.,2.)-4.*pow(5.280,2.)*pow(0.135,2.));
double lcsrpoint=globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*lcsr_fp,2.);

double lcsr_errlow=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp-0.023),2.));
double lcsr_errhigh=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp+0.020),2.));

double nul=0.;

TGraphAsymmErrors* lcsr=new TGraphAsymmErrors(1,&nul,&lcsrpoint,&nul,&nul,&lcsr_errlow, &lcsr_errhigh);

lcsr->SetLineColor(kBlue);
lcsr->SetMarkerColor(kBlue);
lcsr->SetMarkerStyle(26);
lcsr->Draw("P");

gr2shade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.5);
gr2shade->Draw("f");
grshade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.75);
grshade->Draw("f");

//graph->SetMarkerStyle(21);
graph->Draw("P");

   //TLegend* leg = new TLegend(0.1,0.7,0.48,0.9,"Fit","C");
   TLegend* leg = new TLegend(0.6,0.8,0.9,0.9);
   leg->SetBorderSize(0);
   leg->AddEntry(graph,"Average Belle + BaBar","leP");
   leg->AddEntry(lcsr,"LCSR: Bharucha","leP");
   //leg->AddEntry("theograph","Average Belle + BaBar","l");
   leg->Draw();
   //TBox *b = new TBox(17.25,9.25,19.,8.625);
   TBox *b = new TBox(15.25,9.25,17.,8.625);
   //TBox *b = new TBox(15.25,9.125,17.,8.5); 
   b->SetFillColorAlpha(col->GetColor("#a1d99b"),0.3); 
   b->Draw();
   TBox *b2 = new TBox(15.25,9.125,17.,8.75);
   //TBox *b2 = new TBox(15.25,9.,17.,8.625); 
   //TBox *b2 = new TBox(17.25,9.125,19.,8.75);
   b2->SetFillColorAlpha(col->GetColor("#a1d99b"),0.6); 
   b2->Draw();
   //TLine *lin = new TLine(17.25,8.9375,19.,8.9375);
   //TLine* lin=new TLine(15.25,8.8125,17.,8.8125);
   TLine* lin=new TLine(15.25,8.9375,17.,8.9375);
   lin->SetLineColor(col->GetColor("#31a354"));
   lin->SetLineWidth(2);
   lin->Draw();
   double vubgran=1e3;
   DrawText(Form("    #scale[1.0]{#bf{ BCL ( 3 + 1 parameter )}}"), kBlack,0.779,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ |V_{ub}|: ( %1.2f #pm %1.2f ) x 10^{-3}}}",pars[0]*vubgran, sqrt(Cov(0,0))*vubgran,vubgran), kBlack,0.95-(0.275+2)*0.045,0.14);
   //DrawText(Form("   #scale[1.0]{#it{ #chi^{2} / ndf:  %2.2f / %i}}",chi2, bins.size()+3-ndf-1), kBlack,0.95-(4+2)*0.045,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ Fit prob.:  %3.0f%}}",TMath::Prob(chi2, bins.size()+4-ndf-1)*100.), kBlack,0.95-(1.725+2)*0.045,0.14);


//	for(int l=0;l<pars.size();l++){
//		DrawText(Form("-  #scale[1.0]{#it{ %s: %1.2e #pm %1.1e}}",FitPars[l].Data(),pars[l], sqrt(Cov(l,l))), kBlack,0.95-(l+2)*0.045,0.675);
//	}

    //DrawText(Form("-  #scale[1.0]{#it{ average Belle + BaBar  }}"), kBlack,0.95-(0+2)*0.045,0.44);

//this is the HFAG box
   //TBox *box = new TBox(0.1585227,-8.476744,0.2011364,-7.534884);
   TBox *box = new TBox(1,1,6,2.5);
   box->SetFillColor(19);
   box->SetFillStyle(0);
   box->SetLineWidth(2);
   box->Draw();
   //box = new TBox(0.1585227,-8.162791,0.2011364,-7.534884);
   box = new TBox(1,1.5,6,2.5);
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
   text = new TText(3.5,1.3,"Summer 2016");
   text->SetTextAlign(22);
   text->SetTextFont(53);
   text->SetTextSize(15);
   text->Draw();

   	cout << "The total branching fraction is: "<<sum <<" +- "<<sumerr<<endl;



	_can->Print(_ps);
	_can->Clear();
}


void PlotClass::PlotFitResult(VecD pars, TMatrixDSym Cov, double chi2, int ndf, double exerr, double theoerr){

	double globfac=1e6;
	//double globfac=1.;	

 _theory=new TheoryClass(_set);
 double* x=new double(pars.size());
 TVectorD* param=new TVectorD(pars.size());
 for(int i=0;i<pars.size();i++){x[i]=pars[i]; (*param)(i)=pars[i];};
 _theory->setAllPars(x);

    StrV FitPars=Vectorize(_set->GetValue("FitPars","")," ");

	VecD bins = VectorizeD(_set->GetValue("Bins","")," ");  
            if(bins.size() == 0) 
                Fatal("No binvalues could be found.","");

    VecD low = VectorizeD(_set->GetValue("LowerBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No lower bin edge could be found.","");
    
    VecD high = VectorizeD(_set->GetValue("UpperBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No upper bin edge could be found.",""); 

	VecD values = VectorizeD(_set->GetValue("Data","")," ");  
            if(values.size() == 0) 
                Fatal("No values could be found.","");

Int_t MatDim=atoi(_set->GetValue("MatrixDimension",""));

TMatrixDSym datCov(MatDim);

VecVecD covvalues;

for(int i=0;i<MatDim;i++){
    covvalues.push_back(VectorizeD(_set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
  	for(int j=0;j<MatDim;j++){
        datCov(i,j)=covvalues[i][j];
    	}
    }

TVectorD data(values.size());
TVectorD errs(values.size());
TVectorD bin(values.size());
TVectorD xerr(values.size());

for(int i=0;i<values.size();i++){
	data(i)=values[i];
	errs(i)=sqrt(datCov(i,i));
	bin(i)=bins[i];
	xerr(i)=(high[i]-low[i])/2.;
}

	TGraphErrors *graph=new TGraphErrors(bin,data*globfac, xerr, errs*globfac);
/*	graph->SetMinimum(0.);
	graph->SetMaximum(11.);
	graph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
	graph->GetYaxis()->SetTitle("dBF(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
	graph->Draw("AP");*/

int sampling=100;
TVectorD theobin(sampling);
TVectorD theoval(sampling);
TVectorD theo1su(sampling);
TVectorD theo1sd(sampling);
TVectorD theo2su(sampling);
TVectorD theo2sd(sampling);

double* y=new double(*x);
double q2=low[0]+0.000000261; 

TGraph *grshade = new TGraph(2*sampling-1);
TGraph *gr2shade = new TGraph(2*sampling-1);

for(int i=0;i<sampling;i++){
	theobin(i)=q2;
	theoval(i)=_theory->dBF(q2);
	double temp1su=0;
	double temp1sd=0;
	double temp2su=0;
	double temp2sd=0;
	for(int j=0;j<pars.size();j++){
		this->calcVariatedFF(j,1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};
	}
	theo1su(i)=theoval(i)+sqrt(temp1su);
	theo1sd(i)=theoval(i)-sqrt(temp1sd);
	theo2su(i)=theoval(i)+sqrt(temp2su);
	theo2sd(i)=theoval(i)-sqrt(temp2sd);
	//cout << "testing theovals: "<<q2<<" "<<_theory->dBF(q2)<<endl;
	//cout << theoval(i)<<" "<<theo1su(i)<<" "<<theo1sd(i)<<" "<<theo2su(i)<<" "<<theo2sd(i)<<endl;
	_theory->setAllPars(x);
   	q2+=((high[bins.size()-1]-low[0])/(sampling-1));
}

for(int i=0;i<sampling;i++){
	grshade->SetPoint(i,theobin(i),theo1su(i)*globfac);
    grshade->SetPoint(sampling+i,theobin(sampling-i-1),theo1sd(sampling-i-1)*globfac);
    gr2shade->SetPoint(i,theobin(i),theo2su(i)*globfac);
    gr2shade->SetPoint(sampling+i,theobin(sampling-i-1),theo2sd(sampling-i-1)*globfac);
}

TColor* col=new TColor();

//TGraphAsymmErrors* theograph=new TGraphAsymmErrors(theobin, theoval, xerr, xerr, theo1sd, theo1su);
TGraph* theograph=new TGraph(theobin,theoval*globfac);
theograph->SetLineWidth(2);
theograph->SetLineColor(col->GetColor("#31a354"));

theograph->SetMinimum(0.);
theograph->SetMaximum(11.);
theograph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
//theograph->GetYaxis()->SetTitle("d#it{BF}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
theograph->GetYaxis()->SetTitle("d#it{B}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
//theograph->GetXaxis()->SetRangeUser(-0.5,26.42);
theograph->GetXaxis()->SetLimits(-0.9,26.42);
theograph->Draw("Ac");

double lcsr_fp=0.261;
double z=(sqrt(29.322225-0.)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-0.)+sqrt(29.322225-20.178728518));
double ppi=1./(2.*5.280)*sqrt(pow(pow(5.280,2.)+pow(0.135,2.)-0.,2.)-4.*pow(5.280,2.)*pow(0.135,2.));
double lcsrpoint=globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*lcsr_fp,2.);

double lcsr_errlow=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp-0.023),2.));
double lcsr_errhigh=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp+0.020),2.));

double nul=0.;

TGraphAsymmErrors* lcsr=new TGraphAsymmErrors(1,&nul,&lcsrpoint,&nul,&nul,&lcsr_errlow, &lcsr_errhigh);

lcsr->SetLineColor(kBlue);
lcsr->SetMarkerColor(kBlue);
lcsr->SetMarkerStyle(26);
//lcsr->Draw("P");

gr2shade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.5);
gr2shade->Draw("f");
grshade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.75);
grshade->Draw("f");

//graph->SetMarkerStyle(21);
graph->Draw("P");

   //TLegend* leg = new TLegend(0.1,0.7,0.48,0.9,"Fit","C");
   TLegend* leg = new TLegend(0.6,0.8,0.9,0.9);
   leg->SetBorderSize(0);
   leg->AddEntry(graph,"Average Belle + BaBar","leP");
   //leg->AddEntry(lcsr,"LCSR: Bharucha","leP");
   //leg->AddEntry("theograph","Average Belle + BaBar","l");
   leg->Draw();
   //TBox *b = new TBox(17.25,9.25,19.,8.625);
   TBox *b = new TBox(15.25,9.25,17.,8.625);
   //TBox *b = new TBox(15.25,9.125,17.,8.5); 
   b->SetFillColorAlpha(col->GetColor("#a1d99b"),0.3); 
   b->Draw();
   TBox *b2 = new TBox(15.25,9.125,17.,8.75);
   //TBox *b2 = new TBox(15.25,9.,17.,8.625); 
   //TBox *b2 = new TBox(17.25,9.125,19.,8.75);
   b2->SetFillColorAlpha(col->GetColor("#a1d99b"),0.6); 
   b2->Draw();
   //TLine *lin = new TLine(17.25,8.9375,19.,8.9375);
   //TLine* lin=new TLine(15.25,8.8125,17.,8.8125);
   TLine* lin=new TLine(15.25,8.9375,17.,8.9375);
   lin->SetLineColor(col->GetColor("#31a354"));
   lin->SetLineWidth(2);
   lin->Draw();
   double vubgran=1e3;
   DrawText(Form("    #scale[1.0]{#bf{ BCL ( 3 + 1 parameter )}}"), kBlack,0.772,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ |V_{ub}|: ( %1.2f #pm %1.2f (exp.) #pm %1.2f (theo.)) x 10^{-3}}}",pars[0]*vubgran, exerr*vubgran,theoerr*vubgran), kBlack,0.95-(0.275+2)*0.045,0.14);
   //DrawText(Form("   #scale[1.0]{#bf{ |V_{ub}|: ( %1.2f #pm %1.2f ) x 10^{-3}}}",pars[0]*vubgran, sqrt(Cov(0,0))*vubgran,vubgran), kBlack,0.95-(0.275+2)*0.045,0.14);
   //DrawText(Form("   #scale[1.0]{#it{ #chi^{2} / ndf:  %2.2f / %i}}",chi2, bins.size()+3-ndf), kBlack,0.95-(4+2)*0.045,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ Fit prob.:  %3.0f%}}",TMath::Prob(chi2, bins.size()+4-ndf)*100.), kBlack,0.95-(1.725+2)*0.045,0.14);
   //DrawText(Form("   #scale[1.0]{#bf{ Fit prob.:  %3.0f%}}",TMath::Prob(chi2, bins.size()+4-ndf-1)*100.), kBlack,0.95-(1.725+2)*0.045,0.14);


//	for(int l=0;l<pars.size();l++){
//		DrawText(Form("-  #scale[1.0]{#it{ %s: %1.2e #pm %1.1e}}",FitPars[l].Data(),pars[l], sqrt(Cov(l,l))), kBlack,0.95-(l+2)*0.045,0.675);
//	}

    //DrawText(Form("-  #scale[1.0]{#it{ average Belle + BaBar  }}"), kBlack,0.95-(0+2)*0.045,0.44);
//this is the HFAG box
   //TBox *box = new TBox(0.1585227,-8.476744,0.2011364,-7.534884);
   TBox *box = new TBox(1,1,6,2.5);
   box->SetFillColor(19);
   box->SetFillStyle(0);
   box->SetLineWidth(2);
   box->Draw();
   //box = new TBox(0.1585227,-8.162791,0.2011364,-7.534884);
   box = new TBox(1,1.5,6,2.5);
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
   text = new TText(3.5,1.3,"Summer 2016");
   text->SetTextAlign(22);
   text->SetTextFont(53);
   text->SetTextSize(15);
   text->Draw();

	_can->Print(_ps);
	_can->Clear();
}

void PlotClass::calcVariatedFF(int variedPar, double NsigmaFF, TVectorD* FFval, TMatrixDSym Cov){

  // Covariance matrix of FF using world avarages
  TMatrixD FFcov=Cov;

  // eigenvalues: elements of the diagonal form of covariance matrix
  TVectorD lambda;

  // matrix with eigenvector: transformation matrix between original covariance
  // and diagonal form: C' = P^-1 C P
  TMatrixD P = FFcov.EigenVectors(lambda);
  TMatrixD Pinvert = P; Pinvert.Invert();

  // Transform FF parameters in new basis
  (*FFval) *= Pinvert;

  // Vary new FF parameters inside errors given by sqrt of eigenvalues
  (*FFval)(abs(variedPar)) = (*FFval)(abs(variedPar)) + NsigmaFF*sqrt(lambda(abs(variedPar)));
  // Rotate back into original basis to apply event weights
  (*FFval) *= P;
}

void PlotClass::calcVariatedFF(int variedPar, double NsigmaFF){

  // Vector containing world averages of FF parameters
  TVectorD FFval(3);

  FFval(0) = 0.421;
  FFval(1) = -0.35;
  FFval(2) = -0.41;
  
if(abs(variedPar<4))
{

  // Covariance matrix of FF using world avarages
  TMatrixD FFcov(3,3);

//0.000169, 0.0003978, 0.00069888, 0.0003978, 0.01,0.054784, 0.00069888,0.054784,0.4096

  FFcov(0,0) = 0.000169;    FFcov(0,1) = 0.0003978;     FFcov(0,2) = 0.00069888;
  FFcov(1,0) = 0.0003978;    FFcov(1,1) = 0.01;    FFcov(1,2) = 0.054784;
  FFcov(2,0) = 0.00069888;   FFcov(2,1) = 0.054784;    FFcov(2,2) = 0.4096;
  
  // eigenvalues: elements of the diagonal form of covariance matrix
  TVectorD lambda;

  // matrix with eigenvector: transformation matrix between original covariance
  // and diagonal form: C' = P^-1 C P
  TMatrixD P = FFcov.EigenVectors(lambda);
  TMatrixD Pinvert = P; Pinvert.Invert();

  // Transform FF parameters in new basis
  FFval *= Pinvert;

  // Vary new FF parameters inside errors given by sqrt of eigenvalues
  FFval(abs(variedPar)-1) = FFval(abs(variedPar)-1) + NsigmaFF*sqrt(lambda(abs(variedPar)-1));

  // Rotate back into original basis to apply event weights
  FFval *= P;
} 

else cout << "more than 3 is not atm supported, nothing will be changed."<<endl;
  //TODO: encode FF in a vector so the manipulation of only 1 entry is possible!
}

void PlotClass::Drawpval(VecD chi2){
		//gStyle->SetOptStat(2201);
		TH1D* hist = new TH1D("pval","pval", chi2.size()/100.+1, 0,1);
		for(int i=0;i<chi2.size();i++){
			//hist->Fill(TMath::Prob(chi2[i], 13+4-4));
			//hist->Fill(TMath::Prob(chi2[i], 13+4-4-1));
			hist->Fill(TMath::Prob(chi2[i], 13+4-4-1-1));
		}
		hist->Draw();
		_can->Print(_ps);
		_can->Clear();
		delete hist;
	}

void PlotClass::Plotfpvalues(double mean, VecD pars){
		gStyle->SetOptStat(2201);
		TH1D* hist = new TH1D("fp0","fp0", pars.size()/100.+1, mean-0.1,mean+0.1);
		for(int i=0;i<pars.size();i++){
			hist->Fill(pars[i]);
		}
		TLine mid(mean,0.,mean,5000.);
		TLine bha(0.261,0.,0.261,5000.);
		mid.SetLineColor(kGreen);
		bha.SetLineColor(kBlue);
		hist->Draw();
		mid.Draw("same");
		bha.Draw("same");
		_can->Print(_ps);
		_can->Clear();
		delete hist;
	}

void PlotClass::PlotToyPulls(VecD pars, VecVecD parerrs, VecVecD toypars){
	for(int i=0; i<pars.size(); i++){
		gStyle->SetOptStat(2201);
		TH1D* hist = new TH1D("hist","hist", toypars.size()/100.+1, -3.5,3.5);
		TH1D* hist2 = new TH1D("hist2","hist2", toypars.size()/100.+1, -0.003, 0.003);
		for(int j=0; j< toypars.size(); j++){
			hist->Fill((pars[i]-toypars[j][i])/parerrs[j][i]);
			hist2->Fill(pars[i]-toypars[j][i]);
		}
		hist2->Draw();
		_can->Print(_ps);
		_can->Clear();
		hist->Draw();
		_can->Print(_ps);
		_can->Clear();
		//cout << "plotted toy pull, fill hist with "<<hist->GetStdDev()<<endl;
		_pull.push_back(hist->GetStdDev());
		//cout << "filled hist with "<<_pull[i]<<endl;
		delete hist;
		delete hist2;
	}
}

void PlotClass::PlotToyPulls(VecD pars, VecD parerrs, VecVecD toypars){
	for(int i=0; i<pars.size(); i++){
		gStyle->SetOptStat(2201);
		TH1D* hist = new TH1D("hist","hist", toypars.size()/100.+1, -3.5,3.5);
		TH1D* hist2 = new TH1D("hist2","hist2", toypars.size()/100.+1, -0.003, 0.003);
		for(int j=0; j< toypars.size(); j++){
			hist->Fill((pars[i]-toypars[j][i])/parerrs[i]);
			hist2->Fill(pars[i]-toypars[j][i]);
		}
		hist2->Draw();
		_can->Print(_ps);
		_can->Clear();
		hist->Draw();
		_can->Print(_ps);
		_can->Clear();
		_pull.push_back(hist->GetStdDev());
		delete hist;
		delete hist2;
	}
}

void PlotClass::PlotToyPulls(VecD pars, TMatrixDSym* Parcov, VecVecD toypars){
	TVectorD lmbd;
	TMatrixD P = Parcov->EigenVectors(lmbd);
	for(int i=0; i<pars.size(); i++){
		gStyle->SetOptStat(2201);
		TH1D* hist = new TH1D("hist","hist", toypars.size()/100.+1, -3.5,3.5);
		TH1D* hist2;
		if(i==0) hist2 = new TH1D("hist2","hist2", toypars.size()/100.+1, -0.003, 0.003);
		else hist2 = new TH1D("hist2","hist2", toypars.size()/100.+1, -0.3, 0.3);
		for(int j=0; j< toypars.size(); j++){
			hist->Fill((pars[i]-toypars[j][i])/sqrt(lmbd(i)));
			hist2->Fill(pars[i]-toypars[j][i]);
		}
		hist2->Draw();
		_can->Print(_ps);
		_can->Clear();
		hist->Draw();
		_can->Print(_ps);
		_can->Clear();
		_pull.push_back(hist->GetStdDev());
		delete hist;
		delete hist2;
	}
}

void PlotClass::PlotToyResult(TVectorD datvec, VecD pars, TMatrixDSym Cov, double chi2, int ndf){
cout << "Plotting toy"<<endl;
	double globfac=1e6;	
	double sum=0., sumerr=0.;

 _theory=new TheoryClass(_set);
 double* x=new double(pars.size());
 TVectorD* param=new TVectorD(pars.size());
 for(int i=0;i<pars.size();i++){x[i]=pars[i]; (*param)(i)=pars[i];};
 _theory->setAllPars(x);

    StrV FitPars=Vectorize(_set->GetValue("FitPars","")," ");

	VecD bins = VectorizeD(_set->GetValue("Bins","")," ");  
            if(bins.size() == 0) 
                Fatal("No binvalues could be found.","");

    VecD low = VectorizeD(_set->GetValue("LowerBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No lower bin edge could be found.","");
    
    VecD high = VectorizeD(_set->GetValue("UpperBinEdge","")," ");  
            if(bins.size() == 0) 
                Fatal("No upper bin edge could be found.",""); 

	VecD values = VectorizeD(_set->GetValue("Data","")," ");  
            if(values.size() == 0) 
                Fatal("No values could be found.","");

Int_t MatDim=atoi(_set->GetValue("MatrixDimension",""));

TMatrixDSym datCov(MatDim);

VecVecD covvalues;

for(int i=0;i<MatDim;i++){
    covvalues.push_back(VectorizeD(_set->GetValue(Form("CovarianceMatrix.row%i",i+1),"")," ")); 
  	for(int j=0;j<MatDim;j++){
        datCov(i,j)=covvalues[i][j];
    	}
    }

TVectorD data(values.size());
TVectorD errs(values.size());
TVectorD bin(values.size());
TVectorD xerr(values.size());

  // eigenvalues: elements of the diagonal form of covariance matrix
  TVectorD lambda;

  // matrix with eigenvector: transformation matrix between original covariance
  // and diagonal form: C' = P^-1 C P

//use this line for the estimated projection for Belle II 50 ab^-1
  //datCov*=0.281;

  TMatrixD P = datCov.EigenVectors(lambda);
  TMatrixD Pinvert = P; Pinvert.Invert();
  TMatrixD diag=Pinvert*datCov*P;



for(int i=0;i<values.size();i++){
	data(i)=datvec(i);
	//errs(i)=sqrt(datCov(i,i));
	errs(i)=sqrt(diag(i,i));
	bin(i)=bins[i];
	xerr(i)=(high[i]-low[i])/2.;
	sum+=values[i]*(high[i]-low[i]);
	for(int h=0; h<values.size(); h++){
		sumerr+=datCov(i,h);
	}
}
	sumerr=sqrt(sumerr);

	TGraphErrors *graph=new TGraphErrors(bin,data*globfac, xerr, errs*globfac);
/*	graph->SetMinimum(0.);
	graph->SetMaximum(11.);
	graph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
	graph->GetYaxis()->SetTitle("dBF(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
	graph->Draw("AP");*/

int sampling=100;
TVectorD theobin(sampling);
TVectorD theoval(sampling);
TVectorD theo1su(sampling);
TVectorD theo1sd(sampling);
TVectorD theo2su(sampling);
TVectorD theo2sd(sampling);

double* y=new double(*x);
double q2=low[0]+0.000000261; 
//double q2=low[0];

TGraph *grshade = new TGraph(2*sampling-1);
TGraph *gr2shade = new TGraph(2*sampling-1);

for(int i=0;i<sampling;i++){
	theobin(i)=q2;
	theoval(i)=_theory->dBF(q2);
	double temp1su=0;
	double temp1sd=0;
	double temp2su=0;
	double temp2sd=0;
	for(int j=0;j<pars.size();j++){
		this->calcVariatedFF(j,1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-1,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp1sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2su+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};

		this->calcVariatedFF(j,-2,param,Cov);
		for(int k=0;k<pars.size();k++){y[k]=(*param)(k);};
		_theory->setAllPars(y);
		temp2sd+=pow(theoval(i)-_theory->dBF(q2),2.);
		for(int k=0;k<pars.size();k++){(*param)(k)=x[k];};
	}
	theo1su(i)=theoval(i)+sqrt(temp1su);
	theo1sd(i)=theoval(i)-sqrt(temp1sd);
	theo2su(i)=theoval(i)+sqrt(temp2su);
	theo2sd(i)=theoval(i)-sqrt(temp2sd);
	//cout << "testing theovals: "<<q2<<" "<<_theory->dBF(q2)<<endl;
	//cout << theoval(i)<<" "<<theo1su(i)<<" "<<theo1sd(i)<<" "<<theo2su(i)<<" "<<theo2sd(i)<<endl;
	_theory->setAllPars(x);
   	q2+=((high[bins.size()-1]-low[0])/(sampling-1));
}

for(int i=0;i<sampling;i++){
	grshade->SetPoint(i,theobin(i),theo1su(i)*globfac);
    grshade->SetPoint(sampling+i,theobin(sampling-i-1),theo1sd(sampling-i-1)*globfac);
    gr2shade->SetPoint(i,theobin(i),theo2su(i)*globfac);
    gr2shade->SetPoint(sampling+i,theobin(sampling-i-1),theo2sd(sampling-i-1)*globfac);
}

TColor* col=new TColor();

//TGraphAsymmErrors* theograph=new TGraphAsymmErrors(theobin, theoval, xerr, xerr, theo1sd, theo1su);
TGraph* theograph=new TGraph(theobin,theoval*globfac);
theograph->SetLineWidth(2);
theograph->SetLineColor(col->GetColor("#31a354"));

theograph->SetMinimum(0.);
theograph->SetMaximum(11.);
theograph->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
theograph->GetYaxis()->SetTitle("d#it{BF}(B^{0}#rightarrow #pi^{-} l^{+} #nu_{l})/dq^{2} [10^{6} GeV^{-2}]");
//theograph->GetXaxis()->SetRangeUser(-0.5,26.42);
theograph->GetXaxis()->SetLimits(-0.9,26.42);
theograph->Draw("Ac");

double lcsr_fp=0.261;
double z=(sqrt(29.322225-0.)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-0.)+sqrt(29.322225-20.178728518));
double ppi=1./(2.*5.280)*sqrt(pow(pow(5.280,2.)+pow(0.135,2.)-0.,2.)-4.*pow(5.280,2.)*pow(0.135,2.));
double lcsrpoint=globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*lcsr_fp,2.);

double lcsr_errlow=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp-0.023),2.));
double lcsr_errhigh=abs(lcsrpoint-globfac*1.52E-12/6.58211899E-25*8*ppi/3.*1.3604189769E-10*pow(pars[0],2.)/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*ppi*(lcsr_fp+0.020),2.));

double nul=0.;

TGraphAsymmErrors* lcsr=new TGraphAsymmErrors(1,&nul,&lcsrpoint,&nul,&nul,&lcsr_errlow, &lcsr_errhigh);

lcsr->SetLineColor(kBlue);
lcsr->SetMarkerColor(kBlue);
lcsr->SetMarkerStyle(26);
lcsr->Draw("P");

gr2shade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.5);
gr2shade->Draw("f");
grshade->SetFillColorAlpha(col->GetColor("#a1d99b"),0.75);
grshade->Draw("f");

//graph->SetMarkerStyle(21);
graph->Draw("P");

   //TLegend* leg = new TLegend(0.1,0.7,0.48,0.9,"Fit","C");
   TLegend* leg = new TLegend(0.6,0.8,0.9,0.9);
   leg->SetBorderSize(0);
   leg->AddEntry(graph,"Average Belle + BaBar","leP");
   leg->AddEntry(lcsr,"LCSR: Bharucha","leP");
   //leg->AddEntry("theograph","Average Belle + BaBar","l");
   leg->Draw();
   //TBox *b = new TBox(17.25,9.25,19.,8.625);
   TBox *b = new TBox(15.25,9.25,17.,8.625);
   //TBox *b = new TBox(15.25,9.125,17.,8.5); 
   b->SetFillColorAlpha(col->GetColor("#a1d99b"),0.3); 
   b->Draw();
   TBox *b2 = new TBox(15.25,9.125,17.,8.75);
   //TBox *b2 = new TBox(15.25,9.,17.,8.625); 
   //TBox *b2 = new TBox(17.25,9.125,19.,8.75);
   b2->SetFillColorAlpha(col->GetColor("#a1d99b"),0.6); 
   b2->Draw();
   //TLine *lin = new TLine(17.25,8.9375,19.,8.9375);
   //TLine* lin=new TLine(15.25,8.8125,17.,8.8125);
   TLine* lin=new TLine(15.25,8.9375,17.,8.9375);
   lin->SetLineColor(col->GetColor("#31a354"));
   lin->SetLineWidth(2);
   lin->Draw();
   double vubgran=1e3;
   DrawText(Form("    #scale[1.0]{#bf{ BCL ( 3 + 1 parameter )}}"), kBlack,0.779,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ |V_{ub}|: ( %1.2f #pm %1.2f ) x 10^{-3}}}",pars[0]*vubgran, sqrt(Cov(0,0))*vubgran,vubgran), kBlack,0.95-(0.275+2)*0.045,0.14);
   //DrawText(Form("   #scale[1.0]{#it{ #chi^{2} / ndf:  %2.2f / %i}}",chi2, bins.size()+3-ndf-1), kBlack,0.95-(4+2)*0.045,0.6475);
   DrawText(Form("   #scale[1.0]{#bf{ Fit prob.:  %3.0f%}}",TMath::Prob(chi2, bins.size()+4-ndf-1)*100.), kBlack,0.95-(1.725+2)*0.045,0.14);


//	for(int l=0;l<pars.size();l++){
//		DrawText(Form("-  #scale[1.0]{#it{ %s: %1.2e #pm %1.1e}}",FitPars[l].Data(),pars[l], sqrt(Cov(l,l))), kBlack,0.95-(l+2)*0.045,0.675);
//	}

    //DrawText(Form("-  #scale[1.0]{#it{ average Belle + BaBar  }}"), kBlack,0.95-(0+2)*0.045,0.44);

   	cout << "The total branching fraction is: "<<sum <<" +- "<<sumerr<<endl;

	_can->Print(_ps);
	_can->Clear();
}

void PlotClass::PlotFitSteps(TVectorD bin, TVectorD data,TVectorD theo,TVectorD errvec,TVectorD erryvec, double *x, double chi2){
	TFormula *zpar=new TFormula("zpar","(sqrt(29.322225-x)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-x)+sqrt(29.322225-20.178728518))");
	TFormula *fp=new TFormula("fp","1./(1-(x/pow(5.325,2.)))*(([0]*pow(zpar,0.)-pow(-1,-4.)*pow(zpar,4.)*0./4.)+([1]*pow(zpar,1.)-pow(-1,-3.)*pow(zpar,4.)*1./4.)+([2]*pow(zpar,2.)-pow(-1,-2.)*pow(zpar,4.)*2./4.)+([3]*pow(zpar,3.)-pow(-1,-1.)*pow(zpar,4.)*3./4.))");
	fp->SetParameter(0,x[1]);
	fp->SetParameter(1,x[2]);
	fp->SetParameter(2,x[3]);
	fp->SetParameter(3,x[4]);
	TFormula *pX=new TFormula("pX","1./(2.*5.280)*sqrt(pow(pow(5.280,2.)+pow(0.135,2.)-x,2.)-4.*pow(5.280,2.)*pow(0.135,2.))");
	TFormula *form=new TFormula("form","1.52E-12/6.58211899E-25*8*pX/3.*1.3604189769E-10*pow([Vub],2.)*x/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*pX*fp/sqrt(x),2.)");
	form->SetParameter("Vub",x[0]);
	TF1 *func=new TF1("func","form",0,26.4);
	//TF1 *func=new TF1("func","1.52E-12/6.58211899E-25*8*pX/3.*1.3604189769E-10*pow([Vub],2.)*x/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*pX*fp/sqrt(x),2.)",0.,26.4);
	func->SetParameter(0,x[1]);
	func->SetParameter(1,x[2]);
	func->SetParameter(2,x[3]);
	func->SetParameter(3,x[4]);
	func->SetParameter("Vub",x[0]);
//	TF1 *fun1=new TF1("test","3E-6",0,26.4);
	TGraphErrors *graph=new TGraphErrors(bin,data, erryvec, errvec);
	TGraph *theograph=new TGraph(bin,theo);

	
	graph->Draw("AP*");
	theograph->Draw("same");
	//func->Draw("same");
	
	for(int l=0;l<5;l++){
		DrawText(Form("+  #scale[1.0]{#it{ Parameter %i: %f}}",l,x[l]), kBlack,0.875-(l+2)*0.045,0.69);
	}
	DrawText(Form("+  #scale[1.0]{#it{ Chi2: %f}}",chi2), kBlack,0.875-(5+2)*0.045,0.69);
//	fun1->Draw("same");
	_can->Print(_ps);
	_can->Clear();
}

void PlotClass::PlotTheoryPred(TVectorD bin, double *x){
	TFormula *zpar=new TFormula("zpar","(sqrt(29.322225-x)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-x)+sqrt(29.322225-20.178728518))");
	TF1 *fp=new TF1("fp","1./(1-(x/pow(5.325,2.)))*(([0]*pow(zpar,0.)-pow(-1,-4.)*pow(zpar,4.)*0./4.)+([1]*pow(zpar,1.)-pow(-1,-3.)*pow(zpar,4.)*1./4.)+([2]*pow(zpar,2.)-pow(-1,-2.)*pow(zpar,4.)*2./4.)+([3]*pow(zpar,3.)-pow(-1,-1.)*pow(zpar,4.)*3./4.))",bin(0),bin(12));
	fp->SetParameter(0,x[1]);
	fp->SetParameter(1,x[2]);
	fp->SetParameter(2,x[3]);
	fp->SetParameter(3,x[4]);

	TF1 *z=new TF1("z","zpar",bin(0),bin(12));

	TheoryClass* _theory=new TheoryClass(_set);
	_theory->setAllPars(x);

	TVectorD pred(13);
	TVectorD zpred(13);

	for(int i=0;i<13;i++){
		  pred(i)=_theory->fp(bin(i));
		  zpred(i)=_theory->z(bin(i));
	}

	TGraph *graph=new TGraph(bin,pred);

	graph->SetMaximum( 8.5 );
	graph->Draw("AC*");
	fp->Draw("same");
	//fp->Draw();


	for(int l=0;l<5;l++){
		DrawText(Form("+  #scale[1.0]{#it{ Parameter %i: %f}}",l,x[l]), kBlack,0.875-(l+2)*0.045,0.69);
	}
//	fun1->Draw("same");
	_can->Print(_ps);
	_can->Clear();

	TGraph *zgraph=new TGraph(bin,zpred);
	zgraph->Draw("AC*");
	z->Draw("same");



	for(int l=0;l<5;l++){
		DrawText(Form("+  #scale[1.0]{#it{ Parameter %i: %f}}",l,x[l]), kBlack,0.875-(l+2)*0.045,0.69);
	}
//	fun1->Draw("same");
	_can->Print(_ps);
	_can->Clear();
}

void PlotClass::PlotRatePred(TVectorD bin, double *x){
	cout <<"x0 "<<x[0]<<endl;
	TFormula *zpar=new TFormula("zpar","(sqrt(29.322225-x)-sqrt(29.322225-20.178728518))/(sqrt(29.322225-x)+sqrt(29.322225-20.178728518))");
	TFormula *fp=new TFormula("fp","1./(1-(x/pow(5.325,2.)))*(([0]*pow(zpar,0.)-pow(-1,-4.)*pow(zpar,4.)*0./4.)+([1]*pow(zpar,1.)-pow(-1,-3.)*pow(zpar,4.)*1./4.)+([2]*pow(zpar,2.)-pow(-1,-2.)*pow(zpar,4.)*2./4.)+([3]*pow(zpar,3.)-pow(-1,-1.)*pow(zpar,4.)*3./4.))");
	fp->SetParameter(0,x[1]);
	fp->SetParameter(1,x[2]);
	fp->SetParameter(2,x[3]);
	fp->SetParameter(3,x[4]);
	TFormula *pX=new TFormula("pX","1./(2.*5.280)*sqrt(pow(pow(5.280,2.)+pow(0.135,2.)-x,2.)-4.*pow(5.280,2.)*pow(0.135,2.))");
	TFormula *form=new TFormula("form","8*pX/3.*1.3604189769E-10*pow([Vub],2.)*x/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*pX*fp/sqrt(x),2.)");
	form->SetParameter("Vub",x[0]);
	//TF1 *func=new TF1("func","8*pX/3.*1.3604189769E-10*pow([4],2.)*x/(256*pow(TMath::Pi(),3.)*pow(5.280,2.))*pow(2.*5.280*pX*fp/sqrt(x),2.)",0.,26.4);
	//func->SetParameter(0,x[1]);
	//func->SetParameter(1,x[2]);
	//func->SetParameter(2,x[3]);
	//func->SetParameter(3,x[4]);
	//func->SetParameter(4,0.00000001);

	TF1 *func = new TF1("func","form",0,26.4);	

	func->SetParameter(0,x[1]);
	func->SetParameter(1,x[2]);
	func->SetParameter(2,x[3]);
	func->SetParameter(3,x[4]);
	func->SetParameter("Vub",x[0]);

	//TF1 *z=new TF1("z","zpar",bin(0),bin(12));

	TheoryClass* _theory=new TheoryClass(_set);
	_theory->setAllPars(x);

	TVectorD pred(13);

	for(int i=0;i<13;i++){
		  pred(i)=_theory->dGamma(bin(i));
	}

	TGraph *graph=new TGraph(bin,pred);

	graph->Draw("AC*");
	func->Draw("same");
	//func->Draw();


	for(int l=0;l<5;l++){
		DrawText(Form("+  #scale[1.0]{#it{ Parameter %i: %f}}",l,x[l]), kBlack,0.875-(l+2)*0.045,0.69);
	}
//	fun1->Draw("same");
	_can->Print(_ps);
	_can->Clear();
}

void PlotClass::PlotMatrix(TMatrixD* mat){
	TH2D *matplot=new TH2D(*mat);
	matplot->Draw("colztext45");
	_can->Print(_ps);
	_can->Clear();
}

void PlotClass::PlotMatrix(TMatrixDSym* mat){
	TH2D *matplot=new TH2D(*mat);
	matplot->Draw("colztext45");
	_can->Print(_ps);
	_can->Clear();
}

