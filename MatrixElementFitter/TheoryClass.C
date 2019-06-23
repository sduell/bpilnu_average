/**************************************************************************
 * HAMMER Tool for Form Factor Reweighting                                *
 *                                                                        *
 * Authors: F. Bernlochner, S.Duell, Z.Ligeti, M. Papucci, D.Robinson     *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

//#include <analysis/ReweightingTool/include/TheoryClass.h>
#include "TheoryClass.h"

TheoryClass::TheoryClass(TEnv *set)
{
m_set = set;

string matel="Vub"; string m_mq="u";

cons.mB   = m_set->GetValue("mass.B",1.);
cons.mBs  = m_set->GetValue("mass.Bs",1.);
cons.mX  = m_set->GetValue("mass.X",1.);
cons.tauB = m_set->GetValue("tauB",1.);
cons.hbar = m_set->GetValue("hbar",1.);
cons.GF   = m_set->GetValue("GF",1.);
cons.GF2  = pow(cons.GF,2.);
cons.mtau = m_set->GetValue("mtau",1.);
cons.mb   = m_set->GetValue("mass.b",1.);
cons.mq   = m_set->GetValue(Form("mass.%s",m_mq.c_str()),1.);
cons.MatrixElement   = m_set->GetValue(Form("MatrixElement.%s", matel.c_str()),1.);
cons.q2max = pow(cons.mB,2.) + pow(cons.mX,2.) - 2 * cons.mB*cons.mX;
    
cons.tp   = pow(cons.mB + cons.mX, 2.);
cons.t0   = (cons.mB+cons.mX)*pow(sqrt(cons.mB)-sqrt(cons.mX),2.);


//cons.bpPars = VectorizeD(m_set->GetValue("bpPars","")," ");
cons.b0Pars = VectorizeD(m_set->GetValue("b0Pars","")," ");
cons.bpPars = VectorizeD(m_set->GetValue("apPars","")," ");

cons.fpPars = VectorizeD(m_set->GetValue("fpPars","")," ");
cons.f0Pars = VectorizeD(m_set->GetValue("f0Pars","")," ");
}

double TheoryClass::BF(double ml, double qmin, double qmax) {
  double value =0.;
  value = Gamma(ml,qmin,qmax);
  value*=cons.tauB;
  value/=cons.hbar;
  //cout << "FF Parameters: bp0 "<< cons.bpPars[0]<<" bp1 "<<cons.bpPars[1]<<" bp2 "<<cons.bpPars[2]<<" bp3 "<<cons.bpPars[3]<<endl;
  //cout << "The constants are: GF2 "<<cons.GF2<<" MatrixElement "<<cons.MatrixElement<<" mB "<<cons.mB<<" cons.mBs "<<cons.mBs<<" tp "<<cons.tp<<" t0 "<<cons.t0<<" mX "<<cons.mX<<" GammaB "<<cons.hbar/cons.tauB<<endl;
  //cout <<"Rate value from "<<qmin <<" to "<<qmax<<" is: "<< Gamma(ml,qmin,qmax) <<" Branching Fraction Value is: "<<value<<endl;
  //cout <<"fp is "<<fp((qmin+qmax)/2.)<<" z is "<<z((qmin+qmax)/2.)<<" Ppi is "<< pX((qmin+qmax)/2.) <<" and dGamma is (on point) "<< dGamma((qmin+qmax)/2.)<<endl;
  //cout <<"f0 is "<<f0((qmin+qmax)/2.)<<endl;
  return value;
}

double TheoryClass::dBF(double q2) {
  double value =0.;
  value = dGamma(q2);
  value*=cons.tauB;
  value/=cons.hbar;
  return value;
}

double TheoryClass::Gamma(double ml, double qmin, double qmax) {

     double value = 0.;  enum { n = 19 };

     double xi[n]; double wi[n];

     xi[0] = 0.0;                  wi[0] = 0.16105444984878362;
     xi[1] = 0.16035862813102394;  wi[1] = 0.15896882598519965;
     xi[2] = -xi[1];               wi[2] = wi[1];
     xi[3] = 0.3165642266117781;   wi[3] = 0.1527663007284602;
     xi[4] = -xi[3];               wi[4] = wi[3];
     xi[5] = 0.46457038687717755;  wi[5] = 0.1426055640976579;
     xi[6] = -xi[5];               wi[6] = wi[5];
     xi[7] = 0.6005459122435153;   wi[7] = 0.1287567549009922;
     xi[8] = -xi[7];               wi[8] = wi[7];
     xi[9] = 0.7209655043062628;   wi[9] = 0.11156236184219351;
     xi[10] = -xi[9];              wi[10] = wi[9];
     xi[11] = 0.8227150489105877;  wi[11] = 0.09149349482553046;
     xi[12] = -xi[11];             wi[12] = wi[11];
     xi[13] = 0.9031559802279142;  wi[13] = 0.06904552773838758;
     xi[14] = -xi[13];             wi[14] = wi[13];
     xi[15] = 0.9602078316107117;  wi[15] = 0.044807508218039215;
     xi[16] = -xi[15];             wi[16] = wi[15];
     xi[17] = 0.9924070086164836;  wi[17] = 0.019469784466159833;
     xi[18] = -xi[17];             wi[18] = wi[17];

     for(int i=0; i < n; i++) {
      double yi  = (qmax-qmin)*xi[i]/2.0 + (qmax+qmin)/2.0;
      //cout << "The variables for q2="<<yi<< " are: pX "<<pX(yi)<<" fp "<<fp(yi)<<" z "<<z(yi)<<endl;
      double dGi = dGamma(yi);
      value +=  wi[i]*dGi;
     }

     value *= (qmax-qmin)/2.0;
     //cout << "FF Parameters: bp0 "<< cons.bpPars[0]<<" bp1 "<<cons.bpPars[1]<<" bp2 "<<cons.bpPars[2]<<" bp3 "<<cons.bpPars[3]<<endl;
     //cout << "The constants are: GF2 "<<cons.GF2<<" MatrixElement "<<cons.MatrixElement<<" mB "<<cons.mB<<" cons.mBs "<<cons.mBs<<" tp "<<cons.tp<<" t0 "<<cons.t0<<" mX "<<cons.mX<<endl;

     return value;

}

double TheoryClass::dGamma(double q2){
	double rate=0;
  //cout << "debugging dgamma : "<<q2<<" "<<cons.GF2<<" "<<cons.MatrixElement<<" "<<cons.mB<<" "<<cons.mBs<<" "<<cons.mX<<" "<<pX(q2)<<" "<<fp(q2)<<" "<<endl;
  /*rate=pow(fp(q2),2.);
  rate*=sqrt(pow(pow(pow(cons.mB,2.)+pow(cons.mX,2.)-q2,2.)-4.*pow(cons.mB,2.)*pow(cons.mX,2.),3.));
	rate*=cons.GF2*pow(cons.MatrixElement,2.);
	rate/=192*pow(TMath::Pi(),3.)*pow(cons.mB,3.);*/
  //rate=cons.GF2*pow(cons.MatrixElement,2.)*pow(pX(q2),3.)*pow(fp(q2),2.)/(24.*pow(TMath::Pi(),3.));
	rate=8*pX(q2)/3.;
  rate*=cons.GF2*pow(cons.MatrixElement,2.)*q2/(256*pow(TMath::Pi(),3.)*pow(cons.mB,2.));
  rate*=pow(2.*cons.mB*pX(q2)*fp(q2)/sqrt(q2),2.);
  return rate;
}

double TheoryClass::pX(double q2)
{
  return 1./(2.*cons.mB)*sqrt(pow(pow(cons.mB,2.)+pow(cons.mX,2.)-q2,2.)-4.*pow(cons.mB,2.)*pow(cons.mX,2.));  
}

double TheoryClass::f0(double q2)
{
 double value = 0.;
 for(int ipar = 0; ipar < cons.b0Pars.size(); ipar++)
 value += cons.b0Pars[ipar] * ( pow(z(q2),ipar) );
 return value;
}

double TheoryClass::fp(double q2)
{
  double value = 0.;
  for(int ipar = 0; ipar < cons.bpPars.size(); ipar++)
  {
   //value += cons.bpPars[ipar] * ( pow(z(q2),ipar) - (pow(-1,ipar-cons.bpPars.size())*(ipar)*pow(z(q2),cons.bpPars.size())/cons.bpPars.size()));
   //value += cons.bpPars[ipar] * ( pow(z(q2),ipar) - (pow(-1,ipar-4.)*(ipar)*pow(z(q2),4.)/(4.)));
    value += cons.bpPars[ipar] * ( pow(z(q2),ipar) - (pow(-1,ipar-3.)*(ipar)*pow(z(q2),3.)/(3.)));
   //cout << "bp parameter no. "<<ipar<< " is "<<cons.bpPars[ipar]<<endl;
   //value += cons.bpPars[ipar] * ( pow(z(q2),ipar) - pow(-1,ipar)*(ipar)/4*pow(z(q2),cons.bpPars.size()) );
  }
  value *=  1/(1-(q2/pow(cons.mBs,2.)));
  return value;
}

double TheoryClass::z(double q2)
{
 return (sqrt(cons.tp-q2)-sqrt(cons.tp-cons.t0))/(sqrt(cons.tp-q2)+sqrt(cons.tp-cons.t0)); 
}
