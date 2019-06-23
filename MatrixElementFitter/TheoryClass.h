/**************************************************************************
 * Theory Class                                                           *
 *                                                                        *
 * Authors: F. Bernlochner, S.Duell                                       *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

//#include <analysis/ReweightingTool/FFTheoryClass.h>
#include "../utils/Utils.h"
#include "TMatrix.h"

class TheoryClass{
 public: 
 TheoryClass(TEnv *set);
// ~TheoryClass();

  double BF(double ml, double qmin, double qmax);
  double dBF(double q2);
  double Gamma(double ml, double qmin, double qmax);
  double dGamma(double q2);

 public:
  double pX(double q2);
  double f0(double q2);
  double fp(double q2);
  double z(double q2);

 public:
   //void setAllPars(double *x){cons.MatrixElement=x[0]; for(int i=1;i<5;i++){cons.bpPars[i-1]=x[i];}; for(int i=5;i<9;i++){cons.b0Pars[i-5]=x[i];};}; 
   void setFFpars(TVectorD pars){for(int i=0;i<cons.bpPars.size();i++){cons.bpPars[i]=pars(i);};};
   void setAllPars(double *x){cons.MatrixElement=x[0]; for(int i=1;i<4;i++){cons.bpPars[i-1]=x[i];};}; 
   void setfp(double *x){for(int i=0;i<cons.bpPars.size();i++){cons.bpPars[i]=x[i];};};
   void setf0(double *x){for(int i=0;i<cons.b0Pars.size();i++){cons.b0Pars[i]=x[i];};};
   void setMatrixElement(double* x){cons.MatrixElement=x[0];};
 private: 
 TEnv *m_set;
 string m_dec;
 
struct constants {
  double mB, mBs, mX, tauB, GammaB, GF, GF2, hbar, mtau, mb, mq, q2max;
  std::vector<double> bpPars, b0Pars, fpPars, f0Pars;
  double tp, t0;
  double MatrixElement; } cons;      

};
