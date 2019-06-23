#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include "../utils/Utils.h"

using namespace std;

int main(){

    StrV measurements{"HaWon","Lees","Sanchez"};

    for(Str name : measurements){
        cout<< "Reading "<<name.Data()<<endl;
        ifstream statin;
        ifstream sysin;
        statin.open(Form("%s_stat.dat", name.Data()));
        string parse;
        Str line;
        VecVecD a;
        while(!statin.eof()){
            getline(statin,parse);
            line=Form("%s", parse.c_str());
            a.push_back(VectorizeD(line," "));
        }

        TMatrixDSym *StatMatrix = new TMatrixDSym(a.size());
        for(int j=0; j<a.size(); j++){
            for(int k=0; k<a[0].size(); k++){
                (*StatMatrix)(j,k)=a[j][k];
            }
        }
        statin.close();

        cout<<"Statistical correlation matrix read in"<<endl;
        StatMatrix->Print();

        sysin.open(Form("%s_sys.dat", name.Data()));
        VecVecD b;
        while(!sysin.eof()){
            getline(sysin,parse);
            line=Form("%s", parse.c_str());
            b.push_back(VectorizeD(line," "));
        }
        
        TMatrixDSym *SysMatrix = new TMatrixDSym(b.size());
        for(int j=0; j<b.size(); j++){
            for(int k=0; k<b[0].size(); k++){
                (*SysMatrix)(j,k)=b[j][k];
            }
        }
        sysin.close();
        
        cout<<"Systematic correlation matrix read in"<<endl;
        SysMatrix->Print();

        ifstream staterrin;
        ifstream syserrin;
        VecD staterr;
        VecD syserr;

        staterrin.open(Form("%s_staterr.dat", name.Data()));
        while(!staterrin.eof()){
            getline(staterrin,parse);
            line=Form("%s", parse.c_str());
            staterr=VectorizeD(line," ");
        }
        staterrin.close();

        syserrin.open(Form("%s_syserr.dat", name.Data()));
        while(!syserrin.eof()){
            getline(syserrin,parse);
            line=Form("%s", parse.c_str());
            syserr=VectorizeD(line," ");
        }
        syserrin.close();
        cout<<"calculating covariance matrices"<<endl;
        //Now calculate the covariance matrices
        TMatrixDSym *StatCovMatrix = new TMatrixDSym(a.size());
        for(int j=0; j<b.size(); j++){
            for(int k=0; k<b[0].size(); k++){
                (*StatCovMatrix)(j,k)=(*StatMatrix)(j,k)*staterr[j]*staterr[k];
            }
        }
        cout<< "Statistical covariance matrix constructed"<<endl;
        StatCovMatrix->Print();
        TMatrixDSym *SysCovMatrix = new TMatrixDSym(b.size());
        for(int j=0; j<b.size(); j++){
            for(int k=0; k<b[0].size(); k++){
                (*SysCovMatrix)(j,k)=(*SysMatrix)(j,k)*syserr[j]*syserr[k];
            }
        }
        cout<< "Systematic covariance matrix constructed"<<endl;
        SysCovMatrix->Print();
        //Now determine the total covariance and then correlation matrix
        TMatrixDSym *TotCovMatrix=new TMatrixDSym(a.size());
        //(*TotCovMatrix)= (*StatCovMatrix)+(*SysCovMatrix);
        TotCovMatrix->Plus((*StatCovMatrix),(*SysCovMatrix));
        cout<< "Total covariance matrix computed"<<endl;
        TMatrixDSym *TotMatrix = new TMatrixDSym((*TotCovMatrix));
        TMatrixDSym *DiagMatrix = new TMatrixDSym(a.size());
        
        for(int i=0;i<a.size(); i++){
            (*DiagMatrix)(i,i)=sqrt((*TotCovMatrix)(i,i));
        }
        TMatrixDSym *InversMatrix=new TMatrixDSym(b.size());
        (*InversMatrix)=(*DiagMatrix).Invert();
        TotMatrix->SimilarityT((*InversMatrix));
        cout<< "Total correlation matrix computed"<<endl;
        ofstream totout;
        totout.open(Form("%s_tot.dat", name.Data()));
        for(int j=0; j<a.size(); j++){
            for(int k=0; k<a[0].size(); k++){
                totout<<Form("%2.3f ", (*TotMatrix)(j,k));
            }
            totout<<endl;
        }
        totout.close();

        ofstream toterrout;
        toterrout.open(Form("%s_toterr.dat",name.Data()));
        for(int i=0; i<a.size();i++){
            toterrout<<Form("%2.3f ",sqrt((*TotCovMatrix)(i,i)));
        }
    }
    return 0;
}
