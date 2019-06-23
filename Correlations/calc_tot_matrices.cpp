#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include "../utils/Utils.h"

using namespace std;

int main(){

    StrV measurements{"HaWon","Lees","Sanchez"};

    for(Str name : measurements){
        infile statin;
        infile sysin;
        statin.open(Form("%s_stat.dat", name.Data()));
        string parse;
        Str line;
        VecVecD a;
        while(!statin.eof()){
            getline(statin,parse);
            line=Form("%s", parse.c_str());
            a.push_back(VectorizeD(line," "));
        }
cout<<"test"<<endl;
        TMatrixDSym *StatMatrix = new TMatrixDSym(a.size());
        for(int j=0; j<a.size(); j++){
            for(int k=0; k<a[0].size(); k++){
                (*StatMatrix)(j,k)=a[j][k];
            }
        }
        statin.close();
cout<<"test"<<endl;
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
        

        infile staterrin;
        infile syserrin;
        VecD staterr;
        VecD syserr;

        staterrin.open(Form("%s_staterr.dat", name.Data()));
        while(!staterrin.eof()){
            getline(staterrin,parse);
            line=Form("%s", parse.c_str());
            staterr.push_back(VectorizeD(line," "));
        }
        staterrin.close();
cout<<"test"<<endl;cout<<"test"<<endl;
        syserrin.open(Form("%s_syserr.dat", name.Data()));
        while(!syserrin.eof()){
            getline(syserrin,parse);
            line=Form("%s", parse.c_str());
            syserr.push_back(VectorizeD(line," "));
        }
        syserrin.close();

        //Now calculate the covariance matrices
        TMatrixDSym *StatCovMatrix = new TMatrixDSym(a.size());
        for(int j=0; j<b.size(); j++){
            for(int k=0; k<b[0].size(); k++){
                (*StatCovMatrix)(j,k)=(*StatMatrix)(j,k)*staterr[j]*staterr[k];
            }
        }

        TMatrixDSym *SysCovMatrix = new TMatrixDSym(b.size());
        for(int j=0; j<b.size(); j++){
            for(int k=0; k<b[0].size(); k++){
                (*SysCovMatrix)(j,k)=(*SysMatrix)(j,k)*syserr[j]*syserr[k];
            }
        }
cout<<"test"<<endl;
        //Now determine the total covariance and then correlation matrix
        TMatrixDSym *TotCovMatrix = (*StatCovMatrix)+(*SysCovMatrix);
        TMatrixDSym *TotMatrix = new TMatrixDSym((*TotCovMatrix));
        TMatrixDSym *DiagMatrix = new TMatrixDSym(a.size());

        for(int i=0;i<a.size(); i++){
            (*DiagMatrix)(i,i)=sqrt((*TotCovMatrix)(i,i));
        }
        TMatrixDSym *InversMatrix=(*DiagMatrix).Invert();
        TotMatrix->SimilarityT((*InversMatrix));

        outfile totout;
        totout.open(Form("%s_tot.dat", name.Data()));
        for(int j=0; j<a.size(); j++){
            for(int k=0; k<a[0].size(); k++){
                Printf("%2.2f", (*TotMatrix)(j,k));
            }
            Printf("\n");
        }
        totout.close();
        cout<<"test"<<endl;cout<<"test"<<endl;
        outfile toterrout;
        toterrout.open(Form("%s_toterr.dat",name.Data()));
        for(int i=0; i<a.size();i++){
            Printf("%2.2f",sqrt((*TotCovMatrix)(i,i)));
        }
    }
    return 0;
}