#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include "../utils/Utils.h"

using namespace std;

int main(){

	ifstream infile;
	string namae;
	Str line;
	infile.open("errorsources");
	cout << "read in file"<<endl;
	vector<VecD> a;
	int k=0;
	while(!infile.eof()){
		//cout << "read in test"<<endl;
		getline(infile,namae);
		line=Form("%s", namae.c_str());
		a.push_back(VecD());
		a[k]=VectorizeD(line,"	");
		//cout << line << endl;
		k++;
	}
	//Calculate all Covariance matrices

	TMatrixDSym *TotMatrix = new TMatrixDSym(a[0].size());

	for(int i=0; i<a.size(); i++){
		TMatrixDSym *SyMatrix = new TMatrixDSym(a[i].size());
		for(int j=0; j<a[i].size(); j++){
			for(int k=0; k<a[i].size(); k++){
				(*SyMatrix)(j,k)=a[i][j]*a[i][k]/10000.;
				(*TotMatrix)(j,k)=(*TotMatrix)(j,k)+(*SyMatrix)(j,k);
				//TotMatrix->Plus(*TotMatrix,*SyMatrix);
			}
		}
	}

	cout << "Total Matrix is "<<endl;
	//TotMatrix->Print();

	TMatrixDSym *TestMatrix = new TMatrixDSym(a[0].size());

	double array[]={7.8,4.8,6.0,6.8,7.1,7.6};
	for(int j=0; j<6; j++){
		for(int k=0; k<6; k++){
			(*TestMatrix)(j,k)=array[j]*array[k]/10000.;
		}
	}
	cout << "Matrix from Total is "<<endl;
	//TestMatrix->Print();

	TMatrixDSym *DiffMatrix = new TMatrixDSym(a[0].size());
	DiffMatrix->Minus((*TestMatrix),(*TotMatrix));
	//DiffMatrix->Print();


	//Now calculate the correlation matrix from the covariance matrix!
	TMatrixDSym *DiagMatrix = new TMatrixDSym(a[0].size());
	TMatrixDSym *InvertMatrix = new TMatrixDSym(a[0].size());
	for(int i=0; i<6; i++){
		(*DiagMatrix)(i,i)=sqrt((*TotMatrix)(i,i));
		//(*DiagMatrix)(i,i)=sqrt((*DiffMatrix)(i,i));
	}
	//DiagMatrix->Print();
	(*InvertMatrix)=(*DiagMatrix).Invert();
	//TMatrixDSym *SysCorMatrix = new TMatrixDSym(a[0].size());
	TMatrixDSym *SysCorMatrix = new TMatrixDSym(*TotMatrix);
	//SysCorMatrix=DiffMatrix;
	//SysCorMatrix=TotMatrix;
	SysCorMatrix->SimilarityT((*InvertMatrix));
//	SysCorMatrix=InvertMatrix;
//	(*SysCorMatrix)*=(*TotMatrix);
//	(*SysCorMatrix)*=(*InvertMatrix);
	//SysCorMatrix->TMult(*TotMatrix);
	//SysCorMatrix->Print();
	//SysCorMatrix->TMult(*InvertMatrix);
	cout << "Systematic Correlation Matrix is "<<endl;
	//SysCorMatrix->Print();
//	cout << "Read in done"<<endl;
//	for(int i=0; i<a.size(); i++){
//		//cout << "vector size "<<a.size()<<endl;
//		for(int j=0; j<a[i].size(); j++){
//			cout << "output is: "<<endl;
//			cout << a[i][j]<<endl;
//		}
//	}

	//Now read in the statistical correlation matrix and the standard deviations
	
	ifstream full;
	full.open("Sanchez_full.dat");
	vector<VecD> b;
	int n=0;
	while(!full.eof()){
		//cout << "read in test"<<endl;
		getline(full,namae);
		line=Form("%s", namae.c_str());
		b.push_back(VecD());
		b[n]=VectorizeD(line,"	");
		//cout << line << endl;
		n++;
	}
	cout << "matrix size is "<<b.size()<<endl;
	TMatrixDSym *FullCorMatrix = new TMatrixDSym(b.size());
	for(int i=0; i<b.size(); i++){
		for(int j=0; j<b[i].size(); j++){
			(*FullCorMatrix)(i,j)=b[i][j];
		}
	}
	cout << "FullCorrelationmatrix is: "<<endl;
	//FullCorMatrix->Print();

	//Now calculate the covariance matrix from the correlation matrix
	vector<VecD> c;
	n=0;

	ifstream stddev;
	stddev.open("Sanchez_totalerrors_1mode");

	while(!stddev.eof()){
		//cout << "read in test"<<endl;
		getline(stddev,namae);
		line=Form("%s", namae.c_str());
		c.push_back(VecD());
		c[n]=VectorizeD(line,"	");
		//cout << line << endl;
		n++;
	}

	TMatrixDSym *ErrorDiagMatrix = new TMatrixDSym(b.size());

	for(int i=0; i<b.size(); i++){
		(*ErrorDiagMatrix)(i,i)=c[0][i];
	}

	TMatrixDSym *FullCovMatrix = new TMatrixDSym(*FullCorMatrix);

	//FullCovMatrix=FullCorMatrix;
	FullCovMatrix->SimilarityT((*ErrorDiagMatrix));
	cout << "FullCovMatrix is: "<<endl;
	FullCovMatrix->Print();
	TotMatrix->Print();
	TMatrixDSym *StatCovMatrix = new TMatrixDSym(b.size());

	StatCovMatrix->Minus((*FullCovMatrix), (*TotMatrix));
	cout << "Statistical Covariance Matrix is "<<endl;
	StatCovMatrix->Print();

	TMatrixDSym *DiagMat = new TMatrixDSym(a[0].size());
	for(int i=0; i<6; i++){
		(*DiagMat)(i,i)=sqrt((*StatCovMatrix)(i,i));
	}
	TMatrixDSym *InversMatrix = new TMatrixDSym(a[0].size());
	(*InversMatrix)=(*DiagMat).Invert();

	TMatrixDSym *StatCorMatrix = new TMatrixDSym(*StatCovMatrix);
	StatCorMatrix->SimilarityT((*InversMatrix));
	cout << "Statistical correlation matrix is: "<<endl;
	StatCorMatrix->Print();

	
	return 0;
}

