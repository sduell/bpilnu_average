#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include "../utils/Utils.h"

using namespace std;

VecD VectorizeD(Str str, Str sep) {
  VecD result; StrV vecS = Vectorize(str,sep);
  for (uint i=0;i<vecS.size();++i) 
    result.push_back(atof(vecS[i]));
  return result;
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

int calc_sysCov(){

	double globfac=1e-06;
	ifstream infile;
	string namae;
	Str line;
	infile.open("fit_covmat.dat");
	cout << "read in file"<<endl;
	VecD covelements;
	int k=0;
	while(!infile.eof()){
		//cout << "read in test"<<endl;
		getline(infile,namae);
		line=Form("%s", namae.c_str());
		covelements=VectorizeD(line," ");
		//a[k]=VectorizeD(line,"	");
		//cout << line << endl;
		k++;
	}

	TMatrixDSym totalcov(sqrt(covelements.size()));
	k=0;
	for(int i =0; i<sqrt(covelements.size());i++){
		for(int j =0; j<sqrt(covelements.size());j++){
			totalcov(i,j)=covelements[k]*pow(globfac,2.);
			k++;
		}
	}

	ifstream statfile;
	statfile.open("stat_covmat.dat");
	cout << "Reading statistical covariance matrix..."<<endl;
	VecD statelements;
	k=0;
	while(!statfile.eof()){
		//cout << "read in test"<<endl;
		getline(statfile,namae);
		line=Form("%s", namae.c_str());
		statelements=VectorizeD(line," ");
	}
	TMatrixDSym statcov(sqrt(statelements.size()));
	k=0;
	for(int i =0; i<sqrt(statelements.size());i++){
		for(int j =0; j<sqrt(statelements.size());j++){
			statcov(i,j)=statelements[k]*pow(globfac,2.);
			k++;
		}
	}

	TMatrixDSym syscov=totalcov-statcov;
	
	totalcov.Print();
	statcov.Print();
	syscov.Print();

	ofstream offile;
	offile.open("sys_covmat.dat", ios::app);
	for(int i =0; i<sqrt(statelements.size());i++){
		for(int j =0; j<sqrt(statelements.size());j++){
			offile<< statcov(i,j)<<" ";
		}
	}

	infile.close();
	statfile.close();
	offile.close();

	VecD binwidth{2,2,2,2,2,2,2,2,2,2,2,2,2.4};

	double toterr=0., staterr=0., syserr=0.;
	for(int i =0; i<sqrt(statelements.size());i++){
		for(int j =0; j<sqrt(statelements.size());j++){
			toterr+=totalcov(i,j)*binwidth[i]*binwidth[j];
			staterr+=statcov(i,j)*binwidth[i]*binwidth[j];
			syserr+=syscov(i,j)*binwidth[i]*binwidth[j];
		}
	}

//	VecD vals{7.02962,7.12261,6.54852,7.37524,6.64364,7.17296,6.73647,6.26924,6.10849,4.40233,4.30536,3.38906,0.384453};
// Updated central values
//    VecD vals{7.27093, 7.10882, 6.74689, 7.57324, 6.4565, 7.20329, 6.6818, 6.35628, 6.21034, 4.33805, 4.25997, 3.40446, 1.18116};
    VecD vals{7.20444, 7.14098, 6.70338, 7.55861, 6.43738, 7.17247, 6.66887, 6.32951, 6.20395, 4.32124, 4.2493, 3.40174, 1.17164};
    double sum=0;
	double sumerr=0;
	double sumstaterr=0;
	double sumsyserr=0;

	for(int i=0; i<vals.size(); i++){
		sum+=vals[i]*binwidth[i]*globfac;
	}
	sumerr+=sqrt(toterr);
	sumstaterr+=sqrt(staterr);
	sumsyserr+=sqrt(syserr);
	//double fac=1e-02;
	double fac=1;
	cout << "value is "<<sum*fac<<" +- "<<sumstaterr*fac<<" (stat) +- "<<sumsyserr*fac << " (sys) "<< " ["<<sumerr*fac<<"]"<<endl;
//	cout << "total error is: " << sqrt(toterr)<<" Statistical error is: "<<sqrt(staterr)<<" Systematic error is: "<<sqrt(syserr)<<endl;
	

	return 0;
}
//this is form statistical only
//7.02992 pm 0.483668
//7.12233 pm 0.317555
//6.54837 pm 0.295979
//7.37508 pm 0.260857
//6.64319 pm 0.26725
//7.17254 pm 0.271499
//6.73628 pm 0.313261
//6.26856 pm 0.335593
//6.10826 pm 0.330218
//4.40207 pm 0.457463
//3.38912 pm 0.44776
//0.384678 pm 0.50651

//this is from the total fit
//7.02962 0.710376
//7.12261 0.412527
//6.54852 0.370186
//7.37524 0.387713
//6.64364 0.388042
//7.17296 0.415649
//6.73647 0.426105
//6.26924 0.445884
//6.10849 0.423992
//4.40233 0.515887
//4.30536 0.434968
//3.38906 0.48481
//0.384453        0.541021
