/*
 *  HiggsCombinerCommon: Florian Bernlochner
 */

#ifndef _HiggsCombinerCommon
#define _HiggsCombinerCommon
//UseTotCorMat was added for the total covariances!!!!
struct settings {
    Str PlotFileName, InputCard, Variable; StrV Measurements; bool UseTotCorMat ;bool IncludeNP; VecD CrossSectionOptions; Str Label, Unit; StrV MeasLabel;
    Str CombinationMethod, NdataModel;
    double ScanNPoints, ScanNSigma;
    Str ManualScan; bool DoPerformIndividualScan;
};




struct input {
    StrV Bkg_Uncertainties, BF_Vub_Uncertainties, BF_Vcb_Uncertainties, FF_Vub_Uncertainties, FF_Vcb_Uncertainties, general_Uncertainties;
    VecD Bins, Ndata, Nbkg, Nsig, NsigError, NsigErrorLow, NsigErrorHigh;
    VecVecD  Bkg_Uncert, BF_Vub_Uncert, BF_Vcb_Uncert, FF_Vub_Uncert, FF_Vcb_Uncert, general_Uncert;
    bool hasBkg, hasAsymUncert, hasLikelihoodContour; Str Name; double Luminosity, BR;
    Str ContourFileName; TGraphAsymmErrors *nll_contour; VecD NsigSM;
};


#endif
