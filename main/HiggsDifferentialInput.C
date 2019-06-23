/*
 *  HiggsCombiner: Florian Bernlochner
 */

#include "HiggsDifferentialInput.h"


HiggsInput::HiggsInput() { }

HiggsInput::HiggsInput(settings set) {
    
    Printf("> Initializing HiggsInput Class");
    
    _set = set;
    
    // Open input card and fill input vector
    var_input = OpenSettingsFile(_set.InputCard);
    FillInputVector();
    
}

void HiggsInput::FillInputVector() {
    
    Printf("\n  Found %lu measurements to be combined, reading in input", _set.Measurements.size());
    // Loop over all measurements
    for(int i = 0; i < _set.Measurements.size(); i++) {
        
        input in; Str var = _set.Measurements[i]+"."+_set.Variable;
        
        string nm=_set.Measurements[i].Data();

        cout << "\n\n \t" << _set.Measurements[i] << ":\n\n";
        
        // Read in Uncertainty sources for Fiducial acceptance, unfolding factor and Background
        
        in.Bkg_Uncertainties    = Vectorize(var_input->GetValue(var+".Bkg.Uncertainties","")," ");
        in.BF_Vub_Uncertainties = Vectorize(var_input->GetValue(var+".BF.Vub.Uncertainties","")," ");
        in.BF_Vcb_Uncertainties = Vectorize(var_input->GetValue(var+".BF.Vcb.Uncertainties","")," ");
        in.FF_Vub_Uncertainties = Vectorize(var_input->GetValue(var+".FF.Vub.Uncertainties","")," ");
        in.FF_Vcb_Uncertainties = Vectorize(var_input->GetValue(var+".FF.Vcb.Uncertainties","")," ");
        in.general_Uncertainties = Vectorize(var_input->GetValue(var+".general.Uncertainties","")," ");                
        // Read in the Bins, # of data / background / signal + error
        in.Bins = VectorizeD(var_input->GetValue(var+".Bins","")," ");
        in.Ndata = VectorizeD(var_input->GetValue(var+".Ndata","")," ");
        in.Nbkg = VectorizeD(var_input->GetValue(var+".Nbkg","")," ");
        in.Nsig = VectorizeD(var_input->GetValue(var+".Nsig","")," ");
        //rescale the measurement to newer tauB values!
        //if(nm=="Sanchez") {for(int oh=0; oh<in.Nsig.size();oh++){in.Nsig[oh]=in.Nsig[oh]*0.99535316;}}
        //else if(nm=="Sanchez") {for(int oh=0; oh<in.Nsig.size();oh++){in.Nsig[oh]=in.Nsig[oh]*1.002137546;}}
        in.NsigError = VectorizeD(var_input->GetValue(var+".NsigError","")," ");
        //if(nm=="Sanchez") {for(int oh=0; oh<in.NsigError.size();oh++){in.NsigError[oh]=in.NsigError[oh]*0.99535316;}}
        //if(nm=="Sanchez") {for(int oh=0; oh<in.NsigError.size();oh++){in.NsigError[oh]=in.NsigError[oh]*1.002137546;}}
        in.NsigErrorLow  = VectorizeD(var_input->GetValue(var+".NsigErrorLow","")," ");
        in.NsigErrorHigh = VectorizeD(var_input->GetValue(var+".NsigErrorHigh","")," ");
        in.ContourFileName = var_input->GetValue(var+".ContourFileName","");
        if(in.ContourFileName != "") {
            TFile *f = new TFile(in.ContourFileName);
            in.nll_contour = (TGraphAsymmErrors*) f->Get("Graph");
            in.hasLikelihoodContour = true;
            in.NsigSM = VectorizeD(var_input->GetValue(var+".NsigSM","")," ");
        } else in.hasLikelihoodContour = false;
        if(in.NsigErrorLow.size() != 0 && in.NsigErrorLow.size() == in.NsigErrorHigh.size()) in.hasAsymUncert = true;
        else in.hasAsymUncert = false;
        in.BR = var_input->GetValue(_set.Measurements[i]+".BR",1.);
        // Does this measurement have background?
        in.Nbkg.size() > 0 ? in.hasBkg = true : in.hasBkg = false;
        in.Luminosity = var_input->GetValue(var+".Luminosity",-99.);
        if(in.Luminosity < 0)
            in.Luminosity = var_input->GetValue("Luminosity",-99.);
        // Make sure we have either Nsig or Ndata, but not both
        if(in.Nsig.size() > 0 && in.Ndata.size() > 0)
            Fatal("HiggsInput: Check config file, Ndata and Nsig defined for a given measurement; only specify either","");
        // Read in the unfolding & acceptance factor
        
        // Loop over all Bkg uncertainty sources and read in the uncertainties
        Printf("\t Found %lu Bkg subtraction related error sources ", in.Bkg_Uncertainties.size());
        for(int i = 0; i < in.Bkg_Uncertainties.size(); i++) {
            Str s_uncert = in.Bkg_Uncertainties[i];
            VecD uncert = VectorizeD(var_input->GetValue(var+".Bkg."+s_uncert,"")," "); 
            if(uncert.size() == 0) 
                Fatal("Cannot find source: "+var+".Bkg."+s_uncert,"");
            //if(nm=="Sanchez") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*0.99535316;}}
            //else if(nm=="Lees") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*1.002137546;}} 
            in.Bkg_Uncert.push_back(uncert);
        }
        // Loop over all Vub BF factor uncertainty sources and read in the uncertainties
        Printf("\t Found %lu Vub BF Factor (BF_Vub) related error sources ", in.BF_Vub_Uncertainties.size());
        for(int i = 0; i < in.BF_Vub_Uncertainties.size(); i++) {
            Str s_uncert = in.BF_Vub_Uncertainties[i];
            VecD uncert = VectorizeD(var_input->GetValue(var+".BF.Vub."+s_uncert,"")," ");
            if(uncert.size() == 0)
                Fatal("Cannot find source: "+var+".BF.Vub."+s_uncert,"");
            //if(nm=="Sanchez") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*0.99535316;}}
            //else if(nm=="Lees") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*1.002137546;}}         
            in.BF_Vub_Uncert.push_back(uncert);
        }
        // Loop over all Vcb BF factor uncertainty sources and read in the uncertainties
        Printf("\t Found %lu Vcb BF Factor (BF_Vcb) related error sources ", in.BF_Vcb_Uncertainties.size());
        for(int i = 0; i < in.BF_Vcb_Uncertainties.size(); i++) {
            Str s_uncert = in.BF_Vcb_Uncertainties[i];
            VecD uncert = VectorizeD(var_input->GetValue(var+".BF.Vcb."+s_uncert,"")," ");
            if(uncert.size() == 0)
                Fatal("Cannot find source: "+var+".BF.Vcb."+s_uncert,"");
            //if(nm=="Sanchez") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*0.99535316;}}
            //else if(nm=="Lees") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*1.002137546;}}
            in.BF_Vcb_Uncert.push_back(uncert);
        }
        // Loop over all Vub FF uncertainty sources and read in the uncertainties
        Printf("\t Found %lu Vub Form Factor (FF_Vub) related error sources ", in.FF_Vub_Uncertainties.size());
        for(int i = 0; i < in.FF_Vub_Uncertainties.size(); i++) {
            Str s_uncert = in.FF_Vub_Uncertainties[i];
            VecD uncert = VectorizeD(var_input->GetValue(var+".FF.Vub."+s_uncert,"")," ");
            if(uncert.size() == 0)
                Fatal("Cannot find source: "+var+".FF.Vub."+s_uncert,"");
            //if(nm=="Sanchez") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*0.99535316;}}
            //else if(nm=="Lees") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*1.002137546;}} 
            in.FF_Vub_Uncert.push_back(uncert);
        }
        // Loop over all Vcb BF factor uncertainty sources and read in the uncertainties
        Printf("\t Found %lu Vcb Form Factor (FF_Vcb) related error sources ", in.FF_Vcb_Uncertainties.size());
        for(int i = 0; i < in.FF_Vcb_Uncertainties.size(); i++) {
            Str s_uncert = in.FF_Vcb_Uncertainties[i];
            VecD uncert = VectorizeD(var_input->GetValue(var+".FF.Vcb."+s_uncert,"")," ");
            if(uncert.size() == 0)
                Fatal("Cannot find source: "+var+".FF.Vcb."+s_uncert,"");
            //if(nm=="Sanchez") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*0.99535316;}}
            //else if(nm=="Lees") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*1.002137546;}}
            in.FF_Vcb_Uncert.push_back(uncert);
        }
        // Loop over all general factor uncertainty sources and read in the uncertainties
        Printf("\t Found %lu general Factor (general) related error sources ", in.general_Uncertainties.size());
        for(int i = 0; i < in.general_Uncertainties.size(); i++) {
            Str s_uncert = in.general_Uncertainties[i];
            VecD uncert = VectorizeD(var_input->GetValue(var+".general."+s_uncert,"")," ");
            if(uncert.size() == 0)
                Fatal("Cannot find source: "+var+".FF.Vcb."+s_uncert,"");
            //if(nm=="Sanchez") {cout << "inside!"<<endl;for(int oh=0; oh<uncert.size();oh++){cout << "round " <<oh<< " with "<< uncert[oh]<<endl;uncert[oh]=uncert[oh]*0.99535316;}}
            //else if(nm=="Lees") {for(int oh=0; oh<uncert.size();oh++){uncert[oh]=uncert[oh]*1.002137546;}}
            in.general_Uncert.push_back(uncert);
        }
        // Save name
        in.Name = _set.Measurements[i] ; 
        // Store the measurement for later use
        _v_input.push_back( in );   
    }
    
    Printf("\n");
    
    
}
