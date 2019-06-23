/*
 *  DifferentialCombiner: Florian Bernlochner
 */

#include "../utils/Utils.h"

#include <TCanvas.h>
#include "HiggsDifferentialCommon.h"
#include "HiggsDifferentialCombiner.h"

// ---------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    
    
    Printf("\n\n> DifferentialCombiner -- combine Higgs to Boson differential information \n\n");
    
    settings set; SetAtlasStyle(); gStyle->SetPalette(1); gStyle->SetHistMinimumZero();
    
    // Reads in the settings, locations of files, etc
    Str SettingsFile = argc == 1 ? "combine.config" : argv[1];
    TEnv *global_settings = OpenSettingsFile(SettingsFile);
    // -------------------------------------------------------------------------------------
    // Default global settings
    set.PlotFileName     = global_settings->GetValue("PlotFileName","combination.ps");
    set.InputCard        = global_settings->GetValue("InputCard","");
    set.Variable         = global_settings->GetValue("Variable","");
    set.Measurements     = Vectorize(global_settings->GetValue("Measurements","")," ");
    set.IncludeNP        = global_settings->GetValue("IncludeNP",true);
    //This was added for the total covariances!!!
    set.UseTotCorMat     = global_settings->GetValue("UseTotCorMat",false);
    set.CrossSectionOptions = VectorizeD(global_settings->GetValue("CrossSectionOptions","")," ");
    if(set.CrossSectionOptions.size() != 3) Fatal("Specify CrossSectionOptions, exiting","");
    set.Label             = global_settings->GetValue(set.Variable+".Label","");
    set.Unit              = global_settings->GetValue(set.Variable+".Unit","NONE");
    set.CombinationMethod = global_settings->GetValue("CombinationMethod","full");
    set.NdataModel        = global_settings->GetValue("NdataModel","Gaussian");
    for(int i = 0; i < set.Measurements.size(); i++)
        set.MeasLabel.push_back( global_settings->GetValue(set.Measurements[i]+".Label",set.Measurements[i]) );
    // Manual scan; number of sigma and points to sample
    set.ManualScan = global_settings->GetValue("ManualScan","false");
    set.ScanNSigma = global_settings->GetValue("ScanNSigma",4.);
    set.ScanNPoints = global_settings->GetValue("ScanNPoints",40.);
    // Create a unique output file
    set.PlotFileName.ReplaceAll("VAR",set.Variable);
    set.PlotFileName.ReplaceAll("METHOD",set.CombinationMethod);

    set.DoPerformIndividualScan = global_settings->GetValue("DoPerformIndividualScan",false);
    
    // -------------------------------------------------------------------------------------
    // Create a combiner object
    HiggsCombiner comb(set);
    // perform the actual combination
    comb.Combine();
    
    // -------------------------------------------------------------------------------------
    
    
}


