/*
 *  HiggsCombiner: Florian Bernlochner
 */

#include "HiggsDifferentialWorkspace.h"


HiggsWorkspace::HiggsWorkspace() { }

HiggsWorkspace::HiggsWorkspace(settings set, HiggsInput *input) {
    
    Printf("> Initializing HiggsWorkspace Class\n");
    
    _set = set; _input = input; _workspace = new RooWorkspace("DifferentialCombination");
    
    _SyscovMatrix = _StatcovMatrix = _TotcovMatrix = 0; _n_meas = 0; _finest_index = 0; _max_bins = 0;
    
    // Determine number of measurements
    for(int i = 0; i < _input->GetInput().size(); i++)
        _n_meas += _input->GetInput()[i].Bins.size()-1;
    
    // Determine offset between covariance elements
    _start_bin.push_back(0);
    for(int i = 0; i < _input->GetInput().size(); i++){
        _start_bin.push_back(_input->GetInput()[i].Bins.size()-1+_start_bin[i]);
        cout<<"Testing starting bin "<<_start_bin[i]<<endl;
    }
    
    // Construct the total covariance matrix (used for plotting & for the simplified workspace)
    ConstructCovariance();
    
    // Setup the desired workspace
    if(_set.CombinationMethod == "full" || _set.CombinationMethod == "shape")
        SetupWorkspace();
    else if(_set.CombinationMethod == "Multi")
        SetupMultiWorkspace();
    else
        SetupSimplifiedWorkspace();

    //Fatal("","");
            
}

void HiggsWorkspace::SetupMultiWorkspace() {
    
    // Determine the highest granularity and add bin cross sections into the workspace as free parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        if( _max_bins < in.Bins.size() ) {
            _max_bins = in.Bins.size(); _finest_index = i;
        }
    }
    for(int i = 0; i < _max_bins; i++) {
        _workspace->factory( Form("sigma_%d[%f,%f,%f]",i+1,_set.CrossSectionOptions[0],_set.CrossSectionOptions[1],_set.CrossSectionOptions[2]) );
    }
    
//Now read in Nuisance Parameters
        // Read in measurements, unfolding & acceptance factors, background yields, Nuisance parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        // Branching Fraction
        _workspace->factory( Form("BR_%s[%f]",in.Name.Data(),in.BR));
        // Luminosity
        _workspace->factory( Form("Luminosity_%s[%f]",in.Name.Data(),in.Luminosity));
        
        // Bkg yields
        for(int j = 0; j < in.Nbkg.size(); j++) {
            _workspace->factory( Form("Nbkg_%s_%d[%f]",in.Name.Data(),j+1,in.Nbkg[j]));
        }




        // Include the three classes of Nuisance parameters and constraints (CFactors, Acceptance, Nbkg) into workspace if required
        if(_set.IncludeNP) {
            cout <<"IncludeNP in Workspace activated"<<endl;
            // ------------------------------------------------------------------------------------
            // Constraints & NP
            
            for(int j = 0; j < in.Bkg_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.Bkg_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.Bkg_Uncertainties[j].Data(),in.Bkg_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.Bkg_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.BF_Vub_Uncertainties.size(); j++) {
                // Check if constrained already exists
                //cout << "object tbc is: "<<Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vub_Uncertainties[j].Data(),in.BF_Vub_Uncertainties[j].Data())<<endl;
                if( !_workspace->obj( Form("Uncert_%s_const", in.BF_Vub_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vub_Uncertainties[j].Data(),in.BF_Vub_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.BF_Vub_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.BF_Vcb_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.BF_Vcb_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vcb_Uncertainties[j].Data(),in.BF_Vcb_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.BF_Vcb_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.FF_Vub_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.FF_Vub_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.FF_Vub_Uncertainties[j].Data(),in.FF_Vub_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.FF_Vub_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.FF_Vcb_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.FF_Vcb_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.FF_Vcb_Uncertainties[j].Data(),in.FF_Vcb_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.FF_Vcb_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.general_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.general_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.general_Uncertainties[j].Data(),in.general_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.general_Uncertainties[j].Data()) );
                }
            }
            // ------------------------------------------------------------------------------------
            // Errors
            for(int j = 0; j < in.Bkg_Uncertainties.size(); j++)
                for(int k = 0; k < in.Bkg_Uncert[j].size(); k++) {
                    Str s_Bkg_Uncert_Value = Form("Bkg_Uncert_%s_%s_Value_%d", in.Name.Data(), in.Bkg_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_Bkg_Uncert_Value.Data(), in.Bkg_Uncert[j][k]) );
                    _workspace->factory( Form("expr::Bkg_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.Bkg_Uncertainties[j].Data(), k+1, in.Bkg_Uncertainties[j].Data(),s_Bkg_Uncert_Value.Data(),in.Bkg_Uncertainties[j].Data(), s_Bkg_Uncert_Value.Data()) );
                }
            for(int j = 0; j < in.BF_Vub_Uncertainties.size(); j++)
                for(int k = 0; k < in.BF_Vub_Uncert[j].size(); k++) {
                    //cout << "error tbc is: "<<Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vub_Uncertainties[j].Data(),in.BF_Vub_Uncertainties[j].Data())<<endl;
                    Str s_BF_Vub_Uncert_Value = Form("BF_Vub_Uncert_%s_%s_Value_%d", in.Name.Data(), in.BF_Vub_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_BF_Vub_Uncert_Value.Data(), in.BF_Vub_Uncert[j][k]) );
                    _workspace->factory( Form("expr::BF_Vub_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.BF_Vub_Uncertainties[j].Data(), k+1, in.BF_Vub_Uncertainties[j].Data(),s_BF_Vub_Uncert_Value.Data(),in.BF_Vub_Uncertainties[j].Data(), s_BF_Vub_Uncert_Value.Data()) ); 
                }
            for(int j = 0; j < in.BF_Vcb_Uncertainties.size(); j++)
                for(int k = 0; k < in.BF_Vcb_Uncert[j].size(); k++) {
                    Str s_BF_Vcb_Uncert_Value = Form("BF_Vcb_Uncert_%s_%s_Value_%d", in.Name.Data(), in.BF_Vcb_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_BF_Vcb_Uncert_Value.Data(), in.BF_Vcb_Uncert[j][k]) );
                    _workspace->factory( Form("expr::BF_Vcb_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.BF_Vcb_Uncertainties[j].Data(), k+1, in.BF_Vcb_Uncertainties[j].Data(),s_BF_Vcb_Uncert_Value.Data(),in.BF_Vcb_Uncertainties[j].Data(), s_BF_Vcb_Uncert_Value.Data()) );
                }
                //cout << "After BF Uncertainties!"<<endl;
                //cout << "Next is Vub FF Uncertainty with size "<< in.FF_Vub_Uncertainties.size() << endl;
            for(int j = 0; j < in.FF_Vub_Uncertainties.size(); j++){
                for(int k = 0; k < in.FF_Vub_Uncert[j].size(); k++) {
                    //cout << "FF error tbc is: "<<Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.FF_Vub_Uncertainties[j].Data(),in.FF_Vub_Uncertainties[j].Data())<<endl;
                    Str s_FF_Vub_Uncert_Value = Form("FF_Vub_Uncert_%s_%s_Value_%d", in.Name.Data(), in.FF_Vub_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_FF_Vub_Uncert_Value.Data(), in.FF_Vub_Uncert[j][k]) );
                    _workspace->factory( Form("expr::FF_Vub_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.FF_Vub_Uncertainties[j].Data(), k+1, in.FF_Vub_Uncertainties[j].Data(),s_FF_Vub_Uncert_Value.Data(),in.FF_Vub_Uncertainties[j].Data(), s_FF_Vub_Uncert_Value.Data()) );
                }}
                //cout <<"test to see if NP correctly initiated "<<endl;
            for(int j = 0; j < in.FF_Vcb_Uncertainties.size(); j++)
                for(int k = 0; k < in.FF_Vcb_Uncert[j].size(); k++) {
                    Str s_FF_Vcb_Uncert_Value = Form("FF_Vcb_Uncert_%s_%s_Value_%d", in.Name.Data(), in.FF_Vcb_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_FF_Vcb_Uncert_Value.Data(), in.FF_Vcb_Uncert[j][k]) );
                    _workspace->factory( Form("expr::FF_Vcb_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.FF_Vcb_Uncertainties[j].Data(), k+1, in.FF_Vcb_Uncertainties[j].Data(),s_FF_Vcb_Uncert_Value.Data(),in.FF_Vcb_Uncertainties[j].Data(), s_FF_Vcb_Uncert_Value.Data()) );
                }

            for(int j = 0; j < in.general_Uncertainties.size(); j++)
                for(int k = 0; k < in.general_Uncert[j].size(); k++) {
                    Str s_general_Uncert_Value = Form("general_Uncert_%s_%s_Value_%d", in.Name.Data(), in.general_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_general_Uncert_Value.Data(), in.general_Uncert[j][k]) );
                    _workspace->factory( Form("expr::general_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.general_Uncertainties[j].Data(), k+1, in.general_Uncertainties[j].Data(),s_general_Uncert_Value.Data(),in.general_Uncertainties[j].Data(), s_general_Uncert_Value.Data()) );
                //cout << "testing error propagation in program - error is: "<< Form("expr::general_Uncert_%s_%s_%d('(1 + Uncert_%s * %s )',Uncert_%s, %s)", in.Name.Data(), in.general_Uncertainties[j].Data(), k+1, in.general_Uncertainties[j].Data(),s_general_Uncert_Value.Data(),in.general_Uncertainties[j].Data(), s_general_Uncert_Value.Data())<<endl;
                }
        //_workspace->Print();
        //(Form("Fatal called!"), " ");
        } // end of NP if statement
        //cout <<"end of loop no. "<<i<<endl;
    } // end of loop over inputs




    // Read in measurements, unfolding & acceptance factors, background yields, Nuisance parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        
        // Signal or Data yields
        for(int j = 0; j < in.Nsig.size(); j++)
            _workspace->factory( Form("Ndata_%s_%d[%f]",in.Name.Data(),j+1,in.Nsig[j]));
        for(int j = 0; j < in.Ndata.size(); j++)
            _workspace->factory( Form("Ndata_%s_%d[%f]",in.Name.Data(),j+1,in.Ndata[j]));
        for(int j = 0; j < in.NsigError.size(); j++)
            _workspace->factory( Form("Ndata_Error_%s_%d[%f]",in.Name.Data(),j+1,in.NsigError[j]));
        
        if( in.NsigError.size() == 0 )
            for(int j = 0; j < in.Ndata.size(); j++)
                _workspace->factory( Form("Ndata_Error_%s_%d[%f]",in.Name.Data(),j+1,sqrt(in.Ndata[j])));
    } // end of loop over inputs
    
    // Unfold the measured yield and add it into the workspace
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Get input and bins
        input in = _input->GetInput()[i]; VecD Bins = _input->GetInput()[_finest_index].Bins; int current_index (0);
        // Here we construct the relation between cross section and observed yield, if we have a coarser
        // granularity we add all fine cross sections until we reach the bin boundary

        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            Str s_sigma  = Form("(Ndata_%s_%d",in.Name.Data(),j+1), s_sigma_Error = Form("Ndata_Error_%s_%d",in.Name.Data(),j+1);
            Str s_Nargs = Form("Ndata_%s_%d",in.Name.Data(),j+1), s_Nargs_Error = Form("Ndata_Error_%s_%d",in.Name.Data(),j+1);
            // Subtract Bkg if needed
            if(in.hasBkg) {
                s_sigma += Form("- Nbkg_%s_%d)",in.Name.Data(),j+1);
                s_Nargs += Form(",Nbkg_%s_%d",in.Name.Data(),j+1);
                } else s_sigma += ")";
            // Now transform the cross section with the fiducial acceptance and the folding factor into an observed yield
            //s_sigma += Form(" / BR_%s / Luminosity_%s * CFactor_%s_%d / FidAcc_%s_%d",in.Name.Data(),in.Name.Data(),in.Name.Data(),j+1,in.Name.Data(),j+1);
            s_sigma += Form(" / BR_%s / Luminosity_%s",in.Name.Data(),in.Name.Data());
            s_sigma_Error += Form(" / BR_%s / Luminosity_%s",in.Name.Data(),in.Name.Data());
            s_Nargs += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            s_Nargs_Error += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            _workspace->factory( Form(" expr::sigma_meas_%s_%d('%s',%s)",in.Name.Data(),j+1,s_sigma.Data(),s_Nargs.Data()) );
            _workspace->factory( Form(" expr::sigma_meas_Error_%s_%d('%s',%s)",in.Name.Data(),j+1,s_sigma_Error.Data(),s_Nargs_Error.Data()) );
        } // end of loop over bins
    } // end of loop over inputs
    
    // Create a prediction for each measured cross section
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Get input and bins
        input in = _input->GetInput()[i]; VecD Bins = _input->GetInput()[_finest_index].Bins; int current_index (0);
        // Here we construct the relation between cross section and observed yield, if we have a coarser
        // granularity we add all fine cross sections until we reach the bin boundary
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            Str s_Nsig = "", s_Nargs = "";
            while( Bins[current_index] < in.Bins[j+1] ) {
                s_Nsig += s_Nsig == "" ? Form("((sigma_%d",1+current_index) : Form(" + sigma_%d",1+current_index);
                s_Nargs += s_Nargs == "" ? Form("sigma_%d",1+current_index) : Form(",sigma_%d",1+current_index);
                current_index++;
            }
            // Now transform the cross section with the fiducial acceptance and the folding factor into an observed yield
            if(_set.IncludeNP) {
            Str s_NP = "", s_NP_args = "";

                for(int k = 0; k < in.BF_Vub_Uncertainties.size(); k++) {
                    s_NP += Form(" * BF_Vub_Uncert_%s_%s_%d ",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("BF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1) : Form(",BF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.BF_Vcb_Uncertainties.size(); k++) {
                    s_NP += Form(" * BF_Vcb_Uncert_%s_%s_%d ",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("BF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1) : Form(",BF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.FF_Vub_Uncertainties.size(); k++) {
                    s_NP += Form(" * FF_Vub_Uncert_%s_%s_%d ",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("FF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1) : Form(",FF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.FF_Vcb_Uncertainties.size(); k++) {
                    s_NP += Form(" * FF_Vcb_Uncert_%s_%s_%d ",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("FF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1) : Form(",FF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.general_Uncertainties.size(); k++) {
                    s_NP += Form(" * general_Uncert_%s_%s_%d ",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("general_Uncert_%s_%s_%d",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1) : Form(",general_Uncert_%s_%s_%d",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1);
                }
                s_Nsig += Form(") * BR_%s * Luminosity_%s %s)",in.Name.Data(),in.Name.Data(),s_NP.Data());
                s_Nargs += Form(",BR_%s,Luminosity_%s,%s",in.Name.Data(),in.Name.Data(),s_NP_args.Data());
            } else {
                s_Nsig += Form(") * BR_%s * Luminosity_%s )",in.Name.Data(),in.Name.Data());
                s_Nargs += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            }

            // Save the expression for eventual later use for the shape combination
            Str s_Nsig_shape = s_Nsig;
            // If we have Background, add it to match the measured yield
            if(in.hasBkg) {
                if(_set.IncludeNP) {
                    Str s_NP = "", s_NP_args = "";
                    for(int k = 0; k < in.Bkg_Uncertainties.size(); k++) {
                        s_NP += Form(" * Bkg_Uncert_%s_%s_%d ",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1);
                        s_NP_args += s_NP_args == "" ? Form("Bkg_Uncert_%s_%s_%d",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1) : Form(",Bkg_Uncert_%s_%s_%d",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1);
                    }
                    s_Nsig += Form(" + Nbkg_%s_%d %s",in.Name.Data(),j+1,s_NP.Data());
                    s_Nargs += Form(",Nbkg_%s_%d,%s",in.Name.Data(),j+1,s_NP_args.Data());
                } else {
                    s_Nsig += Form(" + Nbkg_%s_%d",in.Name.Data(),j+1);
                    s_Nargs += Form(",Nbkg_%s_%d",in.Name.Data(),j+1);
                }
            }
            _workspace->factory( Form(" expr::sigma_%s_%d('%s',%s)",in.Name.Data(),j+1,s_Nsig.Data(),s_Nargs.Data()) );
        } // end of loop over bins
    } // end of loop over inputs
    
    //_workspace->Print();
    //Fatal("Here","");

    // Construct an argument set for the observables and the prediction:
    for(int i = 0; i < _input->GetInput().size(); i++) {
        RooArgSet observables, predictions;
        input in = _input->GetInput()[i];
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            observables.add( *((RooFormulaVar*) _workspace->obj(Form("sigma_meas_%s_%d",in.Name.Data(),j+1))) );
            predictions.add( *((RooFormulaVar*) _workspace->obj(Form("sigma_%s_%d",in.Name.Data(),j+1))) );
        }
        //Now read in the statistical correlation matrices
        ifstream infile;
        string namae;
        Str line;
        //before only infile.open(Form("/home/sduell/bpilnu_new_average/CombinationCode/Correlations/%s_stat.dat", in.Name.Data()));
        cout << "opening file "<< Form("/home/sduell/sources/bpilnu_new_average/CombinationCode/Correlations/%s_stat.dat", in.Name.Data())<<endl;
        if(_set.UseTotCorMat && ( in.Name!="SibidanovBp" && in.Name!="Sibidanov" )) {
            cout<<"Using total correlationmatrix for "<<in.Name.Data()<<endl;
            infile.open(Form("/home/sduell/sources/bpilnu_new_average/CombinationCode/Correlations/%s_tot.dat", in.Name.Data()));
        }
        else infile.open(Form("/home/sduell/sources/bpilnu_new_average/CombinationCode/Correlations/%s_stat.dat", in.Name.Data()));

        vector<VecD> a;
        int n=0;
        cout << "begin reading file"<<endl;
        while(!infile.eof()){
            getline(infile,namae);
            line=Form("%s", namae.c_str());
            a.push_back(VecD());
            a[n]=VectorizeD(line," ");
            cout << "line: "<<line<<endl;
            cout << "element: "<<a[n][1]<<endl;
            n++;
        }
//Comment this line in to modify the cov matrix with a factor to account for the change in tauB0/tauB+
        //if(in.Name.Data()=="SibidanovBp") covfactor=1.0;
        cout << "file successfully read"<<endl;

        TMatrixDSym *CorMatrix = new TMatrixDSym(a[0].size());
        //normal Correlation matrix
        for(int j=0; j<a.size(); j++){
            for(int k=0; k<a[j].size(); k++){
                    (*CorMatrix)(j,k)=a[j][k];
                    //cout << "Matrix is: "<<a[j][k]<<endl;
            }   
        }


        cout << "Correlation Matrix: "<<endl;
        (CorMatrix)->Print();

        //Now build the (statistical) covariance matrix from the correlation matrix
        VecD errs;
        
        for(int j = 0; j < in.NsigError.size(); j++)
            errs.push_back(in.NsigError[j]);
        
        TMatrixDSym *ErrMatrix = new TMatrixDSym(errs.size());
        
        for(int j=0; j<a.size(); j++){
            (*ErrMatrix)(j,j)=errs[j];
        }

        TMatrixDSym *CovMatrix = new TMatrixDSym(*CorMatrix);
        CovMatrix->SimilarityT(*ErrMatrix);
        CovMatrix->Print();
        cout << "Covariance matrix built"<<endl;
        // Now build the total likelihood
        RooMultiVarGaussian *pdf = new RooMultiVarGaussian(Form("pdf_%s", in.Name.Data()), Form("pdf_%s",in.Name.Data()), predictions, observables, (*CovMatrix));
        _workspace->import(*pdf);
        cout << "constructed pdf"<<endl;
    }
    

    // Multiply in the NP if requested
    Str NPComponentNames = "";
    if(_set.IncludeNP) {
        for(int j = 0; j < _sv_source_const.size(); j++)
            NPComponentNames += Form(",%s",_sv_source_const[j].Data());
    }
    
    // Now let us construct the likelihood function, we loop over all inputs and combine pdfs by multiplying them one by one
    Str pdfTotName = "";
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Grab the input
        input in = _input->GetInput()[i];
        //Str pdfComponentNames = "";
        pdfTotName += i == 0 ? Form("pdf_%s", in.Name.Data()) : Form(",pdf_%s", in.Name.Data());
        
        // Build pdf per decay channel (with NPs if they are included)
        if( !_workspace->pdf( Form("pdf_%s_NPs(%s)",in.Name.Data(), (Form("pdf_%s",in.Name.Data())+NPComponentNames).Data()) ) )
            _workspace->factory( Form("PROD::pdf_%s_NPs(%s)", in.Name.Data(), (Form("pdf_%s",in.Name.Data())+NPComponentNames).Data()) );
        
    }
    
    // Sum PDFs of all components and we are done!
    if( !_workspace->pdf( Form("pdf_tot(%s)", (pdfTotName+NPComponentNames).Data()) ) )
        _workspace->factory( Form("PROD::pdf_tot(%s)",(pdfTotName+NPComponentNames).Data()) );

    
    Printf("Print workspace for Multi combination:\n\n");
    _workspace->Print("v");
    Printf("\n\n");
    
    //Fatal("here","");

}

// Create the combination workspace
void HiggsWorkspace::SetupWorkspace() {
    
    // Determine the highest granularity and add bin cross sections into the workspace as free parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        if( _max_bins < in.Bins.size() ) {
            _max_bins = in.Bins.size(); _finest_index = i;
        }
    }
    
    // If we do a shape combination, we have one degree of freedom less
    if(_set.CombinationMethod == "shape")
        _workspace->factory( Form("sigma_%d[%f,%f,%f]",1,1.0,1.0,1.0) );
    else
        _workspace->factory( Form("sigma_%d[%f,%f,%f]",1,_set.CrossSectionOptions[0],_set.CrossSectionOptions[1],_set.CrossSectionOptions[2]) );
    // Setup the remaining parameters of interests
    for(int i = 1; i < _max_bins-1; i++)
        _workspace->factory( Form("sigma_%d[%f,%f,%f]",i+1,_set.CrossSectionOptions[0],_set.CrossSectionOptions[1],_set.CrossSectionOptions[2]) );
    
    // Read in measurements, unfolding & acceptance factors, background yields, Nuisance parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        // Branching Fraction
        _workspace->factory( Form("BR_%s[%f]",in.Name.Data(),in.BR));
        // Luminosity
        _workspace->factory( Form("Luminosity_%s[%f]",in.Name.Data(),in.Luminosity));
        // Bkg yields
        for(int j = 0; j < in.Nbkg.size(); j++) {
            _workspace->factory( Form("Nbkg_%s_%d[%f]",in.Name.Data(),j+1,in.Nbkg[j]));
        }
        // Include the three classes of Nuisance parameters and constraints (CFactors, Acceptance, Nbkg) into workspace if required
        if(_set.IncludeNP) {
            cout <<"IncludeNP in Workspace activated"<<endl;
            // ------------------------------------------------------------------------------------
            // Constraints & NP

            for(int j = 0; j < in.Bkg_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.Bkg_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.Bkg_Uncertainties[j].Data(),in.Bkg_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.Bkg_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.BF_Vub_Uncertainties.size(); j++) {
                // Check if constrained already exists
                //cout << "object tbc is: "<<Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vub_Uncertainties[j].Data(),in.BF_Vub_Uncertainties[j].Data())<<endl;
                if( !_workspace->obj( Form("Uncert_%s_const", in.BF_Vub_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vub_Uncertainties[j].Data(),in.BF_Vub_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.BF_Vub_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.BF_Vcb_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.BF_Vcb_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vcb_Uncertainties[j].Data(),in.BF_Vcb_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.BF_Vcb_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.FF_Vub_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.FF_Vub_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.FF_Vub_Uncertainties[j].Data(),in.FF_Vub_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.FF_Vub_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.FF_Vcb_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.FF_Vcb_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.FF_Vcb_Uncertainties[j].Data(),in.FF_Vcb_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.FF_Vcb_Uncertainties[j].Data()) );
                }
            }
            for(int j = 0; j < in.general_Uncertainties.size(); j++) {
                // Check if constrained already exists
                if( !_workspace->obj( Form("Uncert_%s_const", in.general_Uncertainties[j].Data()) ) ) {
                    _workspace->factory( Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.general_Uncertainties[j].Data(),in.general_Uncertainties[j].Data()) );
                    _sv_source_const.push_back( Form("Uncert_%s_const",in.general_Uncertainties[j].Data()) );
                }
            }
            // ------------------------------------------------------------------------------------
            // Errors
            for(int j = 0; j < in.Bkg_Uncertainties.size(); j++)
                for(int k = 0; k < in.Bkg_Uncert[j].size(); k++) {
                    Str s_Bkg_Uncert_Value = Form("Bkg_Uncert_%s_%s_Value_%d", in.Name.Data(), in.Bkg_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_Bkg_Uncert_Value.Data(), in.Bkg_Uncert[j][k]) );
                    _workspace->factory( Form("expr::Bkg_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.Bkg_Uncertainties[j].Data(), k+1, in.Bkg_Uncertainties[j].Data(),s_Bkg_Uncert_Value.Data(),in.Bkg_Uncertainties[j].Data(), s_Bkg_Uncert_Value.Data()) );
                }
            for(int j = 0; j < in.BF_Vub_Uncertainties.size(); j++)
                for(int k = 0; k < in.BF_Vub_Uncert[j].size(); k++) {
                    //cout << "error tbc is: "<<Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.BF_Vub_Uncertainties[j].Data(),in.BF_Vub_Uncertainties[j].Data())<<endl;
                    Str s_BF_Vub_Uncert_Value = Form("BF_Vub_Uncert_%s_%s_Value_%d", in.Name.Data(), in.BF_Vub_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_BF_Vub_Uncert_Value.Data(), in.BF_Vub_Uncert[j][k]) );
                    _workspace->factory( Form("expr::BF_Vub_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.BF_Vub_Uncertainties[j].Data(), k+1, in.BF_Vub_Uncertainties[j].Data(),s_BF_Vub_Uncert_Value.Data(),in.BF_Vub_Uncertainties[j].Data(), s_BF_Vub_Uncert_Value.Data()) ); 
                }
            for(int j = 0; j < in.BF_Vcb_Uncertainties.size(); j++)
                for(int k = 0; k < in.BF_Vcb_Uncert[j].size(); k++) {
                    Str s_BF_Vcb_Uncert_Value = Form("BF_Vcb_Uncert_%s_%s_Value_%d", in.Name.Data(), in.BF_Vcb_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_BF_Vcb_Uncert_Value.Data(), in.BF_Vcb_Uncert[j][k]) );
                    _workspace->factory( Form("expr::BF_Vcb_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.BF_Vcb_Uncertainties[j].Data(), k+1, in.BF_Vcb_Uncertainties[j].Data(),s_BF_Vcb_Uncert_Value.Data(),in.BF_Vcb_Uncertainties[j].Data(), s_BF_Vcb_Uncert_Value.Data()) );
                }
                //cout << "After BF Uncertainties!"<<endl;
                //cout << "Next is Vub FF Uncertainty with size "<< in.FF_Vub_Uncertainties.size() << endl;
            for(int j = 0; j < in.FF_Vub_Uncertainties.size(); j++){
                for(int k = 0; k < in.FF_Vub_Uncert[j].size(); k++) {
                    //cout << "FF error tbc is: "<<Form("Gaussian::Uncert_%s_const(Uncert_%s[0,-10,10],0,1)",in.FF_Vub_Uncertainties[j].Data(),in.FF_Vub_Uncertainties[j].Data())<<endl;
                    Str s_FF_Vub_Uncert_Value = Form("FF_Vub_Uncert_%s_%s_Value_%d", in.Name.Data(), in.FF_Vub_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_FF_Vub_Uncert_Value.Data(), in.FF_Vub_Uncert[j][k]) );
                    _workspace->factory( Form("expr::FF_Vub_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.FF_Vub_Uncertainties[j].Data(), k+1, in.FF_Vub_Uncertainties[j].Data(),s_FF_Vub_Uncert_Value.Data(),in.FF_Vub_Uncertainties[j].Data(), s_FF_Vub_Uncert_Value.Data()) );
                }}
                //cout <<"test to see if NP correctly initiated "<<endl;
            for(int j = 0; j < in.FF_Vcb_Uncertainties.size(); j++)
                for(int k = 0; k < in.FF_Vcb_Uncert[j].size(); k++) {
                    Str s_FF_Vcb_Uncert_Value = Form("FF_Vcb_Uncert_%s_%s_Value_%d", in.Name.Data(), in.FF_Vcb_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_FF_Vcb_Uncert_Value.Data(), in.FF_Vcb_Uncert[j][k]) );
                    _workspace->factory( Form("expr::FF_Vcb_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.FF_Vcb_Uncertainties[j].Data(), k+1, in.FF_Vcb_Uncertainties[j].Data(),s_FF_Vcb_Uncert_Value.Data(),in.FF_Vcb_Uncertainties[j].Data(), s_FF_Vcb_Uncert_Value.Data()) );
                }

            for(int j = 0; j < in.general_Uncertainties.size(); j++)
                for(int k = 0; k < in.general_Uncert[j].size(); k++) {
                    Str s_general_Uncert_Value = Form("general_Uncert_%s_%s_Value_%d", in.Name.Data(), in.general_Uncertainties[j].Data(), k+1 );
                    _workspace->factory( Form("%s[%f]", s_general_Uncert_Value.Data(), in.general_Uncert[j][k]) );
                    _workspace->factory( Form("expr::general_Uncert_%s_%s_%d('(1 + (Uncert_%s / 100.) * %s )',Uncert_%s, %s)", in.Name.Data(), in.general_Uncertainties[j].Data(), k+1, in.general_Uncertainties[j].Data(),s_general_Uncert_Value.Data(),in.general_Uncertainties[j].Data(), s_general_Uncert_Value.Data()) );
                cout << "testing error propagation in program - error is: "<< Form("expr::general_Uncert_%s_%s_%d('(1 + Uncert_%s * %s )',Uncert_%s, %s)", in.Name.Data(), in.general_Uncertainties[j].Data(), k+1, in.general_Uncertainties[j].Data(),s_general_Uncert_Value.Data(),in.general_Uncertainties[j].Data(), s_general_Uncert_Value.Data())<<endl;
                }
        //_workspace->Print();
        //(Form("Fatal called!"), " ");
        } // end of NP if statement
        cout <<"end of loop no. "<<i<<endl;
    } // end of loop over inputs
    cout << "after input loop"<<endl;
    // Construct expressions relevant if we wan to do a shape comparison
    if(_set.CombinationMethod == "shape") {
        // Loop over all inputs
        for(int i = 0; i < _input->GetInput().size(); i++) {
            input in = _input->GetInput()[i]; VecD Bins = _input->GetInput()[i].Bins; int current_index (0);
            Str sigma_obs_tot = ""; Str sigma_obs_tot_args = "";
            for(int j = 0; j < Bins.size()-1; j++) {
                Str sigma_obs_i = Form("((%f", !in.hasBkg ? in.Nsig[j] : in.Ndata[j]);
                // Subtract background if needed
                if(in.hasBkg) {
                    if(_set.IncludeNP) {
                        Str s_NP = "";
                        for(int k = 0; k < in.Bkg_Uncertainties.size(); k++) {
                            s_NP += Form(" * Bkg_Uncert_%s_%s_%d ",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1);
                            sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("Bkg_Uncert_%s_%s_%d",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1) : Form(",Bkg_Uncert_%s_%s_%d",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1);
                        }
                        sigma_obs_i += Form(" - Nbkg_%s_%d %s", in.Name.Data(),j+1,s_NP.Data());
                        sigma_obs_tot_args += Form(",Nbkg_%s_%d",in.Name.Data(),j+1);
                    } else {
                        sigma_obs_i += Form(" - Nbkg_%s_%d", in.Name.Data(),j+1);
                        sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("Nbkg_%s_%d",in.Name.Data(),j+1) : Form(",Nbkg_%s_%d",in.Name.Data(),j+1);
                    }
                }
                // Now unfold the observation into a cross section
                Str s_NP = "", s_NP_args = "";
                if(_set.IncludeNP) {

                    for(int k = 0; k < in.BF_Vub_Uncertainties.size(); k++) {
                        s_NP += Form(" / BF_Vub_Uncert_%s_%s_%d ",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1);
                        sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("BF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1) : Form(",BF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1);
                    }
                    for(int k = 0; k < in.BF_Vcb_Uncertainties.size(); k++) {
                        s_NP += Form(" / BF_Vcb_Uncert_%s_%s_%d ",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1);
                        sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("BF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1) : Form(",BF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1);
                    }
                    for(int k = 0; k < in.FF_Vub_Uncertainties.size(); k++) {
                        s_NP += Form(" / FF_Vub_Uncert_%s_%s_%d ",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1);
                        sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("FF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1) : Form(",FF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1);
                    }
                    for(int k = 0; k < in.FF_Vcb_Uncertainties.size(); k++) {
                        s_NP += Form(" / FF_Vcb_Uncert_%s_%s_%d ",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1);
                        sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("FF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1) : Form(",FF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1);
                    }
                    for(int k = 0; k < in.general_Uncertainties.size(); k++) {
                        s_NP += Form(" / general_Uncert_%s_%s_%d ",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1);
                        sigma_obs_tot_args += sigma_obs_tot_args == "" ? Form("general_Uncert_%s_%s_%d",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1) : Form(",general_Uncert_%s_%s_%d",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1);
                    }

                }
                sigma_obs_i += Form(") / BR_%s / Luminosity_%s  %s)",in.Name.Data(),in.Name.Data(),s_NP.Data());
                //sigma_obs_tot_args += Form(",CFactor_%s_%d,FidAcc_%s_%d",in.Name.Data(),j+1,in.Name.Data(),j+1);
                sigma_obs_tot_args += ""; //Check for the other uncertainties
                // Add to the total observed cross section
                if(j != 0) sigma_obs_tot += " + "; sigma_obs_tot += sigma_obs_i;
            } // End loop over all bins
            sigma_obs_tot_args += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            // Store value in workspace
            _workspace->factory( Form(" expr::sigma_obs_tot_%s('%s',%s)",in.Name.Data(),sigma_obs_tot.Data(),sigma_obs_tot_args.Data()) );
            // Now construct the new expressions for the floated cross section fractions
            for(int i = 0; i < _max_bins-1; i++) {
                Str s_sigma = Form("sigma_%d/(",i+1), s_sigma_Nargs = "";
                for(int k = 0; k < _max_bins-1 ; k++) {
                    s_sigma += k == 0 ? Form("sigma_%d",k+1) : Form(" + sigma_%d",k+1);
                    s_sigma_Nargs += s_sigma_Nargs == "" ? Form("sigma_%d",k+1) : Form(",sigma_%d",k+1);
                }
                s_sigma += Form(")*sigma_obs_tot_%s",in.Name.Data()); s_sigma_Nargs += Form(",sigma_obs_tot_%s",in.Name.Data());
                _workspace->factory( Form(" expr::sigma_norm_%s_%d('%s',%s)",in.Name.Data(),i+1,s_sigma.Data(),s_sigma_Nargs.Data()) );
            }
            
        } // End loop over all inputs
    }  // End if for shape combination  
    cout << "after shape combination"<<endl;
    
    // Connect cross sections with each measurement
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Get input and bins
        input in = _input->GetInput()[i]; VecD Bins = _input->GetInput()[_finest_index].Bins; int current_index (0);
        // Here we construct the relation between cross section and observed yield, if we have a coarser
        // granularity we add all fine cross sections until we reach the bin boundary
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            Str s_Nsig = "", s_Nargs = "";
            if(_set.CombinationMethod == "shape") {
             while( Bins[current_index] < in.Bins[j+1] ) {
                 s_Nsig += s_Nsig == "" ? Form("((sigma_norm_%s_%d",in.Name.Data(),1+current_index) : Form(" + sigma_norm_%s_%d",in.Name.Data(),1+current_index);
                 s_Nargs += s_Nargs == "" ? Form("sigma_norm_%s_%d",in.Name.Data(),1+current_index) : Form(",sigma_norm_%s_%d",in.Name.Data(),1+current_index);
                 current_index++;
             }            
            } else {
             while( Bins[current_index] < in.Bins[j+1] ) {
                 s_Nsig += s_Nsig == "" ? Form("((sigma_%d",1+current_index) : Form(" + sigma_%d",1+current_index);
                 s_Nargs += s_Nargs == "" ? Form("sigma_%d",1+current_index) : Form(",sigma_%d",1+current_index);
                 current_index++;
             }
            }
            // Now transform the cross section with the fiducial acceptance and the folding factor into an observed yield
            if(_set.IncludeNP) {
                Str s_NP = "", s_NP_args = "";
            
                for(int k = 0; k < in.BF_Vub_Uncertainties.size(); k++) {
                    s_NP += Form(" * BF_Vub_Uncert_%s_%s_%d ",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("BF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1) : Form(",BF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vub_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.BF_Vcb_Uncertainties.size(); k++) {
                    s_NP += Form(" * BF_Vcb_Uncert_%s_%s_%d ",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("BF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1) : Form(",BF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.BF_Vcb_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.FF_Vub_Uncertainties.size(); k++) {
                    s_NP += Form(" * FF_Vub_Uncert_%s_%s_%d ",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("FF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1) : Form(",FF_Vub_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vub_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.FF_Vcb_Uncertainties.size(); k++) {
                    s_NP += Form(" * FF_Vcb_Uncert_%s_%s_%d ",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("FF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1) : Form(",FF_Vcb_Uncert_%s_%s_%d",in.Name.Data(),in.FF_Vcb_Uncertainties[k].Data(),j+1);
                }
                for(int k = 0; k < in.general_Uncertainties.size(); k++) {
                    s_NP += Form(" * general_Uncert_%s_%s_%d ",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1);
                    s_NP_args += s_NP_args == "" ? Form("general_Uncert_%s_%s_%d",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1) : Form(",general_Uncert_%s_%s_%d",in.Name.Data(),in.general_Uncertainties[k].Data(),j+1);
                }
                s_Nsig += Form(") * BR_%s * Luminosity_%s %s)",in.Name.Data(),in.Name.Data(),s_NP.Data());
                s_Nargs += Form(",BR_%s,Luminosity_%s,%s",in.Name.Data(),in.Name.Data(),s_NP_args.Data());
            } else {
                s_Nsig += Form(") * BR_%s * Luminosity_%s )",in.Name.Data(),in.Name.Data());
                s_Nargs += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            }
            // Save the expression for eventual later use for the shape combination
            Str s_Nsig_shape = s_Nsig;
            // If we have Background, add it to match the measured yield
            if(in.hasBkg) {
                if(_set.IncludeNP) {
                    Str s_NP = "", s_NP_args = "";
                    for(int k = 0; k < in.Bkg_Uncertainties.size(); k++) {
                        s_NP += Form(" * Bkg_Uncert_%s_%s_%d ",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1);
                        s_NP_args += s_NP_args == "" ? Form("Bkg_Uncert_%s_%s_%d",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1) : Form(",Bkg_Uncert_%s_%s_%d",in.Name.Data(),in.Bkg_Uncertainties[k].Data(),j+1);
                    }
                    s_Nsig += Form(" + Nbkg_%s_%d %s",in.Name.Data(),j+1,s_NP.Data());
                    s_Nargs += Form(",Nbkg_%s_%d,%s",in.Name.Data(),j+1,s_NP_args.Data());
                } else {
                    s_Nsig += Form(" + Nbkg_%s_%d",in.Name.Data(),j+1);
                    s_Nargs += Form(",Nbkg_%s_%d",in.Name.Data(),j+1);
                }
            }
            _workspace->factory( Form(" expr::Ndata_%s_%d('%s',%s)",in.Name.Data(),j+1,s_Nsig.Data(),s_Nargs.Data()) );
        } // end of loop over bins
    } // end of loop over inputs
        
    // Now we are ready to read in the actual measurements
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        // Include all measurements into the workspace
        for(int j = 0; j < in.Nsig.size(); j++) {
            // Use likelihood contour and convert it into a PDF
            if(in.hasLikelihoodContour == true) {
                // Create a histogram from the contour and convert the -2 Log (L) into a proper Likelihood
                double *x = in.nll_contour->GetX(); VecD bins; double shift = (x[1] - x[0])/2;
                for(int ibins = 0; ibins < in.nll_contour->GetN(); ibins++)
                    bins.push_back( (x[ibins] - shift)*in.NsigSM[j] );
                TH1 *h = new TH1D("","",bins.size()-1,&bins[0]), *c = new TH1D("","",bins.size()-1,&bins[0]);
                for(int ibins = 0; ibins < in.nll_contour->GetN(); ibins++) {
                    h->SetBinContent(ibins+1, TMath::Exp(-0.5*in.nll_contour->Eval(h->GetXaxis()->GetBinCenter(ibins+1)/in.NsigSM[j])));
                    c->SetBinContent(ibins+1, in.nll_contour->Eval(h->GetXaxis()->GetBinCenter(ibins+1)/in.NsigSM[j]));
                }
                // Debug plots
                // TCanvas can; c->Draw(); can.Print("c_test.pdf"); h->Draw(); can.Print("h_test.pdf");
                h->SetBinContent(0, TMath::Exp(-0.5*in.nll_contour->Eval(h->GetXaxis()->GetBinCenter(1)/in.NsigSM[j])));
                c->SetBinContent(0, in.nll_contour->Eval(h->GetXaxis()->GetBinCenter(1)/in.NsigSM[j]));
                // Now construct a RooDataHist object and a RooHistPdf and import them into the workspace
                _workspace->factory( Form("Ndata_cons_%s_%d[%f]",in.Name.Data(),j+1,in.NsigSM[j]));
                RooDataHist *dh_Ndata = new RooDataHist( Form("dh_Ndata_%s_%d",in.Name.Data(),j+1), Form("dh_Ndata_%s_%d",in.Name.Data(),j+1), RooArgList(*((RooAbsReal*)_workspace->obj(Form("Ndata_cons_%s_%d",in.Name.Data(),j+1)))),Import(*h));
                RooDataHist *dc_Ndata = new RooDataHist( Form("dc_Ndata_%s_%d",in.Name.Data(),j+1), Form("dc_Ndata_%s_%d",in.Name.Data(),j+1), RooArgList(*((RooAbsReal*)_workspace->obj(Form("Ndata_cons_%s_%d",in.Name.Data(),j+1)))),Import(*c));
                RooHistPdf  *rh_Ndata = new RooHistPdf( Form("rh_Ndata_%s_%d",in.Name.Data(),j+1), Form("rh_Ndata_%s_%d",in.Name.Data(),j+1), *((RooAbsReal*)_workspace->obj(Form("Ndata_cons_%s_%d",in.Name.Data(),j+1))),*dh_Ndata,2);
                RooHistPdf  *rc_Ndata = new RooHistPdf( Form("rc_Ndata_%s_%d",in.Name.Data(),j+1), Form("rc_Ndata_%s_%d",in.Name.Data(),j+1), *((RooAbsReal*)_workspace->obj(Form("Ndata_cons_%s_%d",in.Name.Data(),j+1))),*dc_Ndata,2);
                // Tell RooFit to not change the normalization and import them into the workspace
                rh_Ndata->setUnitNorm(true); rc_Ndata->setUnitNorm(true);
                _workspace->import(*rh_Ndata); _workspace->import(*rc_Ndata);
                // Ok, now replace the fundamental RooRealVar with our non-fundamental RooFormulaVar
                _workspace->factory( Form("EDIT::pdf_Ndata_%s_%d(rh_Ndata_%s_%d,Ndata_cons_%s_%d=Ndata_%s_%d)",in.Name.Data(),j+1,in.Name.Data(),j+1,in.Name.Data(),j+1,in.Name.Data(),j+1));
                _workspace->factory( Form("EDIT::ll_Ndata_%s_%d(rc_Ndata_%s_%d,Ndata_cons_%s_%d=Ndata_%s_%d)",in.Name.Data(),j+1,in.Name.Data(),j+1,in.Name.Data(),j+1,in.Name.Data(),j+1));
            } else {
                in.hasAsymUncert == false ? _workspace->factory( Form("Gaussian::pdf_Ndata_%s_%d(%f,Ndata_%s_%d,%f)",in.Name.Data(),j+1,in.Nsig[j],in.Name.Data(),j+1, in.NsigError[j]) )
                : _workspace->factory( Form("RooBifurGauss::pdf_Ndata_%s_%d(%f,Ndata_%s_%d,%f,%f)",in.Name.Data(),j+1,in.Nsig[j],in.Name.Data(),j+1, in.NsigErrorLow[j],in.NsigErrorHigh[j]) );
                // _workspace->factory( Form("Gaussian::pdf_Ndata_%s_%d(Ndata_%s_%d,%f,%f)",in.Name.Data(),j+1,in.Name.Data(),j+1, in.Nsig[j], in.NsigError[j]) );
            }
        }
        for(int j = 0; j < in.Ndata.size(); j++) {
            if(_set.NdataModel == "Gaussian")
                _workspace->factory( Form("Gaussian::pdf_Ndata_%s_%d(%f,Ndata_%s_%d,%f)",in.Name.Data(),j+1,in.Ndata[j],in.Name.Data(),j+1, sqrt(in.Ndata[j]) ) );
            // _workspace->factory( Form("Gaussian::pdf_Ndata_%s_%d(Ndata_%s_%d,%f,%f)",in.Name.Data(),j+1,in.Name.Data(),j+1, in.Ndata[j], sqrt(in.Ndata[j]) ) );
            else if (_set.NdataModel == "Gamma")
                _workspace->factory( Form("Gamma::pdf_Ndata_%s_%d(%f,expr::Ndata_%s_%d_plus('Ndata_%s_%d+1',Ndata_%s_%d),1,0)",in.Name.Data(),j+1,in.Ndata[j],in.Name.Data(),j+1,in.Name.Data(),j+1,in.Name.Data(),j+1) );
            else
                _workspace->factory( Form("Poisson::pdf_Ndata_%s_%d(%f,Ndata_%s_%d)",in.Name.Data(),j+1,in.Ndata[j],in.Name.Data(),j+1) );
            //_workspace->factory( Form("Poisson::pdf_Ndata_%s_%d(Ndata_%s_%d,%f)",in.Name.Data(),j+1,in.Name.Data(),j+1, in.Ndata[j]) );
            
        }
    } // End Component loop
    
    // Multiply in the NP if requested
    Str NPComponentNames = "";
    if(_set.IncludeNP) {
        for(int j = 0; j < _sv_source_const.size(); j++)
            NPComponentNames += Form(",%s",_sv_source_const[j].Data());
    }
    
    // Now let us construct the likelihood function, we loop over all inputs and combine pdfs by multiplying them one by one
    Str pdfTotName = "";
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Grab the input
        input in = _input->GetInput()[i];
        Str pdfComponentNames = "";
        for(int j = 0; j < in.Nsig.size(); j++)
            pdfComponentNames += j == 0 ? Form("pdf_Ndata_%s_%d",in.Name.Data(),j+1) : Form(",pdf_Ndata_%s_%d",in.Name.Data(),j+1);
        // Build name of combined pdf
        for(int j = 0; j < in.Ndata.size(); j++)
            pdfComponentNames += j == 0 ? Form("pdf_Ndata_%s_%d",in.Name.Data(),j+1) : Form(",pdf_Ndata_%s_%d",in.Name.Data(),j+1);
        cout << pdfComponentNames << endl;
        if( !_workspace->pdf( Form("pdf_%s(%s)",in.Name.Data(), pdfComponentNames.Data()) ) )
            _workspace->factory( Form("PROD::pdf_%s(%s)", in.Name.Data(), (pdfComponentNames).Data()) );
        pdfTotName += i == 0 ? Form("pdf_%s", in.Name.Data()) : Form(",pdf_%s", in.Name.Data());
        
        // Build pdf per decay channel (with NPs if they are included)
        if( !_workspace->pdf( Form("pdf_%s_NPs(%s)",in.Name.Data(), (pdfComponentNames+NPComponentNames).Data()) ) )
            _workspace->factory( Form("PROD::pdf_%s_NPs(%s)", in.Name.Data(), (pdfComponentNames+NPComponentNames).Data()) );
        
    }
    
    // Sum PDFs of all components and we are done!
    if( !_workspace->pdf( Form("pdf_tot(%s)", pdfTotName.Data()) ) )
        _workspace->factory( Form("PROD::pdf_tot(%s)",(pdfTotName+NPComponentNames).Data()) );
    
    // Print workspace information into the console
    Printf("Print workspace for combination:\n\n");
    _workspace->Print("v");
    Printf("\n\n");
    //Fatal("here","");
}

// Create the simplified combination workspace identical with the chi2 combination
void HiggsWorkspace::SetupSimplifiedWorkspace() {
    
    // Determine the highest granularity and add bin cross sections into the workspace as free parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        if( _max_bins < in.Bins.size() ) {
            _max_bins = in.Bins.size(); _finest_index = i;
        }
    }
    for(int i = 0; i < _max_bins; i++) {
        _workspace->factory( Form("sigma_%d[%f,%f,%f]",i+1,_set.CrossSectionOptions[0],_set.CrossSectionOptions[1],_set.CrossSectionOptions[2]) );
    }
    
    // Read in measurements, unfolding & acceptance factors, background yields, Nuisance parameters
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        // Branching Fraction
        _workspace->factory( Form("BR_%s[%f]",in.Name.Data(),in.BR));
        // Luminosity
        _workspace->factory( Form("Luminosity_%s[%f]",in.Name.Data(),in.Luminosity));
        // Unfolding factor
               // Bkg yields
        for(int j = 0; j < in.Nbkg.size(); j++)
            _workspace->factory( Form("Nbkg_%s_%d[%f]",in.Name.Data(),j+1,in.Nbkg[j]));
        // Signal or Data yields
        for(int j = 0; j < in.Nsig.size(); j++)
            _workspace->factory( Form("Ndata_%s_%d[%f]",in.Name.Data(),j+1,in.Nsig[j]));
        for(int j = 0; j < in.Ndata.size(); j++)
            _workspace->factory( Form("Ndata_%s_%d[%f]",in.Name.Data(),j+1,in.Ndata[j]));
        for(int j = 0; j < in.NsigError.size(); j++)
            _workspace->factory( Form("Ndata_Error_%s_%d[%f]",in.Name.Data(),j+1,in.NsigError[j]));
        
        if( in.NsigError.size() == 0 )
            for(int j = 0; j < in.Ndata.size(); j++)
                _workspace->factory( Form("Ndata_Error_%s_%d[%f]",in.Name.Data(),j+1,sqrt(in.Ndata[j])));
    } // end of loop over inputs
    
    // Unfold the measured yield and add it into the workspace
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Get input and bins
        input in = _input->GetInput()[i]; VecD Bins = _input->GetInput()[_finest_index].Bins; int current_index (0);
        // Here we construct the relation between cross section and observed yield, if we have a coarser
        // granularity we add all fine cross sections until we reach the bin boundary
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            Str s_sigma  = Form("(Ndata_%s_%d",in.Name.Data(),j+1), s_sigma_Error = Form("Ndata_Error_%s_%d",in.Name.Data(),j+1);
            Str s_Nargs = Form("Ndata_%s_%d",in.Name.Data(),j+1), s_Nargs_Error = Form("Ndata_Error_%s_%d",in.Name.Data(),j+1);
            // Subtract Bkg if needed
            if(in.hasBkg) {
                s_sigma += Form("- Nbkg_%s_%d)",in.Name.Data(),j+1);
                s_Nargs += Form(",Nbkg_%s_%d",in.Name.Data(),j+1);
                } else s_sigma += ")";
            // Now transform the cross section with the fiducial acceptance and the folding factor into an observed yield
            s_sigma += Form(" / BR_%s / Luminosity_%s",in.Name.Data(),in.Name.Data());
            s_sigma_Error += Form(" / BR_%s / Luminosity_%s",in.Name.Data(),in.Name.Data());
            s_Nargs += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            s_Nargs_Error += Form(",BR_%s,Luminosity_%s",in.Name.Data(),in.Name.Data());
            _workspace->factory( Form(" expr::sigma_meas_%s_%d('%s',%s)",in.Name.Data(),j+1,s_sigma.Data(),s_Nargs.Data()) );
            _workspace->factory( Form(" expr::sigma_meas_Error_%s_%d('%s',%s)",in.Name.Data(),j+1,s_sigma_Error.Data(),s_Nargs_Error.Data()) );
        } // end of loop over bins
    } // end of loop over inputs
    
    // Create a prediction for each measured cross section
    for(int i = 0; i < _input->GetInput().size(); i++) {
        // Get input and bins
        input in = _input->GetInput()[i]; VecD Bins = _input->GetInput()[_finest_index].Bins; int current_index (0);
        // Here we construct the relation between cross section and observed yield, if we have a coarser
        // granularity we add all fine cross sections until we reach the bin boundary
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            Str s_Nsig = "", s_Nargs = "";
            while( Bins[current_index] < in.Bins[j+1] ) {
                s_Nsig += s_Nsig == "" ? Form("sigma_%d",1+current_index) : Form(" + sigma_%d",1+current_index);
                s_Nargs += s_Nargs == "" ? Form("sigma_%d",1+current_index) : Form(",sigma_%d",1+current_index);
                current_index++;
            }
            _workspace->factory( Form(" expr::sigma_%s_%d('%s',%s)",in.Name.Data(),j+1,s_Nsig.Data(),s_Nargs.Data()) );
        } // end of loop over bins
    } // end of loop over inputs
    
    // Construct an argument set for the observables and the prediction:
    RooArgSet observables, predictions;
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            observables.add( *((RooFormulaVar*) _workspace->obj(Form("sigma_meas_%s_%d",in.Name.Data(),j+1))) );
            predictions.add( *((RooFormulaVar*) _workspace->obj(Form("sigma_%s_%d",in.Name.Data(),j+1))) );
        }
    }
    
    // Now build the total likelihood
    RooMultiVarGaussian *pdf = new RooMultiVarGaussian("pdf_tot", "pdf_tot", predictions, observables, (*_TotcovMatrix));
    _workspace->import(*pdf);
    
    Printf("Print workspace for combination:\n\n");
    _workspace->Print("v");
    Printf("\n\n");
    
    
}

// Create a demo workspace
void HiggsWorkspace::SetupDemoWorkspace() {
    
    _workspace->factory("Gaussian::pdf_tot(a[0],a_mu[0,-10,10],a_sig[1])");
    
}

// Merge and return covariance
TMatrixDSym HiggsWorkspace::MergeCovariance(TMatrixDSym first, TMatrixDSym second) {
    for(int i = 0; i < first.GetNrows(); i++)
        for(int j = 0; j < first.GetNrows(); j++)
            second(i,j) = first(i,j) == 0 ? second(i,j) : first(i,j) ;
    return second;
}

// Construct statistical, systematic and overall covariance
void HiggsWorkspace::ConstructCovariance() {
        cout << "Start building cov matrix"<<endl;
    // Now let's build a covariance matrix for the statistical errors
    _StatcovMatrix =  new TMatrixDSym( _n_meas ); int j_bin(0);
    for(int i = 0; i < _input->GetInput().size(); i++) {
        input in = _input->GetInput()[i];
        //cout << "input "<< i << " is "<< _input->GetInput()[i]<< endl;
        for(int j = 0; j < in.Bins.size()-1 ; j++) {
            double sigma_j = in.hasBkg ? sqrt(in.Ndata[j]) / in.BR / in.Luminosity
            : in.NsigError[j] / in.BR / in.Luminosity;
            
            
            // Variance is Error^2
            (*_StatcovMatrix)(j_bin,j_bin) = pow( sigma_j, 2.);
            j_bin++;
        }
    }
    cout << "Cov matrix built"<<endl;
    // Now for the systematic covariance (if requested)
    _SyscovMatrix = new TMatrixDSym( _n_meas );
    
    // If requested, construct systematic covariance based on specified sources
    if(_set.IncludeNP) {
        // Loop over all inputs
        for(int i = 0; i < _input->GetInput().size(); i++) {
            input in = _input->GetInput()[i];

        // Loop over all Vub BF uncertainty sources
            for(int j = 0; j < in.BF_Vub_Uncertainties.size(); j++) {
                // Does this covariance exist already? If yes add linearly to get the proper error (since both sources are correlated)
                
                if(_cov_mult[in.BF_Vub_Uncertainties[j]]) {
                    // _cov_mult[in.FidAcc_Uncertainties[j]]->Print("all");
                    TMatrixDSym SourceSyscovMatrix(_n_meas);
                    // Loop over all bins
                    for(int k = 0; k < in.BF_Vub_Uncert[j].size(); k++) {
                        // Uncorrelated Errors for ith bin, note that in.FidAcc is in this formula,
                        // since in.FidAcc_Uncert is the relative uncertainty:
                        // f = N * F * C / L then sigma_f = N / F * C / L / (DF / F)
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k = ( N_k * in.BF_Vub_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        // Add both uncertainties and correlate them: sqrt( sigma_f'^2 + sigma_f^2 + 2 * sigma_f * sigma_f' )
                        double sigma_kp = sqrt( pow(sigma_k,2.) + (*_cov_mult[in.BF_Vub_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i]) + 2*sigma_k*sqrt((*_cov_mult[in.BF_Vub_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i])) );
                        // Now add the total uncertainty back into the covariance matrix
                        SourceSyscovMatrix(k+_start_bin[i],k+_start_bin[i]) = sigma_kp * sigma_kp;
                    }
                    // Merge the old and new covariance
                    (*_cov_mult[in.BF_Vub_Uncertainties[j]]) = MergeCovariance(SourceSyscovMatrix,(*_cov_mult[in.BF_Vub_Uncertainties[j]]));
                } else { // This is a new source:
                    TMatrixDSym *SourceSyscovMatrix = new TMatrixDSym( _n_meas );
                    // Loop over all bins
                    for(int k = 0; k < in.BF_Vub_Uncert[j].size(); k++) {
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k =   ( N_k * in.BF_Vub_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        (*SourceSyscovMatrix)(k+_start_bin[i],k+_start_bin[i]) = sigma_k * sigma_k;
                    }
                    
                    _cov_mult[in.BF_Vub_Uncertainties[j]] = SourceSyscovMatrix;
                    _sv_source_mul.push_back( in.BF_Vub_Uncertainties[j] );
                }
            }
            cout <<"Loop over all Vub BF uncertainty sources done"<<endl;
                    // Loop over all Vcb BF uncertainty sources
            for(int j = 0; j < in.BF_Vcb_Uncertainties.size(); j++) {
                // Does this covariance exist already? If yes add linearly to get the proper error (since both sources are correlated)
                
                if(_cov_mult[in.BF_Vcb_Uncertainties[j]]) {
                    // _cov_mult[in.FidAcc_Uncertainties[j]]->Print("all");
                    TMatrixDSym SourceSyscovMatrix(_n_meas);
                    // Loop over all bins
                    for(int k = 0; k < in.BF_Vcb_Uncert[j].size(); k++) {
                        // Uncorrelated Errors for ith bin, note that in.FidAcc is in this formula,
                        // since in.FidAcc_Uncert is the relative uncertainty:
                        // f = N * F * C / L then sigma_f = N / F * C / L / (DF / F)
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k = ( N_k * in.BF_Vcb_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        // Add both uncertainties and correlate them: sqrt( sigma_f'^2 + sigma_f^2 + 2 * sigma_f * sigma_f' )
                        double sigma_kp = sqrt( pow(sigma_k,2.) + (*_cov_mult[in.BF_Vcb_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i]) + 2*sigma_k*sqrt((*_cov_mult[in.BF_Vcb_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i])) );
                        // Now add the total uncertainty back into the covariance matrix
                        SourceSyscovMatrix(k+_start_bin[i],k+_start_bin[i]) = sigma_kp * sigma_kp;
                    }
                    // Merge the old and new covariance
                    (*_cov_mult[in.BF_Vcb_Uncertainties[j]]) = MergeCovariance(SourceSyscovMatrix,(*_cov_mult[in.BF_Vcb_Uncertainties[j]]));
                } else { // This is a new source:
                    TMatrixDSym *SourceSyscovMatrix = new TMatrixDSym( _n_meas );
                    // Loop over all bins
                    for(int k = 0; k < in.BF_Vcb_Uncert[j].size(); k++) {
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k =   ( N_k * in.BF_Vcb_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        (*SourceSyscovMatrix)(k+_start_bin[i],k+_start_bin[i]) = sigma_k * sigma_k;
                    }
                    
                    _cov_mult[in.BF_Vcb_Uncertainties[j]] = SourceSyscovMatrix;
                    _sv_source_mul.push_back( in.BF_Vcb_Uncertainties[j] );
                }
            }
            cout <<"Loop over all Vcb BF uncertainty sources done"<<endl;

        // Loop over all Vub FF uncertainty sources
            for(int j = 0; j < in.FF_Vub_Uncertainties.size(); j++) {
                // Does this covariance exist already? If yes add linearly to get the proper error (since both sources are correlated)
                
                if(_cov_mult[in.FF_Vub_Uncertainties[j]]) {
                    // _cov_mult[in.FidAcc_Uncertainties[j]]->Print("all");
                    TMatrixDSym SourceSyscovMatrix(_n_meas);
                    // Loop over all bins
                    for(int k = 0; k < in.FF_Vub_Uncert[j].size(); k++) {
                        // Uncorrelated Errors for ith bin, note that in.FidAcc is in this formula,
                        // since in.FidAcc_Uncert is the relative uncertainty:
                        // f = N * F * C / L then sigma_f = N / F * C / L / (DF / F)
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k = ( N_k * in.FF_Vub_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        // Add both uncertainties and correlate them: sqrt( sigma_f'^2 + sigma_f^2 + 2 * sigma_f * sigma_f' )
                        double sigma_kp = sqrt( pow(sigma_k,2.) + (*_cov_mult[in.FF_Vub_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i]) + 2*sigma_k*sqrt((*_cov_mult[in.FF_Vub_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i])) );
                        // Now add the total uncertainty back into the covariance matrix
                        SourceSyscovMatrix(k+_start_bin[i],k+_start_bin[i]) = sigma_kp * sigma_kp;
                    }
                    // Merge the old and new covariance
                    (*_cov_mult[in.FF_Vub_Uncertainties[j]]) = MergeCovariance(SourceSyscovMatrix,(*_cov_mult[in.FF_Vub_Uncertainties[j]]));
                } else { // This is a new source:
                    TMatrixDSym *SourceSyscovMatrix = new TMatrixDSym( _n_meas );
                    // Loop over all bins
                    for(int k = 0; k < in.BF_Vub_Uncert[j].size(); k++) {
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k =   ( N_k * in.FF_Vub_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        (*SourceSyscovMatrix)(k+_start_bin[i],k+_start_bin[i]) = sigma_k * sigma_k;
                    }
                    
                    _cov_mult[in.FF_Vub_Uncertainties[j]] = SourceSyscovMatrix;
                    _sv_source_mul.push_back( in.FF_Vub_Uncertainties[j] );
                }
            }
            cout <<"Loop over all Vub FF uncertainty sources done"<<endl;
                    // Loop over all Vcb FF uncertainty sources
            for(int j = 0; j < in.FF_Vcb_Uncertainties.size(); j++) {
                // Does this covariance exist already? If yes add linearly to get the proper error (since both sources are correlated)
                
                if(_cov_mult[in.FF_Vcb_Uncertainties[j]]) {
                    // _cov_mult[in.FidAcc_Uncertainties[j]]->Print("all");
                    TMatrixDSym SourceSyscovMatrix(_n_meas);
                    // Loop over all bins
                    for(int k = 0; k < in.FF_Vcb_Uncert[j].size(); k++) {
                        // Uncorrelated Errors for ith bin, note that in.FidAcc is in this formula,
                        // since in.FidAcc_Uncert is the relative uncertainty:
                        // f = N * F * C / L then sigma_f = N / F * C / L / (DF / F)
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k = ( N_k * in.FF_Vcb_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        // Add both uncertainties and correlate them: sqrt( sigma_f'^2 + sigma_f^2 + 2 * sigma_f * sigma_f' )
                        double sigma_kp = sqrt( pow(sigma_k,2.) + (*_cov_mult[in.FF_Vcb_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i]) + 2*sigma_k*sqrt((*_cov_mult[in.FF_Vcb_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i])) );
                        // Now add the total uncertainty back into the covariance matrix
                        SourceSyscovMatrix(k+_start_bin[i],k+_start_bin[i]) = sigma_kp * sigma_kp;
                    }
                    // Merge the old and new covariance
                    (*_cov_mult[in.FF_Vcb_Uncertainties[j]]) = MergeCovariance(SourceSyscovMatrix,(*_cov_mult[in.FF_Vcb_Uncertainties[j]]));
                } else { // This is a new source:
                    TMatrixDSym *SourceSyscovMatrix = new TMatrixDSym( _n_meas );
                    // Loop over all bins
                    for(int k = 0; k < in.FF_Vcb_Uncert[j].size(); k++) {
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k =   ( N_k * in.FF_Vcb_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        (*SourceSyscovMatrix)(k+_start_bin[i],k+_start_bin[i]) = sigma_k * sigma_k;
                    }
                    
                    _cov_mult[in.FF_Vcb_Uncertainties[j]] = SourceSyscovMatrix;
                    _sv_source_mul.push_back( in.FF_Vcb_Uncertainties[j] );
                }
            }
                cout <<"Loop over all Vcb FF uncertainty sources done"<<endl;
                    // Loop over all general uncertainty sources
            for(int j = 0; j < in.general_Uncertainties.size(); j++) {
                // Does this covariance exist already? If yes add linearly to get the proper error (since both sources are correlated)
                
                if(_cov_mult[in.general_Uncertainties[j]]) {
                    // _cov_mult[in.FidAcc_Uncertainties[j]]->Print("all");
                    TMatrixDSym SourceSyscovMatrix(_n_meas);
                    // Loop over all bins
                    for(int k = 0; k < in.general_Uncert[j].size(); k++) {
                        // Uncorrelated Errors for ith bin, note that in.FidAcc is in this formula,
                        // since in.FidAcc_Uncert is the relative uncertainty:
                        // f = N * F * C / L then sigma_f = N / F * C / L / (DF / F)
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k = ( N_k * in.general_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        // Add both uncertainties and correlate them: sqrt( sigma_f'^2 + sigma_f^2 + 2 * sigma_f * sigma_f' )
                        double sigma_kp = sqrt( pow(sigma_k,2.) + (*_cov_mult[in.general_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i]) + 2*sigma_k*sqrt((*_cov_mult[in.general_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i])) );
                        // Now add the total uncertainty back into the covariance matrix
                        SourceSyscovMatrix(k+_start_bin[i],k+_start_bin[i]) = sigma_kp * sigma_kp;
                    }
                    // Merge the old and new covariance
                    (*_cov_mult[in.general_Uncertainties[j]]) = MergeCovariance(SourceSyscovMatrix,(*_cov_mult[in.general_Uncertainties[j]]));
                } else { // This is a new source:
                    TMatrixDSym *SourceSyscovMatrix = new TMatrixDSym( _n_meas );
                    // Loop over all bins
                    for(int k = 0; k < in.general_Uncert[j].size(); k++) {
                        double N_k = in.Nsig.size() > 0 ? in.Nsig[k] : in.Ndata[k] - in.Nbkg[k];
                        double sigma_k =   ( N_k * in.general_Uncert[j][k] / 100. / in.BR / in.Luminosity );
                        (*SourceSyscovMatrix)(k+_start_bin[i],k+_start_bin[i]) = sigma_k * sigma_k;
                    }
                    
                    _cov_mult[in.general_Uncertainties[j]] = SourceSyscovMatrix;
                    _sv_source_mul.push_back( in.general_Uncertainties[j] );
                }
            }
                cout <<"Loop over all general uncertainty sources done"<<endl;



            cout << "Built correction factors"<<endl;

            // Loop over all correction factor uncertainty sources
            for(int j = 0; j < in.Bkg_Uncertainties.size(); j++) {   
                // Does this covariance exist already? If yes add linearly to get the proper error (since both sources are correlated)        
                if(_cov_mult[in.Bkg_Uncertainties[j]]) {
                    //_cov_mult[in.Bkg_Uncertainties[j]]->Print("all");    
                    TMatrixDSym SourceSyscovMatrix(_n_meas);
                    // Loop over all bins
                    for(int k = 0; k < in.Bkg_Uncert[j].size(); k++) {
                        // Uncorrelated Errors for ith bin, note that in.FidAcc is in this formula, 
                        // since in.FidAcc_Uncert is the relative uncertainty: 
                        // f = (N - NBkg) * F * C / L then sigma_bkg = NBkg / F * C / L * (DNBkg / NBkg)
                        double sigma_k = ( in.Bkg_Uncert[j][k] / in.BR / in.Luminosity );
                        // Add both uncertainties and correlate them: sqrt( sigma_f'^2 + sigma_f^2 - 2 * sigma_f * sigma_f' )
                        // Notice, the minus sign comes from (N-NBkg)
                        double sigma_kp = sqrt( pow(sigma_k,2.) + (*_cov_mult[in.Bkg_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i]) - 2*sigma_k*sqrt((*_cov_mult[in.Bkg_Uncertainties[j]])(k+_start_bin[i],k+_start_bin[i])) );
                        // Now add the total uncertainty back into the covariance matrix
                        SourceSyscovMatrix(k+_start_bin[i],k+_start_bin[i]) = sigma_kp * sigma_kp;                      
                    }
                    // Merge the old and new covariance          
                    (*_cov_mult[in.Bkg_Uncertainties[j]]) = MergeCovariance(SourceSyscovMatrix,(*_cov_mult[in.Bkg_Uncertainties[j]]));
                } else { // This is a new source:
                    TMatrixDSym *SourceSyscovMatrix = new TMatrixDSym( _n_meas );
                    // Loop over all bins
                    for(int k = 0; k < in.Bkg_Uncert[j].size(); k++) {
                        double sigma_k = in.Bkg_Uncert[j][k] / in.BR / in.Luminosity ;
                        (*SourceSyscovMatrix)(k+_start_bin[i],k+_start_bin[i]) = sigma_k*sigma_k;
                    } 
                    _cov_mult[in.Bkg_Uncertainties[j]] = SourceSyscovMatrix;
                    _sv_source_mul.push_back( in.Bkg_Uncertainties[j] );
                }
            }     
            
        }
        
    }  
    
    for(int i = 0; i < _sv_source_mul.size(); i++) { 
        // Now calculate proper bin-by-bin and measurement-by-measurement cross correlation and build the the systematic covariance  
        for(int ni = 0; ni < _cov_mult[_sv_source_mul[i]]->GetNcols(); ni++)
            for(int nj = 0; nj < _cov_mult[_sv_source_mul[i]]->GetNrows(); nj++)
                (*_cov_mult[_sv_source_mul[i]])(ni,nj) = sqrt((*_cov_mult[_sv_source_mul[i]])(ni,ni)) * sqrt((*_cov_mult[_sv_source_mul[i]])(nj,nj));
        (*_SyscovMatrix) += (*_cov_mult[_sv_source_mul[i]]);
    }
    cout << "SyscovMatrix built"<<endl;

    _TotcovMatrix = new TMatrixDSym( _n_meas );
    // Construct total covariance and we are done
    (*_TotcovMatrix) = (*_SyscovMatrix) + (*_StatcovMatrix);

    Printf(" Systematic Covariance: ");
    _SyscovMatrix->Print("all");
    
    Printf(" Statistical Covariance: ");
    _StatcovMatrix->Print("all");
    
    Printf(" Total Covariance: ");
    _TotcovMatrix->Print("all");
    
    // Construct correlation matrices 
    _StatcorMatrix = new TMatrixDSym( _n_meas );
    _SyscorMatrix = new TMatrixDSym( _n_meas );
    _TotcorMatrix = new TMatrixDSym( _n_meas );
    
    for(int i = 0; i < _StatcovMatrix->GetNcols(); i++)
        for(int j = 0; j < _StatcovMatrix->GetNrows(); j++)
            (*_StatcorMatrix)(i,j) = 100*(*_StatcovMatrix)(i,j)/sqrt((*_StatcovMatrix)(i,i))/sqrt((*_StatcovMatrix)(j,j));
    
    for(int i = 0; i < _SyscorMatrix->GetNcols(); i++)
        for(int j = 0; j < _SyscorMatrix->GetNrows(); j++)
            (*_SyscorMatrix)(i,j) = 100*(*_SyscovMatrix)(i,j)/sqrt((*_SyscovMatrix)(i,i))/sqrt((*_SyscovMatrix)(j,j));
    
    for(int i = 0; i < _TotcovMatrix->GetNcols(); i++)
        for(int j = 0; j < _TotcovMatrix->GetNrows(); j++)
            (*_TotcorMatrix)(i,j) = 100*(*_TotcovMatrix)(i,j)/sqrt((*_TotcovMatrix)(i,i))/sqrt((*_TotcovMatrix)(j,j));
    
    
}


