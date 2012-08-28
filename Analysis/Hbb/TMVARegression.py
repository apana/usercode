#// @(#)root/tmva $Id: TMVARegression.C 38475 2011-03-17 10:46:00Z evt $
#/**********************************************************************************
# * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
# * Package   : TMVA                                                               *
# * Root Macro: TMVARegression                                                     *
# *                                                                                *
# * This macro provides examples for the training and testing of the               *
# * TMVA classifiers.                                                              *
# *                                                                                *
# * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
# * and linearly correlated input variables.                                       *
# *                                                                                *
# * The methods to be used can be switched on and off by means of booleans, or     *
# * via the prompt command, for example:                                           *
# *                                                                                *
# *    root -l TMVARegression.C\(\"LD,MLP\"\)                                      *
# *                                                                                *
# * (note that the backslashes are mandatory)                                      *
# * If no method given, a default set is used.                                     *
# *                                                                                *
# * The output file "TMVAReg.root" can be analysed with the use of dedicated       *
# * macros (simply say: root -l <macro.C>), which can be conveniently              *
# * invoked through a GUI that will appear at the end of the run of this macro.    *
# **********************************************************************************/

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os
import math
import ROOT

#import TMVA/src/Tools.h
#import TMVA/src/Factory.h
# --------------------------------------------            

# Import ROOT classes
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut, TChain


# functions such as deltaR to be used in BDT
ROOT.gROOT.ProcessLine('.L funcs.C+')
from ROOT import deltaPhi, evalEt, evalMt, resolutionBias, METdeltaPhi



OUTFNAME = "TMVAReg_Dijet.root"
#OUTFNAME = "TMVAReg_Subjet.root"
ANALYSIS = "Dijet"
#ANALYSIS = "Subjet"


# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % OUTFNAME
    print "  -a | --analysis : Dijet or Subjet (default: '%s')" % ANALYSIS
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "


# Main routine
def TMVARegression():

    try:
        # retrieve command line options
        shortopts  = "a:o:vh?"
        longopts   = ["analysis=","outputfile=", "verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )
        
    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)
        
    _outfname   = OUTFNAME
    _analysis   = ANALYSIS
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-o", "--outputfile"):
            _outfname = a
        elif o in ("-a", "--analysis"):
            _analysis = a
        elif o in ("-v", "--verbose"):
            verbose = True

    
            
    # Import TMVA classes from ROOT
    from ROOT import TMVA

    # Output file
    outputFile = TFile( _outfname, 'RECREATE' )

    #   // Create the factory object. Later you can choose the methods
    #   // whose performance you'd like to investigate. The factory will
    #   // then run the performance analysis for you.
    #   //
    #   // The first argument is the base of the name of all the
    #   // weightfiles in the directory weights_Reg/ 
    #   //
    #   // The second argument is the output file for the training results
    #   // All TMVA output can be suppressed by removing the "!" (not) in 
    #   // front of the "Silent" argument in the option string
    factory = TMVA.Factory ("TMVARegression", outputFile, 
                                                   "!V:!Silent:Color:DrawProgressBar" )
    # Set verbosity
    factory.SetVerbose( verbose )
    
    TMVA.gConfig().GetIONames().fWeightFileDir = "weights_Reg_8TeV" + "_" + _analysis

    if _analysis == "Dijet":
  
        factory.AddVariable("hJet_pt", "hJet_pt", "units", 'F')
        factory.AddVariable("hJet_eta", "hJet_eta", "units", 'F')
        factory.AddVariable("hJet_phi", "hJet_phi", "units", 'F')
        factory.AddVariable("hJet_ptRaw*((hJet_ptRaw+resolutionBias(fabs(hJet_eta))*(hJet_ptRaw-hJet_genPt))/hJet_ptRaw)", "hJet_ptRaw*((hJet_ptRaw+resolutionBias(fabs(hJet_eta))*(hJet_ptRaw-hJet_genPt))/hJet_ptRaw)", "units", 'F')
        factory.AddVariable("hJet_Mt:=evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)","hJet_Mt", "units", 'F')
        factory.AddVariable("hJet_Et:=evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)","hJet_Et", "units", 'F')
        factory.AddVariable("hJet_ptLeadTrack", "hJet_ptLeadTrack", "units", 'F')
        factory.AddVariable("hJet_vtxPt", "hJet_vtxPt", "units", 'F')
        factory.AddVariable("hJet_vtx3dL", "hJet_vtx3dL", "units", 'F')
        factory.AddVariable("hJet_vtx3deL", "hJet_vtx3deL", "units", 'F')
        factory.AddVariable("hJet_vtxMass", "hJet_vtxMass", "units", 'F')
        factory.AddVariable("hJet_chf", "hJet_chf", "units", 'F')
        factory.AddVariable("hJet_nch", "hJet_nch", "units", 'F')
        factory.AddVariable("hJet_nconstituents", "hJet_nconstituents", "units", 'F')
        factory.AddVariable("hJet_JECUnc", "hJet_JECUnc", "units", 'F')
        factory.AddVariable("rho25", "rho25", "units", 'F')
        factory.AddVariable("MET.et", "MET.et", "units", 'F')
        factory.AddVariable("METdPhi:=METdeltaPhi(MET.phi, hJet_phi[0], hJet_phi[1])","METdPhi", "units",'F')

        #Add the variable carrying the regression target
        factory.AddTarget( "hJet_genPt" )

    elif _analysis == "Subjet":

        factory.AddVariable("fathFilterJets_pt", "fathFilterJets_pt", "units", 'F')
        factory.AddVariable("fathFilterJets_eta", "fathFilterJets_eta", "units", 'F')
        factory.AddVariable("fathFilterJets_phi", "fathFilterJets_phi", "units", 'F')
        factory.AddVariable("fathFilterJets_ptRaw*((fathFilterJets_ptRaw+resolutionBias(fabs(fathFilterJets_eta))*(fathFilterJets_ptRaw-fathFilterJets_genPt))/fathFilterJets_ptRaw)", "fathFilterJets_ptRaw*((fathFilterJets_ptRaw+resolutionBias(fabs(fathFilterJets_eta))*(fathFilterJets_ptRaw-fathFilterJets_genPt))/fathFilterJets_ptRaw)", "units", 'F')
        factory.AddVariable("fathFilterJets_Mt:=evalMt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)","fathFilterJets_Mt", "units", 'F')
        factory.AddVariable("fathFilterJets_Et:=evalEt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)","fathFilterJets_Et", "units", 'F')
        factory.AddVariable("fathFilterJets_ptLeadTrack", "fathFilterJets_ptLeadTrack", "units", 'F')
        factory.AddVariable("fathFilterJets_vtxPt", "fathFilterJets_vtxPt", "units", 'F')
        factory.AddVariable("fathFilterJets_vtx3dL", "fathFilterJets_vtx3dL", "units", 'F')
        factory.AddVariable("fathFilterJets_vtx3deL", "fathFilterJets_vtx3deL", "units", 'F')
        factory.AddVariable("fathFilterJets_vtxMass", "fathFilterJets_vtxMass", "units", 'F')
        factory.AddVariable("fathFilterJets_chf", "fathFilterJets_chf", "units", 'F')
        factory.AddVariable("rho25", "rho25", "units", 'F')
        factory.AddVariable("MET.et", "MET.et", "units", 'F')
        factory.AddVariable("METdPhi:=METdeltaPhi(MET.phi, fathFilterJets_phi[0], fathFilterJets_phi[1])","METdPhi", "units",'F')

        factory.AddTarget("fathFilterJets_genPt")

    else:
        print "Problem specifying analysis. Please choose Dijet or Subjet."
        sys.exit(1) 

    ## Get the Signal trees
    en7TeV = False
    en8TeV = True

    regWeight = 1.
    chain = TChain("tree")

    if en7TeV: #change the ntuple names later!!
        chain.Add("Step2_output_May11/WH_125_ForRegression.root")
        chain.Add("Step2_output_May11/WH_115_ForRegression.root")
        chain.Add("Step2_output_May11/WH_120_ForRegression.root")
        chain.Add("Step2_output_May11/WH_130_ForRegression.root")
        chain.Add("Step2_output_May11/WH_135_ForRegression.root")
        

    if en8TeV and _analysis == "Dijet":
        chain.Add("dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp.root")
        chain.Add("dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp.root")
        chain.Add("dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp.root")
        chain.Add("dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp.root")
        chain.Add("dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp.root")
        

    if en8TeV and _analysis == "Subjet": 
        chain.Add("/uscmst1b_scratch/lpc1/lpctrig/apana/Higgs/Step2/NtupleV34/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Step2/ZH/ZH_110_summer12_33b.root")
        chain.Add("/uscmst1b_scratch/lpc1/lpctrig/apana/Higgs/Step2/NtupleV34/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Step2/ZH/ZH_115_summer12_33b.root")
        chain.Add("/uscmst1b_scratch/lpc1/lpctrig/apana/Higgs/Step2/NtupleV34/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Step2/ZH/ZH_120_summer12_33b.root")
        chain.Add("/uscmst1b_scratch/lpc1/lpctrig/apana/Higgs/Step2/NtupleV34/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Step2/ZH/ZH_125_summer12_33b.root")
        chain.Add("/uscmst1b_scratch/lpc1/lpctrig/apana/Higgs/Step2/NtupleV34/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Step2/ZH/ZH_130_summer12_33b.root")
        chain.Add("/uscmst1b_scratch/lpc1/lpctrig/apana/Higgs/Step2/NtupleV34/CMSSW_5_2_5/src/VHbbAnalysis/VHbbDataFormats/bin/Step2/ZH/ZH_135_summer12_33b.root")

        
    NEntries = chain.GetEntries()
    print "Number of entries on Chain:",NEntries

    regTree = chain
    
    factory.AddRegressionTree( regTree, regWeight )

    #This would set individual event weights (the variables defined in the 
    #expression need to exist in the original TTree)
    #factory->SetWeightExpression( "var1", "Regression" )


    if _analysis == "Dijet":
        cutString=\
            "Vtype == 0"                            + " && " +\
            "hJet_pt[0] > 20.0"                     + " && " +\
            "hJet_pt[1] > 20.0"                     + " && " +\
            "hJet_eta[0] < 2.4"                     + " && " +\
            "hJet_eta[1] < 2.4"                     + " && " +\
            "max(hJet_csv[0],hJet_csv[1]) > 0.0"    + " && " +\
            "min(hJet_csv[0],hJet_csv[1]) > 0.0"    + " && " +\
            "H.pt > 120"


    elif _analysis == "Subjet":
        cutString=\
            "Vtype == 0"                            + " && " +\
            "fathFilterJets_pt[0] > 20.0"                     + " && " +\
            "fathFilterJets_pt[1] > 20.0"                     + " && " +\
            "fathFilterJets_eta[0] < 2.4"                     + " && " +\
            "fathFilterJets_eta[1] < 2.4"                     + " && " +\
            "max(fathFilterJets_csv[0],fathFilterJets_csv[1]) > 0.0"    + " && " +\
            "min(fathFilterJets_csv[0],fathFilterJets_csv[1]) > 0.0"    + " && " +\
            "FatH.filteredpt > 120"

    else:
        print "Problem specifying analysis. Please choose Dijet or Subjet."
        sys.exit(1)

    print cutString
    mycut = TCut( cutString )
        
    
    # tell the factory to use all remaining events in the trees after training for testing. The number is 25% of the events after cuts:
    if en7TeV:
        factory.PrepareTrainingAndTestTree( mycut, "nTrain_Regression=125000:nTest_Regression=125000:SplitMode=Random:NormMode=NumEvents:!V" )
    if en8TeV:
        factory.PrepareTrainingAndTestTree( mycut, "nTrain_Regression=130500:nTest_Regression=130500:SplitMode=Random:NormMode=NumEvents:!V" )

    #If no numbers of events are given, half of the events in the tree are used 
    #for training, and the other half for testing:
    #factory.PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  

    #---- Book MVA methods
   
    #please lookup the various method configuration options in the corresponding cxx files, eg:
    #src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    #it is possible to preset ranges in the option string in which the cut optimisation should be done:
    #"...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable


    #Boosted Decision Trees
    factory.BookMethod( TMVA.Types.kBDT, "BDT",
                        "!H:!V:NTrees=60:nEventsMin=5:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" )

    
    # -------------------------------------------------------------------------------------------

    #---- Now you can tell the factory to train, test, and evaluate the MVAs

    # Train MVAs using the set of training events
    factory.TrainAllMethods()

    # ---- Evaluate all MVAs using the set of test events
    factory.TestAllMethods()

    # ----- Evaluate and compare performance of all configured MVAs
    factory.EvaluateAllMethods()    

    # --------------------------------------------------------------

    
    NEntries = regTree.GetEntries()
    print "Number of entries on Tree: ",NEntries

    # Save the output
    outputFile.Close()

    print "==> Wrote root file %s\n" % _outfname
    print "==> TMVARegression is done!\n"       

    #Launch the GUI for the root macros
    #if (!gROOT->IsBatch()) TMVARegGui( outfileName );


if __name__ == "__main__":
    TMVARegression()
