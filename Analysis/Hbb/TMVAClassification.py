#!/usr/bin/env python
# @(#)root/tmva $Id: TMVAClassification.py,v 1.59.2.1 2008/12/02 08:51:43 andreas.hoecker Exp $
# ------------------------------------------------------------------------------ #
# Project      : TMVA - a Root-integrated toolkit for multivariate data analysis #
# Package      : TMVA                                                            #
# Python script: TMVAClassification.py                                           #
#                                                                                #
# This python script provides examples for the training and testing of all the   #
# TMVA classifiers through PyROOT. Note that the use PyROOT requires that you    #
# have a python version > 2.2 installed on your computer.                        #
#                                                                                #
# The Application works similarly, please see:                                   #
#    TMVA/macros/TMVAClassificationApplication.C                                 #
# For regression, see:                                                           #
#    TMVA/macros/TMVARegression.C                                                #
#    TMVA/macros/TMVARegressionpplication.C                                      #
# and translate to python as done here.                                          #
#                                                                                #
# As input data is used a toy-MC sample consisting of four Gaussian-distributed  #
# and linearly correlated input variables.                                       #
#                                                                                #
# The methods to be used can be switched on and off via the prompt command, for  #
# example:                                                                       #
#                                                                                #
#    python TMVAClassification.py --methods Fisher,Likelihood                    #
#                                                                                #
# The output file "TMVA.root" can be analysed with the use of dedicated          #
# macros (simply say: root -l <../macros/macro.C>), which can be conveniently    #
# invoked through a GUI that will appear at the end of the run of this macro.    #
#                                                                                #
# for help type "python TMVAClassification.py --help"                            #
# ------------------------------------------------------------------------------ #

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os
import math
# --------------------------------------------

# Default settings for command line arguments
OUTFNAME = "TMVA.root"
OUTDIR = "TMVAout"
ANALYSIS="Dijet"
TREE  = "tree"

INPUTDIR="dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2"
SAMPLES={
    # key      :( Sig or Bckg,    rootfile,                                               cross section)
    "HIGGS"    :("S","DiJetPt_ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp.root",                   22.97),
    "ZJets_ptL":("B","DiJetPt_DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball.root",   62140.0),
    "ZJets_ptH":("B","DiJetPt_DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph.root",               40510.0),
    "ZZ"       :("B","DiJetPt_ZZ_TuneZ2star_8TeV_pythia6_tauola.root",                         8255.61),
    "WW"       :("B","DiJetPt_WW_TuneZ2star_8TeV_pythia6_tauola.root",                         57109.7),
    "WZ"       :("B","DiJetPt_WZ_TuneZ2star_8TeV_pythia6_tauola.root",                         32316.1),
    "TTbar"    :("B","DiJetPt_TTJets_TuneZ2star_8TeV-madgraph-tauola.root",                    225197.),
    "ST-s"     :("B","DiJetPt_T_s-channel_TuneZ2star_8TeV-powheg-tauola.root",                 3893.94),
    "STbar-s"  :("B","DiJetPt_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola.root",              1757.76),
    ## ## "ST-t"           :("hsts/Hbb_Zee__T_t-channel.root"                  ,xx,           30004.2),
    "STbar-t"  :("B","DiJetPt_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola.root",              55531.0),
    "ST-tW"    :("B","DiJetPt_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root",             11177.3),
    "STbar-tW" :("B","DiJetPt_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root" ,         11177.3),
    }

FRIENDS={
    # key      :( Sig or Bckg,    rootfile,                                               cross section)
    "HIGGS"    :("friend.root"),
    "ZJets_ptL":("friendDY.root"),
    }

# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % OUTFNAME
    print "  -d | --outputdir : name of output directory for the xml files (default: '%s')" % OUTDIR
    print "  -a | --analysis : Dijet or Subjet (default: '%s')" % ANALYSIS
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

def deltaPhi(phi1,phi2):
                      
    result=phi1-phi2
    while (result > math.pi): result -= 2*math.pi
    while (result <= -math.pi): result += 2*math.pi

    return result

# Main routine
def main():

    try:
        # retrive command line options
        shortopts  = "a:o:d:vh?"
        longopts   = ["analysis=","outputfile=", "outputdir=","verbose", "help", "usage"]
        opts, args = getopt.getopt( sys.argv[1:], shortopts, longopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    _outfname   = OUTFNAME
    _outdir     = OUTDIR
    _analysis   = ANALYSIS
    verbose     = False
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-o", "--outputfile"):
            _outfname = a
        elif o in ("-d", "--outputdir"):
            _outdir = a
        elif o in ("-a", "--analysis"):
            _analysis = a
        elif o in ("-v", "--verbose"):
            verbose = True

    # Import ROOT classes
    from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut
    
    # check ROOT version, give alarm if 5.18 
    if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
        print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
        print "*** does not run properly (function calls with enums in the argument are ignored)."
        print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
        print "*** or use another ROOT version (e.g., ROOT 5.19)."
        sys.exit(1)
        
    # Import TMVA classes from ROOT
    from ROOT import TMVA


    # Output file
    outputFile = TFile( _outfname, 'RECREATE' )
    
    # Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for more factory options)
    # All TMVA output can be suppressed by removing the "!" (not) in 
    # front of the "Silent" argument in the option string
    factory = TMVA.Factory( "TMVAClassification", outputFile, 
                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose( verbose )
    
    # If you wish to modify default settings 
    # (please check "src/Config.h" to see all available global options)
    #    gConfig().GetVariablePlotting()).fTimesRMS = 8.0
    #    gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory"
    # TMVA.gConfig().GetIONames().fWeightFileDir = "weights_" + _analysis
    # OutDir="weights_" + _analysis

    if os.path.exists(_outdir):
        print "Output directory: ",_outdir," already exists, please specify new directory"
        sys.exit(1)

    os.makedirs(_outdir)
    TMVA.gConfig().GetIONames().fWeightFileDir = _outdir
        

    # Define the input variables that shall be used for the classifier training
    # note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
    # [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
    if _analysis == "Dijet":



        factory.AddVariable("H_mass := H.mass", 'F');
        factory.AddVariable("H_pt :=H.pt", 'F');
        factory.AddVariable("V_pt :=V.pt", 'F');
        factory.AddVariable("hJ12_MaxCsv := max(hJet_csv[0],hJet_csv[1])", 'F');
        factory.AddVariable("hJ12_MinCsv := min(hJet_csv[0],hJet_csv[1])", 'F');
        factory.AddVariable("HV_dPhi := HVdPhi", 'F');
        factory.AddVariable("H_dEta := H.dEta", 'F');

        # You can add so-called "Spectator variables", which are not used in the MVA training, 
        # but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
        # input variables, the response values of all trained MVAs, and the spectator variables
        factory.AddSpectator("hJet_pt1 := hJet_pt[0]", 'F')
        factory.AddSpectator("hJet_pt2 := hJet_pt[1]", 'F')

    elif _analysis == "Subjet":
        factory.AddVariable("H_mass := FatH.filteredmass", 'F');
        factory.AddVariable("H_pt   := FatH.filteredpt", 'F');
        factory.AddVariable("V_pt   := V.pt", 'F');
        # factory.AddVariable("HV_dPhi := math.fabs(deltaPhi(FatH.filteredphi,V.phi))", 'F');
        factory.AddVariable("HV_dPhi := " +\
                             "FatH.filteredphi - V.phi > pi ? " +\
                             "abs(FatH.filteredphi - V.phi - 2*pi) : " +\
                             "FatH.filteredphi - V.phi < -pi ? " +\
                             "abs(FatH.filteredphi - V.phi + 2*pi) : " +\
                             "abs(FatH.filteredphi - V.phi)", 'F' )

        factory.AddVariable("SJ1_csv := fathFilterJets_csv[0]", 'F');
        factory.AddVariable("SJ2_csv := fathFilterJets_csv[1]", 'F');
        factory.AddVariable("SJ3_csv := Alt$(fathFilterJets_csv[2],0)", 'F');

        factory.AddSpectator("SJ1_pt := fathFilterJets_pt[0]", 'F');
        factory.AddSpectator("SJ2_pt := fathFilterJets_pt[1]", 'F');
        factory.AddSpectator("SJ3_pt := Alt$(fathFilterJets_pt[2],0)", 'F');

    else:
        print "Problem specifying analysis. Please choose Dijet or Subjet."
        sys.exit(1)


    ## Get the Signal and Background trees
    for Sample in SAMPLES.keys():
        SampleInfo=SAMPLES[Sample]
        ## FriendInfo=FRIENDS[Sample]

        SampleType=SampleInfo[0] # signal or background
        infile=os.path.join(INPUTDIR,SampleInfo[1])
        xs=SampleInfo[2]


        f=TFile.Open(infile)# Open the input file 
        h = f.Get("Count")  # get number of step 1 events
        nEVT=int(h.GetBinContent(1))

        wt  =xs/(nEVT)        
        print Sample,": ",infile
        print "XS:nEVT:wt: ", xs,nEVT,wt

        theTree      = f.Get( TREE )


        # theTree.AddFriend("treeFR",FriendInfo)
        if SampleType == "S":
            factory.AddSignalTree    ( theTree, wt )
        elif SampleType == "B":
            factory.AddBackgroundTree( theTree, wt )        
        else:
            print "Trouble extracting SampleType for this sample"
            sys.exit(1)
        

    # table10 AN-2011/430
    if _analysis == "Dijet":
        cutString=\
            "Vtype == 1"             + " && " +\
            "H.HiggsFlag > 0"        + " && " +\
            "V.mass > 75.0"          + " && " +\
            "V.mass < 105.0"         + " && " +\
            "V.pt > 100.0"           + " && " +\
            "hJet_pt[0] > 20.0"      + " && " +\
            "hJet_pt[1] > 20.0"      + " && " +\
            "H.mass > 80.0"          + " && " +\
            "H.mass < 150.0"         + " && " +\
            "max(hJet_csv[0],hJet_csv[1]) > 0.244"  + " && " +\
            "min(hJet_csv[0],hJet_csv[1]) > 0.244" 
    elif _analysis == "Subjet":
        cutString=\
            "Vtype == 1"             + " && " +\
            "FatH.FatHiggsFlag > 0"  + " && " +\
            "V.mass > 75.0"          + " && " +\
            "V.mass < 105.0"         + " && " +\
            "V.pt > 100.0"           + " && " +\
            "nfathFilterJets >= 2"   + " && " +\
            "fathFilterJets_pt[0] > 20.0"      + " && " +\
            "fathFilterJets_pt[1] > 20.0"      + " && " +\
            "FatH.filteredmass > 80.0"          + " && " +\
            "FatH.filteredmass < 150.0"         + " && " +\
            "max(fathFilterJets_csv[0],fathFilterJets_csv[1]) > 0.244"  + " && " +\
            "min(fathFilterJets_csv[0],fathFilterJets_csv[1]) > 0.244" 
    else:
        print "Problem specifying analysis. Please choose Dijet or Subjet."
        sys.exit(1)


    print cutString
    mycutSig = TCut( cutString ) 
    mycutBkg = TCut( cutString ) 
    
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples

    # prepareOptions="nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V"
    prepareOptions="SplitMode=Random:!V"
    factory.PrepareTrainingAndTestTree( mycutSig, mycutBkg, prepareOptions)


    bdtOptions = \
        "!H"                       + ":" +\
        "!V"                       + ":" +\
        "NTrees=850"               + ":" +\
        "nEventsMin=150"           + ":" +\
        "MaxDepth=3"               + ":" +\
        "BoostType=AdaBoost"       + ":" +\
        "AdaBoostBeta=0.3"         + ":" +\
        "SeparationType=GiniIndex" + ":" +\
        "nCuts=20"                 + ":" +\
        "PruneMethod=NoPruning"
    # "PruneMethod=CostComplexity"

    # 

    print bdtOptions

    factory.BookMethod( TMVA.Types.kBDT, "BDT", bdtOptions)
   

    # Train MVAs
    factory.TrainAllMethods()
    
    # Test MVAs
    factory.TestAllMethods()
    
    # Evaluate MVAs
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()
    
    print "=== wrote root file %s\n" % _outfname
    print "=== TMVAClassification is done!\n"
    
    # open the GUI for the result macros    
    #gROOT.ProcessLine( "TMVAGui(\"%s\")" % _outfname )
    
    # keep the ROOT thread running
    #gApplication.Run() 

# ----------------------------------------------------------

if __name__ == "__main__":
    main()
