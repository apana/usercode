#!/usr/bin/env python

import sys,string,math,os
import ConfigParser

from array import array
from optparse import OptionParser

from PhysicsTools.PythonAnalysis import *
from ROOT import *
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()

# import ROOT
# from ROOT import TLorentzVector
# from ROOT import gDirectory, gROOT, gFile, AddressOf


# from ROOT import AutoLibraryLoader


from myRootIOFuncs import *


######################################


Analyses=["Boosted","AK5"]
CandidateDict={"Zmumu":0, "Zee":1, "Wmun":2, "Wen":3, "Znn":4}  # update checkForV is these definitions are changed

## Define mass histogram binning here
nmbins=60
minmass=0.
maxmass=300.

CSVT=0.898

# ROOT.gROOT.LoadMacro('higgsStruct.h+')
ROOT.gROOT.LoadMacro('higgsStruct_new.h+')

########################################

def defineArray(type,size):
    obj = array(type)
    for f in range(0,size):
        obj.append(0)
    return obj

def usage():
    """ Usage: HbbAnalyzer [options] cfgFile outputFile inputFiles ( --hh for options)

    """
    pass

class BDT_READER_AK5:
    def __init__(self,weightFile):
        self.reader=TMVA.Reader("Color:!Silent");

        self.bdt_Hmass = array( 'f', [ 0 ] )
        self.bdt_Hpt = array( 'f', [ 0 ] )
        self.bdt_Vpt = array( 'f', [ 0 ] )
        self.bdt_csv1 = array( 'f', [ 0 ] )
        self.bdt_csv2 = array( 'f', [ 0 ] )
        self.bdt_HVdPhi = array( 'f', [ 0 ] )
        self.bdt_dEta = array( 'f', [ 0 ] )

        self.bdt_hpt1 = array( 'f', [ 0 ] )
        self.bdt_hpt2 = array( 'f', [ 0 ] )

        self.reader.AddVariable("H_mass := H.mass", self.bdt_Hmass );
        self.reader.AddVariable("H_pt := H.pt", self.bdt_Hpt );
        self.reader.AddVariable("V_pt :=V.pt", self.bdt_Vpt );
        self.reader.AddVariable("hJ12_MaxCsv := max(hJet_csv[0],hJet_csv[1])", self.bdt_csv1 );
        self.reader.AddVariable("hJ12_MinCsv := min(hJet_csv[0],hJet_csv[1])", self.bdt_csv2 );
        self.reader.AddVariable("HV_dPhi := HVdPhi", self.bdt_HVdPhi );
        self.reader.AddVariable("H_dEta := H.dEta", self.bdt_dEta );

        self.reader.AddSpectator("hJet_pt[0]", self.bdt_hpt1 );
        self.reader.AddSpectator("hJet_pt[1]", self.bdt_hpt2 );

        self.reader.BookMVA("BDT", weightFile);

    def loadVar(self):

        # load the BDT variables with corresponding values from the Step2 ntuple
        self.bdt_Hmass[0]=H.mass
        self.bdt_Hpt[0]=H.pt
        self.bdt_Vpt[0]=V.pt
        self.bdt_csv1[0]=tree.hJet_csv[0]
        self.bdt_csv2[0]=tree.hJet_csv[1]
        self.bdt_HVdPhi[0]=tree.HVdPhi
        self.bdt_dEta[0]=tree.dEta

    def evaluate(self):
        bdt=self.reader.EvaluateMVA("BDT")
        return bdt


class BDT_READER_BST:
    def __init__(self,weightFile):
        self.reader=TMVA.Reader("Color:!Silent");

        self.bdt_Hmass = array( 'f', [ 0 ] )
        self.bdt_Hpt = array( 'f', [ 0 ] )
        self.bdt_Vpt = array( 'f', [ 0 ] )
        self.bdt_csv1 = array( 'f', [ 0 ] )
        self.bdt_csv2 = array( 'f', [ 0 ] )
        self.bdt_csv3 = array( 'f', [ 0 ] )
        self.bdt_HVdPhi = array( 'f', [ 0 ] )
        self.bdt_dEta = array( 'f', [ 0 ] )

        self.bdt_fjpt1 = array( 'f', [ 0 ] )
        self.bdt_fjpt2 = array( 'f', [ 0 ] )
        self.bdt_fjpt3 = array( 'f', [ 0 ] )

        self.reader.AddVariable("H_mass := FatH.filteredmass", self.bdt_Hmass )
        self.reader.AddVariable("H_pt := FatH.filteredpt", self.bdt_Hpt )
        self.reader.AddVariable("V_pt :=V.pt", self.bdt_Vpt )
        self.reader.AddVariable("HV_dPhi := " +\
                                    "FatH.filteredphi - V.phi > pi ? " +\
                                    "abs(FatH.filteredphi - V.phi - 2*pi) : " +\
                                    "FatH.filteredphi - V.phi < -pi ? " +\
                                    "abs(FatH.filteredphi - V.phi + 2*pi) : " +\
                                    "abs(FatH.filteredphi - V.phi)", self.bdt_HVdPhi)

        self.reader.AddVariable("SJ1_csv := fathFilterJets_csv[0]", self.bdt_csv1);
        self.reader.AddVariable("SJ2_csv := fathFilterJets_csv[1]", self.bdt_csv2);
        self.reader.AddVariable("SJ3_csv := Alt$(fathFilterJets_csv[2],0)", self.bdt_csv3);

        self.reader.AddSpectator("SJ1_pt := fathFilterJets_pt[0]", self.bdt_fjpt1);
        self.reader.AddSpectator("SJ2_pt := fathFilterJets_pt[1]", self.bdt_fjpt2);
        self.reader.AddSpectator("SJ3_pt := Alt$(fathFilterJets_pt[2],0)", self.bdt_fjpt3);

        self.reader.BookMVA("BDT", weightFile);

    def loadVar(self):

        # load the BDT variables with corresponding values from the Step2 ntuple
        self.bdt_Hmass[0]=FatH.filteredmass
        self.bdt_Hpt[0]=FatH.filteredpt
        self.bdt_Vpt[0]=V.pt
        self.bdt_HVdPhi[0]=math.fabs(deltaPhi(FatH.filteredphi,V.phi))

        # print len(tree.fathFilterJets_csv), FatH.FatHiggsFlag
        
        self.bdt_csv1[0]=-999.
        self.bdt_csv2[0]=-999.
        self.bdt_csv3[0]=-999.

        self.bdt_fjpt1[0]=-999.
        self.bdt_fjpt2[0]=-999.
        self.bdt_fjpt3[0]=-999.

        if tree.nfathFilterJets>0:
            self.bdt_csv1[0]=tree.fathFilterJets_csv[0]
            self.bdt_fjpt1[0]=tree.fathFilterJets_pt[0]

        if tree.nfathFilterJets>1:
            self.bdt_csv2[0]=tree.fathFilterJets_csv[1]
            self.bdt_fjpt2[0]=tree.fathFilterJets_pt[1]

        if tree.nfathFilterJets>2:
            self.bdt_csv3[0]=tree.fathFilterJets_csv[2]
            self.bdt_fjpt3[0]=tree.fathFilterJets_pt[2]

    def evaluate(self):
        bdt=self.reader.EvaluateMVA("BDT")
        return bdt

class JobPar:
    def __init__(self):
        getJobopt(self)

def getJobopt(self):


    parser = OptionParser(usage=usage.__doc__)
    parser.add_option("-n", "--nevts", dest="nevts", default=-1, type="int",
                      help="number of events to process (default=-1 (all))")
    parser.add_option("-d", "--data", action="store_true", dest="isData",default=False,
                      help="default=False")
    parser.add_option("--lj", action="store_true", dest="lightJets",default=False,
                      help="Cut on light (udscg) flavour jets -- default=False")
    parser.add_option("--hj", action="store_true", dest="heavyJets",default=False,
                      help="Cut on heavy (b) jets -- default=False")
    parser.add_option("--mc7TeV", action="store_true", dest="MC7TeV",default=False,
                      help="default=False")
    parser.add_option("--hh", action="store_true", dest="help",default=False,
                      help="print help")

    ## parser.add_option("-c", "--cfg", dest="cfgfile",
    ##                   help="job configuration file", metavar="FILE")
    ## parser.add_option("-i", "--infile", dest="inputFile",type="string",
    ##                   help="input files to process", metavar="INFILE")
    ## parser.add_option("-o", "--outfile", dest="outputFile",type="string",
    ##                   help="output histogram file", metavar="OUTFILE")
    ## parser.add_option("-q", "--quiet",
    ##                   action="store_false", dest="verbose", default=True,
    ##                   help="don't print status messages to stdout")

    (options, args) = parser.parse_args()

    narg=len(args)

    if options.help:
        parser.print_help()
        sys.exit(0)

    if narg < 3 :
        print usage.__doc__
        parser.print_help()
        sys.exit(1)


    cfgFile            = args[0]
    self.outFile       = args[1]
    self.inputFiles    = args[2]

    print self.inputFiles
    
    self.nevts         = options.nevts
    self.isData        = options.isData

    self.lightJets     = False
    self.heavyJets     = False

    self.MC7TeV = False

    if not self.isData:
        self.lightJets     = options.lightJets
        self.heavyJets     = options.heavyJets

        self.MC7TeV        = options.MC7TeV

    if self.lightJets and self.heavyJets:
        print "Error: Job setup to cut on both heavy and light flavour quarks"
        sys.exit(1)


    if not os.path.exists(cfgFile):
        print "Configuration file ",cfgFile," does not exist -- Exiting program"
        sys.exit(1)

    Cfg = ConfigParser.ConfigParser()
    Cfg.read(cfgFile)


    # self.nevts         = Cfg.getint("InputOutput", "nevts")
    # self.outFile       = Cfg.get("InputOutput", "outFile")
    # self.inputFiles    = Cfg.get("InputOutput", "inputFiles")
    # self.isData        = Cfg.getboolean("InputOutput", "isData")

    self.WhichAnalysis = Cfg.get("InputOutput", "WhichAnalysis")    
    # self.MC7TeV        = Cfg.getboolean("InputOutput", "MC7TeV")

    self.ApplyJECResidual = Cfg.getboolean("SubjetAnalysis", "ApplyJECResidual")

    self.MinLLpT_BST   = Cfg.getfloat("SubjetAnalysis", "MinLLpT_BST")
    self.MaxLLpT_BST   = Cfg.getfloat("SubjetAnalysis", "MaxLLpT_BST")

    self.MinFatJetpT   = Cfg.getfloat("SubjetAnalysis", "MinFatJetpT")
    self.MinFJ1pT      = Cfg.getfloat("SubjetAnalysis", "MinFJ1pT")
    self.MinFJ2pT      = Cfg.getfloat("SubjetAnalysis", "MinFJ2pT")
    self.MinFJ3pT      = Cfg.getfloat("SubjetAnalysis", "MinFJ3pT")
    self.Min_aJetpt_BST = Cfg.getfloat("SubjetAnalysis", "Min_aJetpt_BST")

    self.BDT_Weights_BST = Cfg.get("SubjetAnalysis", "BDT_Weights_BST")

### BDT
##    self.BDT_MinJetpT1  = Cfg.getfloat("BDT", "MinJetpT1")
##    self.BDT_MinJetpT2  = Cfg.getfloat("BDT", "MinJetpT2")
###
    self.MinLLpTAK5     = Cfg.getfloat("DijetAnalysis", "MinLLpTAK5")
    self.MaxLLpTAK5     = Cfg.getfloat("DijetAnalysis", "MaxLLpTAK5")
    self.MinDiJetpT     = Cfg.getfloat("DijetAnalysis", "MinDiJetpT")
    self.MinJetpT1      = Cfg.getfloat("DijetAnalysis", "MinJetpT1")
    self.MinJetpT2      = Cfg.getfloat("DijetAnalysis", "MinJetpT2")
    self.Min_aJetpt_AK5 = Cfg.getfloat("DijetAnalysis", "Min_aJetpt_AK5")

    self.BDT_Weights_AK5 = Cfg.get("DijetAnalysis", "BDT_Weights_AK5")

### 
    self.MinHVdphi  = Cfg.getfloat("Common", "MinHVdphi")
    self.MinMETJdPhi= Cfg.getfloat("Common", "MinMETJdPhi")

    self.MaxJetEta  = Cfg.getfloat("Common", "MaxJetEta")
    self.MaxAjets   = Cfg.getint("Common", "MaxAjets")
    self.PUWeight   = Cfg.getboolean("Common", "PUWeight")

    self.TriggerBit = Cfg.getint("Common", "TriggerBit")

    self.CSV1       = Cfg.getfloat("Common", "CSV1")
    self.CSV2       = Cfg.getfloat("Common", "CSV2")
    self.CSVL       = Cfg.getfloat("Common", "CSVL")

    self.MinLLMass  = Cfg.getfloat("Common", "MinLLMass")
    self.MaxLLMass  = Cfg.getfloat("Common", "MaxLLMass")

    self.MinMuonPt    = Cfg.getfloat("Common", "MinMuonPt")
    self.MinElePt     = Cfg.getfloat("Common", "MinElePt")
    self.MaxLeptonEta = Cfg.getfloat("Common", "MaxLeptonEta")

def BookHistograms():

    ### Book the histograms

    hists={}

    # currentDir=gFile.pwd()
    # print currentDir

    hName="NStep1"
    hTitle="Number of Step1 Events" 
    hists[hName] = Book1D(hName,hTitle,1,-0.5,0.5,False)

    hName="nPVs"
    hTitle="Number of Primary Vertices" 
    hists[hName] = Book1D(hName,hTitle,30,-0.5,29.5)

    hName="TriggerBits"
    hTitle="Trigger Bits" 
    hists[hName] = Book1D(hName,hTitle,39,-0.5,38.5)

    hName="nFJ"
    hTitle="Number of FilterJets -- V_{pT} > V_{pt}^{cut}" 
    hists[hName] = Book1D(hName,hTitle,6,-0.5,5.5)

    hName="ptFJ1"
    hTitle="FilterJet1 p_{T} -- V_{pT} > V_{pt}^{cut}" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

    hName="ptFJ2"
    hTitle="FilterJet2 p_{T} -- V_{pT} > V_{pt}^{cut}" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

    hName="ptFJ3"
    hTitle="FilterJet3 p_{T} -- V_{pT} > V_{pt}^{cut}" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

    hName="Vtype"
    hTitle="Candidate type" 
    hists[hName] = Book1D(hName,hTitle,6,-0.5,5.5,False)

    hName="HVdPhi"
    hTitle="Delta #phi FatJ and V" 
    hists[hName] = Book1D(hName,hTitle,16,0.,math.pi)

    hName="HMETJdPhi"
    hTitle="Delta #phi FatJ and MET" 
    hists[hName] = Book1D(hName,hTitle,16,0.,math.pi)

    hName="Hmass"
    hTitle="bb Mass" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="HmassPt300"
    hTitle="bb Mass -- p_{T} > 300 GeV" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass2FJ"
    hTitle="bb Mass -- Using top two FilterJets" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass2FJPt300"
    hTitle="bb Mass -- Using top two FilterJets -- p_{T} > 300 GeV" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_2"
    hTitle="bb Mass (from FatH branch)" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="HmassPt300_2"
    hTitle="bb Mass (from FatH branch)  -- p_{T} > 300 GeV" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hpt"
    hTitle="bb p_{T}" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)


    Book_AK5Hists(hists)
    Book_VHists(hists)
    Book_Counters(hists)
    Book_General(hists)
    Book_Comparisons(hists)
    if not jobpar.isData: Book_MCHists(hists)

    outf.cd()

    return hists

def Book_AK5Hists(hists):
    outf.cd()
    outf.mkdir("AK5")
    outf.cd("AK5")

    hName="HVdPhi_AK5"
    hTitle="Delta #phi H(bb) cand and V" 
    hists[hName] = Book1D(hName,hTitle,16,0.,math.pi)

    hName="HVdPhi_AK5_alt"
    hTitle="Delta #phi H(bb) cand and V -- alt calc." 
    hists[hName] = Book1D(hName,hTitle,16,0.,math.pi)

    hName="HMETJdPhi_AK5_alt"
    hTitle="Delta #phi H(bb) cand and MET -- alt calc." 
    hists[hName] = Book1D(hName,hTitle,16,0.,math.pi)

    hName="Hmass_AK5"
    hTitle="bb Mass (from AK5 jet Analysis)" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hpt_AK5"
    hTitle="bb p_{T} (from AK5 jet Analysis)" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

    hName="Hpt1-Hpt2_AK5"
    hTitle="p_{T} Jet1 - p_{T} Jet2 (from AK5 jet Analysis)" 
    hists[hName] = Book1D(hName,hTitle,204,-20.,1000.)


    hName="Hmass_AK5_alt"
    hTitle="bb Mass (from AK5 jet Analysis)" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_pt160_AK5"
    hTitle="bb Mass p_{T} > 160 (from AK5 jet Analysis)" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_pt300_AK5"
    hTitle="bb Mass p_{T} > 300 (from AK5 jet Analysis)" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)


    outf.cd()

def Book_Comparisons(hists):

    outf.cd()
    outf.mkdir("Comparisons")
    outf.cd("Comparisons")

    hName="LeadingElePt"
    hTitle="Leading Electron Pt" 
    hists[hName] = Book1D(hName,hTitle,500, 0., 500.)

    hName="LeadingElePt_ZMass"
    hTitle="Leading Electron Pt -- In Z Mass Window" 
    hists[hName] = Book1D(hName,hTitle,500, 0., 500.)

    hName="LeadingEleEta"
    hTitle="Leading Electron #eta" 
    hists[hName] = Book1D(hName,hTitle,100, -5., 5.)

    hName="LeadingJetPt"
    hTitle="Leading Jet Pt"
    hists[hName] = Book1D(hName,hTitle,1000, 0., 1000.)
        
    hName="2ndLeadingJetPt"
    hTitle="2nd Leading Jet Pt"
    hists[hName] = Book1D(hName,hTitle,1000, 0., 1000.)

    hName="DR_jet_Zelec"
    hTitle="Delta R between leading jet and Z leptons" 
    hists[hName] = Book1D(hName,hTitle,100, 0., 3.)

    hName="NaJets"
    hTitle="Number of Additional Jets  - Jet id applied"
    hists[hName] = Book1D(hName,hTitle,10, 0., 10.)

    for ControlRegion in ["V+udscg", "V+bb", "ttbar"]:
        for Analysis in ["AK5", "Subjet"]:

            hName="Mee_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="e^{+}e^{-} Mass  -- " + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,300,0.,300.)

            hName="Zee_pT_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="Z->e^{+}e^{-} p_{T} -- " + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,500,0.,500.)

            hName="Hmass_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="bb Mass -- " + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

            hName="Hpt_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="bb p_{T} -- " + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

            hName="csv_jet1_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="CSV Jet 1 (Jets ordered by CSV) -- " + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,20,0.,1.)

            hName="csv_jet2_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="CSV Jet 2 (Jets ordered by CSV) -- "  + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,20,0.,1.)

            hName="BDT_" + Analysis + "Analysis_" + ControlRegion + "CR"
            hTitle="BDT Output  -- " + Analysis + "Analysis_" + ControlRegion + "CR"
            hists[hName] = Book1D(hName,hTitle,40,-1.,1.)


    outf.cd()

def Book_MCHists(hists):
    outf.cd()
    outf.mkdir("MC")
    outf.cd("MC")

    hName="genZpt"
    hTitle="Generated Z p{T}" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

    hName="genWpt"
    hTitle="Generated W p{T}" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)

    hName="genZpt_candH"
    hTitle="Generated Z p{T} -- Reconstructed Higgs AK5 Candidate" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.)


    hName="Hmass_50genZpt100"
    hTitle="bb Mass (from AK5 jet Analysis) 50 < Z p_{T} < 100" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_100genZpt150"
    hTitle="bb Mass (from AK5 jet Analysis) 100 < Z p_{T} < 150" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_150genZpt200"
    hTitle="bb Mass (from AK5 jet Analysis) 150 < Z p_{T} < 200" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_200genZpt300"
    hTitle="bb Mass (from AK5 jet Analysis) 200 < Z p_{T} < 300" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)

    hName="Hmass_300genZpt"
    hTitle="bb Mass (from AK5 jet Analysis) Z p_{T} > 300" 
    hists[hName] = Book1D(hName,hTitle,nmbins, minmass, maxmass)


    outf.cd()

def Book_VHists(hists):

    outf.cd()
    outf.mkdir("V")
    outf.cd("V")

    hName="Mll"
    hTitle="ll Mass" 
    hists[hName] = Book1D(hName,hTitle,150,0.,150.)

    hName="Mmumu"
    hTitle="#mu#mu Mass" 
    hists[hName] = Book1D(hName,hTitle,150,0.,150.)

    hName="Mee"
    hTitle="e^{+}e^{-} Mass" 
    hists[hName] = Book1D(hName,hTitle,150,0.,150.)

    hName="Zmumu_pT"
    hTitle="Z->#mu#mu p_{T}" 
    hists[hName] = Book1D(hName,hTitle,500,0.,500.)

    hName="Zee_pT"
    hTitle="Z->e^{+}e^{-} p_{T}" 
    hists[hName] = Book1D(hName,hTitle,500,0.,500.)

    if ( not jobpar.isData):
        hName="Zee_pT_ptGen180"
        hTitle="Z->e^{+}e^{-} p_{T} -- Gen Z p_{T} > 180 GeV" 
        hists[hName] = Book1D(hName,hTitle,500,0.,500.)

        hName="Zee_pT_ptGen200"
        hTitle="Z->e^{+}e^{-} p_{T} -- Gen Z p_{T} > 200 GeV" 
        hists[hName] = Book1D(hName,hTitle,500,0.,500.)

        hName="Zee_pT_ptGen220"
        hTitle="Z->e^{+}e^{-} p_{T} -- Gen Z p_{T} > 220 GeV" 
        hists[hName] = Book1D(hName,hTitle,500,0.,500.)

        hName="Zee_pT_ptGen230"
        hTitle="Z->e^{+}e^{-} p_{T} -- Gen Z p_{T} > 230 GeV" 
        hists[hName] = Book1D(hName,hTitle,500,0.,500.)

        hName="Zee_pT_ptGen240"
        hTitle="Z->e^{+}e^{-} p_{T} -- Gen Z p_{T} > 240 GeV" 
        hists[hName] = Book1D(hName,hTitle,500,0.,500.)

        hName="Zee_pT_ptGen250"
        hTitle="Z->e^{+}e^{-} p_{T} -- Gen Z p_{T} > 250 GeV" 
        hists[hName] = Book1D(hName,hTitle,500,0.,500.)

    outf.cd()

def Book_Counters(hists):

    outf.cd()
    outf.mkdir("Counters")
    outf.cd("Counters")

    hName="Counter"
    hTitle="Event selection counter" 
    hists[hName] = Book1D(hName,hTitle,6,-0.5,5.5,False)
    hists[hName].GetXaxis().SetBinLabel(1,"Events processed")
    hists[hName].GetXaxis().SetBinLabel(2,"Event in JSON")
    hists[hName].GetXaxis().SetBinLabel(3,"Passed Trigger")
    hists[hName].GetXaxis().SetBinLabel(4,"Passed Z mass")
    hists[hName].GetXaxis().SetBinLabel(5,"Passed Light/Heavy jet cut")

    
    for AN in Analyses:

        hName=AN + WhichAnalysis + "Counter" 
        hTitle="Event selection counter -- " + WhichAnalysis + " -- " + AN + " analysis"
        hists[hName] = Book1D(hName,hTitle,11,-0.5,10.5,False)
        if WhichAnalysis == "Zmumu" or WhichAnalysis == "Zee":
            label="Z->ll Mass selection"
        elif WhichAnalysis == "Znn":
            label="Inv. Z selection -- MET passed/ No leptons"
        else:
            label="W selection -- "
            if WhichAnalysis == "Wmun":
                label =label + "Muon pT"
            elif WhichAnalysis == "Wen":
                label =label + "Elec pT and pfMET"
            else:
                print "Undefined analysis type: ", WhichAnalysis
                sys.exit(1)

        hists[hName].GetXaxis().SetBinLabel(1,label)
        if (AN == "AK5"):
            hists[hName].GetXaxis().SetBinLabel(2,"Passed ind. jet pt,eta cuts")
            hists[hName].GetXaxis().SetBinLabel(3,"Passed V pt cut")
            hists[hName].GetXaxis().SetBinLabel(4,"Passed dijet pt cut")
            hists[hName].GetXaxis().SetBinLabel(5,"Passed csv1 cut")
            hists[hName].GetXaxis().SetBinLabel(6,"Passed csv2 cut")
            hists[hName].GetXaxis().SetBinLabel(7,"Passed bb-V dphi cut")
            hists[hName].GetXaxis().SetBinLabel(8,"Passed MET-Jet dphi cut")
            hists[hName].GetXaxis().SetBinLabel(9,"Passed additional jets cut")
        elif (AN == "Boosted"):
            hists[hName].GetXaxis().SetBinLabel(2,"FatHiggs Flag -- 2 filter jets w pt > cut")
            hists[hName].GetXaxis().SetBinLabel(3,"Passed V pt cut")
            hists[hName].GetXaxis().SetBinLabel(4,"Passed fatjet pt cut")
            hists[hName].GetXaxis().SetBinLabel(5,"Passed all filterJet pt cuts")
            hists[hName].GetXaxis().SetBinLabel(6,"Passed csv1 cut")
            hists[hName].GetXaxis().SetBinLabel(7,"Passed csv2 cut")
            hists[hName].GetXaxis().SetBinLabel(8,"Passed bb-V dphi cut")
            hists[hName].GetXaxis().SetBinLabel(9,"Passed MET-Jet dphi cut")
            hists[hName].GetXaxis().SetBinLabel(10,"Passed additional jets cut")


    hName="MC_Counter"
    hTitle="MC Counter" 
    hists[hName] = Book1D(hName,hTitle,5,-0.5,4.5,False)
    hists[hName].GetXaxis().SetBinLabel(1,"Events processed")
    hists[hName].GetXaxis().SetBinLabel(2,"NAK5 > 1")
    hists[hName].GetXaxis().SetBinLabel(3,"NAK5 (pt>20) > 1")

    outf.cd()

def Book_General(hists):
    outf.cd()
    outf.mkdir("General")
    outf.cd("General")

    hName="Json"
    hTitle="Json contains Event"
    hists[hName] = Book1D(hName,hTitle,2,0.,2.,False)
    
    hName="METet"
    hTitle="e_{T} -- MET" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.,False)
    hName="METnoPUet"
    hTitle="e_{T} -- MET" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.,False)
    hName="METnoPUChet"
    hTitle="e_{T} -- MET" 
    hists[hName] = Book1D(hName,hTitle,200,0.,1000.,False)


    outf.cd()

def checkForV():

    ## return string

    VBoson="notV"

    if tree.Vtype == 0 or tree.Vtype == 1: 
        if tree.Vtype==0:
            ltype=13
        else:
            ltype=11
        if InZMassWindow() and nVLepton(ltype) >1 : VBoson="Z"
    elif tree.Vtype == 2 or tree.Vtype == 3: 
        if tree.Vtype == 2: # Wmunu
            if getMaxVLeptonPt(13)> jobpar.MinMuonPt: VBoson="W"
        elif tree.Vtype == 3: # Wenu
            if getMaxVLeptonPt(11)> jobpar.MinElePt and METnoPUCh.et>35: VBoson="W"
    elif tree.Vtype == 4:  # Znn
        if nVLepton(11) ==0 and nVLepton(13) ==0 and MET.et>jobpar.MinMET: VBoson="Z"

    return VBoson


def getMaxVLeptonPt(ltype):

    # ltype is pdg code for lepton
    ptmax=0.
    for i in xrange( tree.nvlep ):
        if tree.vLepton_type[i] == ltype: 
            if tree.vLepton_pt[i] >  ptmax: ptmax=tree.vLepton_pt[i]

    return ptmax

def goodVLepton(ltype,i):
    # ltype is pdg code for lepton
    if (ltype==11):
        MinPt=jobpar.MinElePt
    elif (ltype==13):
        MinPt=jobpar.MinMuonPt
    else:
        MinPt=9999.
        print "Problem with ltype"
        return False

    if tree.vLepton_type[i] != ltype: return False
    if tree.vLepton_pt[i] < MinPt : return False
    if math.fabs(tree.vLepton_eta[i])>1.44 and math.fabs(tree.vLepton_eta[i])<1.57: return False
    if math.fabs(tree.vLepton_eta[i])>jobpar.MaxLeptonEta: return False
    ## if math.fabs(tree.vLepton_id95[i] -7.0) > 0.001: return False
    leptonID=tree.vLepton_id2012tight[i]

    return leptonID

def nVLepton(ltype):

    n=0
    for i in xrange( tree.nvlep ):
        if goodVLepton(ltype,i): n=n+1  
        
    return n

def InZMassWindow():

    ## require ll mass within Z mass window, no pTMin requirement at this point

    if tree.Vtype > 1: return False
    if V.mass < jobpar.MinLLMass or V.mass > jobpar.MaxLLMass : return False

    return True

def fill_ZHists():

    Vtype=tree.Vtype
    if Vtype > 1: return

    h["Mll"].Fill(V.mass)
    if Vtype == 0: 
        h["Mmumu"].Fill(V.mass)
    else:
        h["Mee"].Fill(V.mass)

    if V.mass < jobpar.MinLLMass or V.mass > jobpar.MaxLLMass : return

    if Vtype == 0: 
        h["Zmumu_pT"].Fill(V.pt)
    else:
        h["Zee_pT"].Fill(V.pt)
        if ( not jobpar.isData):
            if tree.genZpt >180.:
                h["Zee_pT_ptGen180"].Fill(V.pt)
            if tree.genZpt >200.:
                h["Zee_pT_ptGen200"].Fill(V.pt)
            if tree.genZpt >220.:
                h["Zee_pT_ptGen220"].Fill(V.pt)
            if tree.genZpt >230.:
                h["Zee_pT_ptGen230"].Fill(V.pt)
            if tree.genZpt >240.:
                h["Zee_pT_ptGen240"].Fill(V.pt)
            if tree.genZpt >250.:
                h["Zee_pT_ptGen250"].Fill(V.pt)


    return

def ApplyAK5SelectionCuts():

    if (tree.Vtype != CandidateDict[WhichAnalysis]):    return False  ## require that Vtype matches analysis
    # if (tree.Vtype > 1):    return False  ## require Z =ee,mm 

    VBoson=checkForV() # apply mass selection cut
    if (VBoson == "notV"): return False

    hCounterName="AK5" + WhichAnalysis + "Counter" 

    h[hCounterName].Fill(0.) ## number of events entering Z analysis


    if (tree.hJet_pt[0]<jobpar.MinJetpT1 or tree.hJet_pt[1]<jobpar.MinJetpT2): return False
    if (math.fabs(tree.hJet_eta[0])>jobpar.MaxJetEta or math.fabs(tree.hJet_eta[1])>jobpar.MaxJetEta): return False
    h[hCounterName].Fill(1.)


    ## Apply V pT cut
    if V.pt < jobpar.MinLLpTAK5 or V.pt > jobpar.MaxLLpTAK5: return False
    h[hCounterName].Fill(2.)

    ## Dijet pT cut
    if H.pt < jobpar.MinDiJetpT: return False
    h[hCounterName].Fill(3.)

    ### Combined secondary vertex cut
    if (tree.hJet_csv[0]<jobpar.CSV1 and tree.hJet_csv[1]<jobpar.CSV1): return False  # at least one btag>0.9
    h[hCounterName].Fill(4.)
    if (tree.hJet_csv[0]>=jobpar.CSV1):
        if (tree.hJet_csv[1]<jobpar.CSV2): return False  # other tag has to be greater than  0.5
    else:
        if (tree.hJet_csv[0]<jobpar.CSV2): return False
    h[hCounterName].Fill(5.)

    h["HVdPhi_AK5"].Fill(tree.HVdPhi)

    if tree.Vtype ==4:
        dPhi=math.fabs(deltaPhi(H.phi,MET.phi))
    else:
        dPhi=math.fabs(deltaPhi(H.phi,V.phi))
    h["HVdPhi_AK5_alt"].Fill(dPhi)

    if tree.HVdPhi< jobpar.MinHVdphi: return False
    h[hCounterName].Fill(6.)

    METJ1dPhi=math.fabs(deltaPhi(MET.phi,tree.hJet_phi[0]))
    METJ2dPhi=math.fabs(deltaPhi(MET.phi,tree.hJet_phi[1]))
    if METJ1dPhi>METJ2dPhi:
         METJdPhi=METJ1dPhi 
    else:
         METJdPhi=METJ2dPhi  
         
    h["HMETJdPhi_AK5_alt"].Fill(METJdPhi)

    ## Cut on METJdPhi for Znnu Analysis
    if METJdPhi<jobpar.MinMETJdPhi: return False
    h[hCounterName].Fill(7.)


    ## Cut on number of additional jets for other analyses
    naJets=0
    for i in xrange( tree.naJets ):
        if tree.aJet_pt[i]>jobpar.Min_aJetpt_AK5 and math.fabs(tree.aJet_eta[i]<jobpar.MaxJetEta): naJets=naJets+1
    if naJets > jobpar.MaxAjets: return False
    h[hCounterName].Fill(8.)

    return True

def findNAK5Jets(eta0,phi0,ptMin):

    ## find number of AK5 jets with pT>ptMin and located drCut away from a given jet
    ## addional jet

    drCut=1.2

    n=0
    # first loop over nhjets
    for i in xrange( tree.nhJets ):
        if tree.hJet_pt[i]>ptMin and \
                math.fabs(tree.hJet_eta[i])<jobpar.MaxJetEta and \
                deltaR(eta0,phi0,tree.hJet_eta[i],tree.hJet_phi[i]) > drCut: 
            n=n+1

    # now over the additional jets
    for i in xrange( tree.naJets ):
        if tree.aJet_pt[i]>ptMin and \
                math.fabs(tree.aJet_eta[i])<jobpar.MaxJetEta and \
                deltaR(eta0,phi0,tree.aJet_eta[i],tree.aJet_phi[i]) > drCut: 
            n=n+1

    return n

def findNAK5Jets_fromFilterJets(ptMin):

    ## find number of AK5 jets with pT>ptMin and located drCut away from the two leading filterJets
    drCut=0.3

    n=0
    # first loop over nhjets
    for i in xrange( tree.nhJets ):
        if tree.hJet_pt[i]>ptMin and \
                math.fabs(tree.hJet_eta[i])<jobpar.MaxJetEta and \
                deltaR(tree.fathFilterJets_eta[0],tree.fathFilterJets_phi[0],tree.hJet_eta[i],tree.hJet_phi[i]) > drCut and \
                deltaR(tree.fathFilterJets_eta[1],tree.fathFilterJets_phi[1],tree.hJet_eta[i],tree.hJet_phi[i]) > drCut:
            n=n+1

    # now over the additional jets
    for i in xrange( tree.naJets ):
        if tree.aJet_pt[i]>ptMin and \
                math.fabs(tree.aJet_eta[i])<jobpar.MaxJetEta and \
                deltaR(tree.fathFilterJets_eta[0],tree.fathFilterJets_phi[0],tree.aJet_eta[i],tree.aJet_phi[i]) > drCut and \
                deltaR(tree.fathFilterJets_eta[1],tree.fathFilterJets_phi[1],tree.aJet_eta[i],tree.aJet_phi[i]) > drCut: 
            n=n+1

    print n
    return n


def ApplyBoostedANSelectionCuts():

    if (tree.Vtype != CandidateDict[WhichAnalysis]):    return False  ## check that the Vtype matches analysis
    # if (tree.Vtype > 1):    return False  ## require Z =ee,mm 

    VBoson=checkForV() # apply mass selection cut

    if (VBoson == "notV"): return False

    hCounterName="Boosted" + WhichAnalysis + "Counter" 

    h[hCounterName].Fill(0.) ## number of events entering Z analysis

    if (not FatH.FatHiggsFlag): return False  # found fatjet with at least 2 filterjets w pt>30
    h[hCounterName].Fill(1.)

    ## Apply V pT cut
    if V.pt < jobpar.MinLLpT_BST or V.pt > jobpar.MaxLLpT_BST: return False
    h[hCounterName].Fill(2.)

    if (FatH.pt< jobpar.MinFatJetpT or math.fabs(FatH.eta) > jobpar.MaxJetEta): return False
    h[hCounterName].Fill(3.)

    nfathFilterJets=int(tree.nfathFilterJets)
    h["nFJ"].Fill(nfathFilterJets)

    if (nfathFilterJets < 2) : return False
    h["ptFJ1"].Fill(tree.fathFilterJets_pt[0])
    h["ptFJ2"].Fill(tree.fathFilterJets_pt[1])
    if (nfathFilterJets > 2):
        h["ptFJ3"].Fill(tree.fathFilterJets_pt[2])

    if (tree.fathFilterJets_pt[0] < jobpar.MinFJ1pT): return False
    if (tree.fathFilterJets_pt[1] < jobpar.MinFJ2pT): return False
    if (nfathFilterJets > 2):
        if (tree.fathFilterJets_pt[2] < jobpar.MinFJ3pT): return False
    h[hCounterName].Fill(4.)

    if (tree.fathFilterJets_pt[0]< tree.fathFilterJets_pt[1]):
        print "A: FilterjetPTs not sorted: ",tree.fathFilterJets_pt[0],tree.fathFilterJets_pt[1]
    if (nfathFilterJets > 2) and (tree.fathFilterJets_pt[0]< tree.fathFilterJets_pt[2]):
        print "B: FilterjetPTs not sorted: ",tree.fathFilterJets_pt[0],tree.fathFilterJets_pt[2]

    ### Combined secondary vertex cut
    if (tree.fathFilterJets_csv[0]<jobpar.CSV1 and tree.fathFilterJets_csv[1]<jobpar.CSV1): return False  # at least one btag>0.9
    h[hCounterName].Fill(5.)
    if (tree.fathFilterJets_csv[0]>=jobpar.CSV1):
        if (tree.fathFilterJets_csv[1]<jobpar.CSV2): return False  # other tag has to be greater than  0.5
    else:
        if (tree.fathFilterJets_csv[0]<jobpar.CSV2): return False
    h[hCounterName].Fill(6.)

    if tree.Vtype ==4:
        dPhi=math.fabs(deltaPhi(FatH.phi,MET.phi))
    else:
        dPhi=math.fabs(deltaPhi(FatH.phi,V.phi))
    h["HVdPhi"].Fill(dPhi)
    if dPhi< jobpar.MinHVdphi: return False
    h[hCounterName].Fill(7.)

    METJ1dPhi=math.fabs(deltaPhi(MET.phi,tree.fathFilterJets_phi[0]))
    METJ2dPhi=math.fabs(deltaPhi(MET.phi,tree.fathFilterJets_phi[1]))
    if METJ1dPhi>METJ2dPhi:
        METJdPhi=METJ1dPhi 
    else:
        METJdPhi=METJ2dPhi   
    h["HMETJdPhi"].Fill(METJdPhi)

    if METJdPhi<jobpar.MinMETJdPhi: return False
    h[hCounterName].Fill(8.)

    ## Number of additional jets
    naJets=0
    for i in xrange( tree.naJets ):
        if tree.aJet_pt[i]>jobpar.Min_aJetpt_BST: naJets=naJets+1
    
    if naJets <= jobpar.MaxAjets: 
        h[hCounterName].Fill(8.)

    naJets=findNAK5Jets(FatH.eta,FatH.phi,jobpar.Min_aJetpt_BST)
    # naJets=findNAK5Jets_fromFilterJets(jobpar.Min_aJetpt_BST)
    if naJets > jobpar.MaxAjets: return False

    h[hCounterName].Fill(9.)


    return True

def deltaPhi(phi1,phi2):
                      
    result=phi1-phi2
    while (result > math.pi): result -= 2*math.pi
    while (result <= -math.pi): result += 2*math.pi

    return result

def deltaR(eta1,phi1,eta2,phi2):
    deta=eta1-eta2
    dphi=deltaPhi(phi1,phi2)
                      
    return math.sqrt(deta*deta + dphi*dphi)

def GetCSVs(Analysis):

    csv1=-999.
    csv2=-999.

    csv=[]
    if Analysis == "AK5":
        for c in tree.hJet_csv:
            csv.append(c)
        nJets=int(tree.nhJets)
    else:
        for c in tree.fathFilterJets_csv:
            csv.append(c)
        nJets=int(tree.nfathFilterJets)

    csv.sort(reverse=True)
    if nJets>0:
        csv1=csv[0]
        if nJets>1:
            csv2=csv[1]

    return csv1,csv2


def ZllComparisonPlots(wt=1.):

    if (tree.Vtype != CandidateDict[WhichAnalysis]):    return ## require that Vtype matches analysis
    VBoson="notZ"
    if tree.Vtype == 0 or tree.Vtype == 1: 
        if tree.Vtype==0:
            ltype=13
        else:
            ltype=11
        if InZMassWindow() and nVLepton(ltype) >1 : VBoson="Z"
    else:
        return
    
    if tree.vLepton_pt[1] > tree.vLepton_pt[0]:
        print "Leading pT Electron is 2nd electron!!!"
        
    if goodVLepton(ltype,0):
        h["LeadingElePt"].Fill(tree.vLepton_pt[0],wt)
        h["LeadingEleEta"].Fill(tree.vLepton_eta[0],wt)
        if InZMassWindow(): h["LeadingElePt_ZMass"].Fill(tree.vLepton_pt[0],wt)


    dr1=deltaR(tree.vLepton_eta[0],tree.vLepton_phi[0],tree.hJet_eta[0],tree.hJet_phi[0])
    dr2=deltaR(tree.vLepton_eta[1],tree.vLepton_phi[1],tree.hJet_eta[0],tree.hJet_phi[0])
    h["DR_jet_Zelec"].Fill(dr1)
    h["DR_jet_Zelec"].Fill(dr2)


    if IsGood_hJet(0):
        h["LeadingJetPt"].Fill(tree.hJet_pt[0],wt)

    if IsGood_hJet(1):
        h["2ndLeadingJetPt"].Fill(tree.hJet_pt[1],wt)

    naJets=0
    for i in xrange( tree.naJets ):
        if IsGood_aJet(i): naJets=naJets+1

    h["NaJets"].Fill(naJets,wt)


    ## VBoson with good lepton ID, in V mass window, and pt>ptmin 
    ## Trigger satified

    # print bdt_Hmass[0], bdt_Hpt[0], bdt_Vpt[0]
    # bdt=reader.EvaluateMVA("BDT")


    # print "BDT: ",bdt,reader_ak5.bdt_Hmass[0], reader_ak5.bdt_Hpt[0], reader_ak5.bdt_Vpt[0]

    for Analysis in ["AK5", "Subjet"]:

        bdt=-999.
        if Analysis == 'AK5':
            naJets=tree.naJets
            dPhi=tree.HVdPhi
            Hpt=H.pt
            Hmass=H.mass
            if (H.HiggsFlag > 0):
                bdt=reader_ak5.evaluate()

        elif Analysis == 'Subjet':
               
            naJets=findNAK5Jets(FatH.eta,FatH.phi,jobpar.Min_aJetpt_BST)
            dPhi=math.fabs(deltaPhi(FatH.phi,V.phi))
            Hpt=FatH.filteredpt
            Hmass=FatH.filteredmass
            if (FatH.FatHiggsFlag > 0 and tree.nfathFilterJets >= 2):
                bdt=reader_bst.evaluate()
        else:
            print "Undefined Control Region"
            sys.exit(1)

        csv1, csv2 = GetCSVs(Analysis)

        for ControlRegion in ["V+udscg", "V+bb", "ttbar"]:

            inControlRegion=False
            if ControlRegion == "V+udscg": ## udscg control region (Table 13 CMS AN-2011/430)
                if VBoson == "Z" and V.pt > jobpar.MinLLpTAK5 and naJets < 2 and dPhi > jobpar.MinHVdphi and csv1 < jobpar.CSV1: 
                    inControlRegion=True

            elif ControlRegion == "V+bb":
                if VBoson == "Z" and V.pt > jobpar.MinLLpTAK5 and naJets < 2 and dPhi > jobpar.MinHVdphi and csv1 > CSVT and csv2 > 0.5 and (Hmass < 90. or Hmass > 145.): 
                    inControlRegion=True

            elif ControlRegion == "ttbar":
                if nVLepton(ltype) > 1 and (not InZMassWindow()) and csv1 > CSVT and csv2 > 0.5 and Hpt> 100.:
                    inControlRegion=True

            if inControlRegion:
                
                ## comparision plots 

                h["Mee_"+ Analysis + "Analysis_" + ControlRegion + "CR"].Fill(V.mass,wt)    
                h["Zee_pT_" + Analysis + "Analysis_" + ControlRegion + "CR"].Fill(V.pt,wt)

                h["BDT_" + Analysis + "Analysis_" + ControlRegion + "CR"].Fill(bdt,wt)

                if (Analysis=="AK5" and IsGood_hJet(0) and IsGood_hJet(1)) or \
                        (Analysis=="Subjet" and IsGood_fJet(0) and IsGood_fJet(1) and IsGood_fJet(2)) :
                        # satifies jet ID and pt requirements

                    h["Hmass_" + Analysis + "Analysis_" + ControlRegion + "CR"].Fill(Hmass,wt)
                    h["Hpt_" + Analysis + "Analysis_" + ControlRegion + "CR"].Fill(Hpt,wt)

                    h["csv_jet1_" + Analysis + "Analysis_" + ControlRegion + "CR"].Fill(csv1,wt)
                    h["csv_jet2_" + Analysis + "Analysis_" + ControlRegion + "CR"].Fill(csv2,wt)


def IsGood_fJet(indx):

    ## satisfies jet pt,eta and id requirements

    if indx >= len(tree.fathFilterJets_pt): return False
    
    if indx==0:
        if (tree.fathFilterJets_pt[0] < jobpar.MinFJ1pT): return False

    elif indx==1:
        if (tree.fathFilterJets_pt[1] < jobpar.MinFJ2pT): return False

    elif indx==2:
        if (tree.nfathFilterJets > 2):
            if (tree.fathFilterJets_pt[2] < jobpar.MinFJ3pT): return False
    else:
        print "Bad filter jet index"
        return False


    if math.fabs(tree.fathFilterJets_eta[indx]) > jobpar.MaxJetEta: return False 
    if tree.fathFilterJets_chf[indx]>0.99:  return False

    return True

def IsGood_hJet(indx):

    
    ## if math.fabs(tree.hJet_nhf[indx]) > 0.99: return False 
    ## if math.fabs(tree.hJet_nef[indx]) > 0.99: return False 
    ## if tree.hJet_nconstituents < 2: return False
    ## 
    ## if (math.fabs(tree.hJet_eta[indx]) < 2.4):
    ##     if (tree.hJet_cef[indx] > 0.99): return False
    ##     if (tree.hJet_chf[indx] == 0): return False
    ##     if (tree.hJet_nch[indx] == 0): return False


    if math.fabs(tree.hJet_eta[indx]) > jobpar.MaxJetEta: return False 
    if tree.hJet_pt[indx]< jobpar.MinJetpT1: return False

    jid=False

    if indx<0 or indx>1:
        print "Illegal index passed"
        return False

    if indx==0:
        if len(tree.hJet_id)>0:
            jid=bool(tree.hJet_id[0])
    elif indx==1:
        if len(tree.hJet_id)>1:
            jid=bool(tree.hJet_id[1])

    return jid


def IsGood_aJet(indx):

    ## if tree.aJet_pt[indx]< jobpar.MinJetpT1: return False
    ## if math.fabs(tree.aJet_eta[indx]) > jobpar.MaxJetEta: return False 
    ## if math.fabs(tree.aJet_nhf[indx]) > 0.99: return False 
    ## if math.fabs(tree.aJet_nef[indx]) > 0.99: return False 
    ## if tree.aJet_nconstituents < 2: return False
    ## 
    ## if (math.fabs(tree.aJet_eta[indx]) < 2.4):
    ##     if (tree.aJet_cef[indx] > 0.99): return False
    ##     if (tree.aJet_chf[indx] == 0): return False
    ##     if (tree.aJet_nch[indx] == 0): return False
    ## 
    ## return True

    jid=False

    if indx<0:
        print "Illegal index passed"
        return False

    if len(tree.aJet_id)<=indx:
        print "Bug in aJet_id variable -- to be fixed. Returning False for now" 
        return False
    else:
        jid=bool(tree.aJet_id[indx])

    return jid

def AK5JetsAnalysis(wt=1.):

    if len(tree.hJet_pt)>1:
        h["Hpt1-Hpt2_AK5"].Fill(tree.hJet_pt[0]-tree.hJet_pt[1],wt)

    EvPassed=ApplyAK5SelectionCuts()
    if not EvPassed : return

    # alternative mass calc.
    B1=TLorentzVector()
    B2=TLorentzVector()
    B1.SetPtEtaPhiE(tree.hJet_pt[0],tree.hJet_eta[0],tree.hJet_phi[0],tree.hJet_e[0])
    B2.SetPtEtaPhiE(tree.hJet_pt[1],tree.hJet_eta[1],tree.hJet_phi[1],tree.hJet_e[1])
    mass=(B1+B2).M()

    h["Hmass_AK5"].Fill(H.mass,wt)
    h["Hpt_AK5"].Fill(H.pt,wt)

    h["Hmass_AK5_alt"].Fill(mass,wt)

    if (H.pt >160.): 
        h["Hmass_pt160_AK5"].Fill(H.mass,wt)

    if (H.pt >300.): 
        h["Hmass_pt300_AK5"].Fill(H.mass,wt)

def ApplyAK3Residual(xx):

    if not jobpar.ApplyJECResidual:
        return 1

    x=xx
    if x<10: x=10
    if x>600: x=600

    par=[-1932.24059306, -517.206840357, 120.326645168, 
          0.0993702488393, -54.7280305025, 184.288787916,
          0.971859290257]

    value1 = par[0]*math.exp(-0.5 * (x - par[1])**2 / par[2]**2 );
    value2 = par[3]*math.exp(-0.5 * (x - par[4])**2 / par[5]**2 );

    return par[6]-(value1+value2)


def BoostedAnalysis(wt=1.):

    EvPassed=ApplyBoostedANSelectionCuts()
    if not EvPassed : return

#   filterJets
    nfathFilterJets=int(tree.nfathFilterJets)
    fathFilterJets_pt= tree.fathFilterJets_pt
    fathFilterJets_eta=tree.fathFilterJets_eta
    fathFilterJets_phi=tree.fathFilterJets_phi
    fathFilterJets_e=  tree.fathFilterJets_e
    fathFilterJets_csv=tree.fathFilterJets_csv

    F1=TLorentzVector()
    F2=TLorentzVector()
    F3=TLorentzVector()
    m=0
    pt=0

    if nfathFilterJets>1:

        corr1=ApplyAK3Residual( fathFilterJets_pt[0] )
        corr2=ApplyAK3Residual( fathFilterJets_pt[1] )

        # F1.SetPtEtaPhiE(fathFilterJets_pt[0],fathFilterJets_eta[0],fathFilterJets_phi[0],fathFilterJets_e[0])
        # F2.SetPtEtaPhiE(fathFilterJets_pt[1],fathFilterJets_eta[1],fathFilterJets_phi[1],fathFilterJets_e[1])
        
        F1.SetPtEtaPhiE(fathFilterJets_pt[0]/corr1,fathFilterJets_eta[0],fathFilterJets_phi[0],fathFilterJets_e[0]/corr1)
        F2.SetPtEtaPhiE(fathFilterJets_pt[1]/corr2,fathFilterJets_eta[1],fathFilterJets_phi[1],fathFilterJets_e[1]/corr2)
        m=(F1+F2).M()
        pt=(F1+F2).Pt()

        m2j=m
        pt2j=pt

        # print fathFilterJets_pt[0],corr,corrpt


        if nfathFilterJets>2:
        ## if nfathFilterJets>2 and fathFilterJets_pt[2]> 20:
            # F3.SetPtEtaPhiE(fathFilterJets_pt[2],fathFilterJets_eta[2],fathFilterJets_phi[2],fathFilterJets_e[2])
            corr3=ApplyAK3Residual( fathFilterJets_pt[2] )
            F3.SetPtEtaPhiE(fathFilterJets_pt[2]/corr3,fathFilterJets_eta[2],fathFilterJets_phi[2],fathFilterJets_e[2]/corr3)

            m=(F1+F2+F3).M()
            pt=(F1+F2+F3).Pt()

    if (m>0): 

        h["Hmass"].Fill(m,wt)
        h["Hpt"].Fill(pt,wt)
        h["Hmass_2"].Fill(FatH.filteredmass,wt)

        h["Hmass2FJ"].Fill(m2j,wt)

        if pt>300:
            h["HmassPt300"].Fill(m,wt)
            h["Hmass2FJPt300"].Fill(m2j,wt)
            h["HmassPt300_2"].Fill(FatH.filteredmass,wt)


def GeneratedInfo():

    # if tree.genZpt<50: return

    h["MC_Counter"].Fill(0.)  # count the number processed
    h["genZpt"].Fill(tree.genZpt)

    h["genWpt"].Fill(tree.genWpt)

    if tree.nhJets>1:
        h["MC_Counter"].Fill(1.)  
        if (tree.hJet_pt[0]>20. and tree.hJet_pt[1]>20):

            if tree.genZpt<100:
                h["Hmass_50genZpt100"].Fill(H.mass)
            elif tree.genZpt<150:
                h["Hmass_100genZpt150"].Fill(H.mass)
            elif tree.genZpt<200:
                h["Hmass_150genZpt200"].Fill(H.mass)
            elif tree.genZpt<300:
                h["Hmass_200genZpt300"].Fill(H.mass)
            else:
                h["Hmass_300genZpt"].Fill(H.mass)

            h["MC_Counter"].Fill(2.)
            if (H.mass>80 and H.mass<150):
                h["MC_Counter"].Fill(3.)  # count the number processed
                h["genZpt_candH"].Fill(tree.genZpt)



def processEvent():

    ## if jobpar.isData and EVENT.run < 175860:return

    h["Counter"].Fill(0.)  # count the number processed
    h["Json"].Fill(float(EVENT.json))
    if not EVENT.json: return
    h["Counter"].Fill(1.)  # count the number with valid json
    # tree.Print()
    
    if (tree.nhJets<2):
        print "njJets<2. -- This should not happen"
        return

    # nt=0
    # s=" "
    for itrig in xrange( len(triggers) ):
        #nt=nt+triggers[itrig]
        #s=s+str(triggers[itrig])
        if triggers[itrig]>0: h["TriggerBits"].Fill(itrig)

    # print "XXXX ",s
    # print "Trig Sum: ",nt
    if jobpar.TriggerBit > -1:
        # print jobpar.TriggerBit,triggers[jobpar.TriggerBit],triggers[5]
        # if triggers[jobpar.TriggerBit]==0: return
        if triggers[5]==0 and triggers[6]==0: return
    h["Counter"].Fill(2.)


    # load the BDT
    reader_ak5.loadVar()
    reader_bst.loadVar()

    
    fill_ZHists()
    if InZMassWindow(): h["Counter"].Fill(3.)

    if jobpar.heavyJets:
        if tree.eventFlav != 5: return

    if jobpar.lightJets:
        if tree.eventFlav == 5: return

    h["Counter"].Fill(4.)

    ## print tree.
    wt=1.  ## overall event weight
    ## if jobpar.PUWeight: wt=wt*tree.PUweight
    if jobpar.PUWeight: 
        if (jobpar.MC7TeV):
            wt=wt*tree.PUweight2011B
        else:
            wt=wt*tree.PUweight

    if jobpar.isData:
        ## apply trigger weight
        if tree.weightTrig2012ADiEle>0.001:
            wt=wt/tree.weightTrig2012ADiEle
        else:
            wt=0.

    nPVs=int(tree.nPVs)
    Vtype=int(tree.Vtype)

    h["nPVs"].Fill(nPVs,wt)
    # print nPVs

    h["Vtype"].Fill(Vtype)

    h["METet"].Fill(MET.et,wt)
    h["METnoPUet"].Fill(METnoPU.et,wt)
    h["METnoPUChet"].Fill(METnoPUCh.et,wt)

    if not jobpar.isData: GeneratedInfo()


    ZllComparisonPlots(wt)
    AK5JetsAnalysis(wt)
    BoostedAnalysis(wt)


def OpenFile(file_in,iodir):
    """  file_in -- Input file name
         iodir   -- 'r' readonly  'r+' read+write """
    try:
        ifile=open(file_in, iodir)
        # print "Opened file: ",file_in," iodir ",iodir
    except:
        print "Could not open file: ",file_in
        sys.exit(1)
    return ifile

def CloseFile(ifile):
    ifile.close()


def ReadFile(file):


    infile=OpenFile(file,'r')
    iline=0

    x = infile.readline()

    totevents=0
    while x != "":
        iline+=1
        xx=string.rstrip(x)
        inputLine=string.split(xx)
        print xx,len(inputLine)
        events=int(inputLine[4])
        totevents=totevents+events
        print inputLine[1], inputLine[4], inputLine[13], totevents
        
        x = infile.readline()


if __name__ == "__main__":

    jobpar=JobPar()
    print "Analysis: ",jobpar.WhichAnalysis
    print "Data: ", jobpar.isData
    print "CSV1: " ,jobpar.CSV1
    print "CSV2: " ,jobpar.CSV2

    if jobpar.isData:
        if jobpar.PUWeight == True:
            print "Turning off vertex weighting for data"
            jobpar.PUWeight=False
    else:
        if jobpar.MC7TeV:
            print "Running 7 TeV MC"
        else:
            print "Running 8 TeV MC"

    WhichAnalysis=jobpar.WhichAnalysis

    AutoLibraryLoader.enable()


    if not CandidateDict.has_key(WhichAnalysis):
        print "Not setup for ",WhichAnalysis," Analysis"
        print "Options are: ",CandidateDict.keys()
        sys.exit(1)
        

    # create output root file and book histograms
    outf=OutputRootFile(jobpar.outFile)

    h=BookHistograms()

    tree, nStep1 =getRootChain(jobpar.inputFiles,"tree")
    h["NStep1"].Fill(0,nStep1)

    if WhichAnalysis != "Zmumu" and WhichAnalysis != "Zee" and \
            WhichAnalysis != "Znn" and WhichAnalysis != "Wen" and WhichAnalysis != "Wmun":
        print "Not setup for ",WhichAnalysis," Analysis"
        sys.exit(1)

    leaves = tree.GetListOfBranches();
    # printLeaves(leaves)

    FatH   = ROOT.FatHiggsInfo()
    H      = ROOT.HiggsInfo()
    V      = ROOT.VInfo()
    EVENT  = ROOT.EventInfo()
    MET      = ROOT.METInfo()
    METnoPU  = ROOT.METInfo()
    METnoPUCh= ROOT.METInfo()

    tree.SetBranchAddress("FatH",  AddressOf(FatH, "FatHiggsFlag") );
    tree.SetBranchAddress("H",     AddressOf(H, "HiggsFlag") );
    tree.SetBranchAddress("V",     AddressOf(V, "mass") );
    tree.SetBranchAddress("EVENT", AddressOf(EVENT, "run") );
    tree.SetBranchAddress("MET"      , AddressOf(MET, "et") );
    tree.SetBranchAddress("METnoPU"  , AddressOf(METnoPU, "et") );
    tree.SetBranchAddress("METnoPUCh", AddressOf(METnoPUCh, "et") );

    triggers = defineArray( 'b' , 39 )
    tree.SetBranchAddress("triggerFlags",triggers)
    # print triggers

    NEntries = tree.GetEntries()
    print "Number of entries on Tree:",NEntries
    nevt=NEntries
    if jobpar.nevts>-1:
        nevt=jobpar.nevts
    print "Processing: ",nevt," events"

    # setup BDT
    reader_ak5 = BDT_READER_AK5(jobpar.BDT_Weights_AK5)
    reader_bst = BDT_READER_BST(jobpar.BDT_Weights_BST)

    ## reader = TMVA.Reader("Color:!Silent");
    ## 
    ## bdt_Hmass = array( 'f', [ 0 ] )
    ## bdt_Hpt = array( 'f', [ 0 ] )
    ## bdt_Vpt = array( 'f', [ 0 ] )
    ## bdt_csv0 = array( 'f', [ 0 ] )
    ## bdt_csv1 = array( 'f', [ 0 ] )
    ## bdt_HVdPhi = array( 'f', [ 0 ] )
    ## bdt_dEta = array( 'f', [ 0 ] )
    ## 
    ## bdt_hpt0 = array( 'f', [ 0 ] )
    ## bdt_hpt1 = array( 'f', [ 0 ] )
    ## 
    ## reader.AddVariable("H_mass := H.mass", bdt_Hmass );
    ## reader.AddVariable("H_pt := H.pt", bdt_Hpt );
    ## reader.AddVariable("V_pt :=V.pt", bdt_Vpt );
    ## reader.AddVariable("hJ12_MaxCsv := max(hJet_csv[0],hJet_csv[1])", bdt_csv0 );
    ## reader.AddVariable("hJ12_MinCsv := min(hJet_csv[0],hJet_csv[1])", bdt_csv1 );
    ## reader.AddVariable("HV_dPhi := HVdPhi", bdt_HVdPhi );
    ## reader.AddVariable("H_dEta := H.dEta", bdt_dEta );
    ## 
    ## reader.AddSpectator("hJet_pt[0]", bdt_hpt0 );
    ## reader.AddSpectator("hJet_pt[1]", bdt_hpt1 );
    ## 
    ## reader.BookMVA("BDT", jobpar.BDT_Weights_AK5);

    decade=0
    for jentry in xrange( nevt ):

        progress = 10.0*jentry/(1.0*nevt);
        k = math.floor(progress); 
        if (k > decade):
            print 10*k," %",jentry, "\tRun:Event:lumi:json: ",EVENT.run, EVENT.event, EVENT.lumi, EVENT.json
        decade = k  

        tree.GetEntry(jentry) 
        processEvent()  ### work is done here


    outf.Write();
    outf.Close();
