import sys,string,math,os
import ConfigParser

import ROOT
from ROOT import TLorentzVector
from ROOT import gDirectory, gROOT, gFile, AddressOf
from ROOT import AutoLibraryLoader

from myRootIOFuncs import *


######################################


## max dist for jet matching
drMatch=0.25

## Define mass histogram binning here
nmbins=60
minmass=0.
maxmass=300.

ROOT.gROOT.LoadMacro('jetStruct.h+')

########################################

def usage():
    """ Usage: pyHiggsAnalyzer configFile

    """
    pass

class JobPar:
    def __init__(self):
        getJobopt(self)

def getJobopt(self):

    narg=len(sys.argv)
    if narg < 2 :
        print usage.__doc__
        sys.exit(1)

    cfgFile= sys.argv[1]

    if not os.path.exists(cfgFile):
        print "Configuration file ",cfgFile," does not exist -- Exiting program"
        sys.exit(1)

    Cfg = ConfigParser.ConfigParser()
    Cfg.read(cfgFile)

    self.nevts         = Cfg.getint("InputOutput", "nevts")
    self.outFile       = Cfg.get("InputOutput", "outFile")
    self.inputFiles    = Cfg.get("InputOutput", "inputFiles")
    self.WhichTree     =Cfg.get("InputOutput", "WhichTree")
    self.WhichTrigger  =Cfg.get("InputOutput", "WhichTrigger")
    self.isData        = Cfg.getboolean("InputOutput", "isData")

###
    self.MaxEta     = Cfg.getfloat("Analysis", "MaxEta")
    self.rhoMin     = Cfg.getfloat("Analysis", "rhoMin")
    self.rhoMax     = Cfg.getfloat("Analysis", "rhoMax")
### 

def BookHistograms():

    ### Book the histograms

    hists={}

    hName="Counter"
    hTitle="Event selection counter" 
    hists[hName] = Book1D(hName,hTitle,5,-0.5,4.5,False)

    hName="nHLT"
    hTitle="Number of HLTJets" 
    hists[hName] = Book1D(hName,hTitle,60,-0.5,59.5)

    hName="nvrt"
    hTitle="Number of Primary Vertices" 
    hists[hName] = Book1D(hName,hTitle,40,-0.5,39.5)

    hName="rho"
    hTitle="Event rho" 
    hists[hName] = Book1D(hName,hTitle,100,0,50.)

    hname="drMin"; htitle="Min dr between jets";
    hists[hname] = Book1D(hname, htitle, 100, 0., 1.); 

    hname="drMin_pt40"; htitle="Min dr between jets -- p_{T} > 40 GeV";
    hists[hname] = Book1D(hname, htitle, 100, 0., 1.); 


    njr=20
    jrmin=0.
    jrmax=2.

    neta=20
    etamin=-5.
    etamax=5.

    hname="HLToverReco_vsReco"; htitle="p_{T} (coll2) / p_{T} (coll1) vs jet collection 1 pt,eta";  
    hists[hname]=Book3D(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);
    hname="RecooverHLT_vsReco"; htitle="p_{T} (coll2) / p_{T} (coll1) vs jet collection 2 pt,eta";
    hists[hname]=Book3D(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);



    n=600
    minx=0.
    maxx=600.

    hName="ptHLTJets"
    hTitle="HLT Jet p_{T}" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoJets"
    hTitle="Reco Jet p_{T}" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptHLTL"
    hTitle="Leading HLT Jet p_{T}" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL"
    hTitle="Leading Reco Jet p_{T}" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT25"
    hTitle="Leading Reco Jet p_{T} -- HLTJet25" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT30"
    hTitle="Leading Reco Jet p_{T} -- HLTJet30" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT35"
    hTitle="Leading Reco Jet p_{T} -- HLTJet35" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT40"
    hTitle="Leading Reco Jet p_{T} -- HLTJet40" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT45"
    hTitle="Leading Reco Jet p_{T} -- HLTJet45" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT50"
    hTitle="Leading Reco Jet p_{T} -- HLTJet50" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT55"
    hTitle="Leading Reco Jet p_{T} -- HLTJet55" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT60"
    hTitle="Leading Reco Jet p_{T} -- HLTJet60" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT70"
    hTitle="Leading Reco Jet p_{T} -- HLTJet70" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT80"
    hTitle="Leading Reco Jet p_{T} -- HLTJet80" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT100"
    hTitle="Leading Reco Jet p_{T} -- HLTJet100" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT140"
    hTitle="Leading Reco Jet p_{T} -- HLTJet140" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT200"
    hTitle="Leading Reco Jet p_{T} -- HLTJet200" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT260"
    hTitle="Leading Reco Jet p_{T} -- HLTJet260" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT320"
    hTitle="Leading Reco Jet p_{T} -- HLTJet320" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)

    hName="ptRecoL_HLT400"
    hTitle="Leading Reco Jet p_{T} -- HLTJet400" 
    hists[hName] = Book1D(hName,hTitle,n,minx,maxx)


    return hists

def getNJets(eta0,phi0,ptMin):

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

def deltaPhi(phi1,phi2):
                      
    result=phi1-phi2
    while (result > math.pi): result -= 2*math.pi
    while (result <= -math.pi): result += 2*math.pi

    return result

def deltaR(eta1,phi1,eta2,phi2):
    deta=eta1-eta2
    dphi=deltaPhi(phi1,phi2)
                      
    return math.sqrt(deta*deta + dphi*dphi)

def processEvent(ievt):

    wt=1.
    h["Counter"].Fill(0.)  # count the numbers processed
    # tree.Print()
    

    h["nHLT"].Fill(tree.nhJets,wt)
    # print nPVs
    h["nvrt"].Fill(tree.NVrtx,wt)

    if (tree.rho<jobpar.rhoMin):return
    if (tree.rho>jobpar.rhoMax):return

    h["rho"].Fill(tree.rho,wt)

    if ievt==0:
        print "Number of trigger bits: ",tree.ntrig,ievt
        # hName="triggerResults"
        # hTitle="HLT Results" 
        # h[hName] = Book1D(hName,hTitle,tree.ntrig,0,tree.ntrig,False)

        for key in trigDict.keys():
            # print key,trigDict.get(key)
            h[hName].GetXaxis().SetBinLabel(trigDict.get(key)+1,key);

    for i in xrange( tree.ntrig ):
        if tree.triggerBits[i]==1: h["triggerResults"].Fill(float(i))

    
    ##  Require trigger to fire ##
    itrig=trigDict.get(jobpar.WhichTrigger)
    yesTrig=tree.triggerBits[itrig]
    if yesTrig==0: return
    h["Counter"].Fill(1.)  # count the numbers passing the trigger

    ptHLTL=0.
    for i in xrange( tree.nhJets ):
        h["ptHLTJets"].Fill(tree.hJet_pt[i])
        if tree.hJet_pt[i]>ptHLTL: ptHLTL=tree.hJet_pt[i]

    h["ptHLTL"].Fill(ptHLTL)

    ptRecoL=0.
    n=tree.nrJets

    for i in xrange( tree.nrJets ):
        h["ptRecoJets"].Fill(tree.rJet_pt[i])
        # print i, tree.nrJets,len(tree.rJet_id),tree.nhJets
        goodID=tree.rJet_nhf[i]< 0.99 and \
            tree.rJet_nef[i]< 0.99 and \
            tree.rJet_nconstituents[i]>1
        if (math.fabs(tree.rJet_eta[i])<jobpar.MaxEta)  and goodID:
                if tree.rJet_pt[i]>ptRecoL: ptRecoL=tree.rJet_pt[i]

        if goodID:
            ptr=tree.rJet_pt[i]
            etar=tree.rJet_eta[i]
            phir=tree.rJet_phi[i]
            #find the best matched jet
            drmin=99.;
            match_h=-1;
            for j in xrange( tree.nhJets ):
                etah=tree.hJet_eta[j]
                phih=tree.hJet_phi[j]
                dr=deltaR(etar,phir,etah,phih);                
                if dr<drmin:
                    drmin=dr;
                    match_h=j;
            h["drMin"].Fill(drmin)
            if ptr>40: h["drMin_pt40"].Fill(drmin)
            # plot the response for matched jets
            if drmin<drMatch:
                pth=tree.hJet_pt[match_h]
                h["HLToverReco_vsReco"].Fill(ptr,etar,pth/ptr,1.)
                h["RecooverHLT_vsReco"].Fill(ptr,etar,ptr/pth,1.)

    ## done with double loop
    if ptRecoL>0:
        h["ptRecoL"].Fill(ptRecoL)
        if (ptHLTL>25.): h["ptRecoL_HLT25"].Fill(ptRecoL)
        if (ptHLTL>30.): h["ptRecoL_HLT30"].Fill(ptRecoL)
        if (ptHLTL>35.): h["ptRecoL_HLT35"].Fill(ptRecoL)
        if (ptHLTL>40.): h["ptRecoL_HLT40"].Fill(ptRecoL)
        if (ptHLTL>45.): h["ptRecoL_HLT45"].Fill(ptRecoL)
        if (ptHLTL>50.): h["ptRecoL_HLT50"].Fill(ptRecoL)
        if (ptHLTL>55.): h["ptRecoL_HLT55"].Fill(ptRecoL)
        if (ptHLTL>60.): h["ptRecoL_HLT60"].Fill(ptRecoL)
        if (ptHLTL>70.): h["ptRecoL_HLT70"].Fill(ptRecoL)
        if (ptHLTL>80.): h["ptRecoL_HLT80"].Fill(ptRecoL)
        if (ptHLTL>100.): h["ptRecoL_HLT100"].Fill(ptRecoL)
        if (ptHLTL>140.): h["ptRecoL_HLT140"].Fill(ptRecoL)
        if (ptHLTL>200.): h["ptRecoL_HLT200"].Fill(ptRecoL)
        if (ptHLTL>260.): h["ptRecoL_HLT260"].Fill(ptRecoL)
        if (ptHLTL>320.): h["ptRecoL_HLT320"].Fill(ptRecoL)
        if (ptHLTL>400.): h["ptRecoL_HLT400"].Fill(ptRecoL)

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

    print "Data: ", jobpar.isData
    print "Input:  " ,jobpar.inputFiles
    print "Output: " ,jobpar.outFile


    AutoLibraryLoader.enable()


    # create output root file and book histograms
    outf=OutputRootFile(jobpar.outFile)

    h=BookHistograms()
    hName="triggerResults"
    hTitle="HLT Results" 
    h[hName] = Book1D(hName,hTitle,500,0,500,False)

    tree=getRootChain(jobpar.inputFiles,jobpar.WhichTree)
    f= TFile(jobpar.inputFiles)

    histname=jobpar.WhichTree[:jobpar.WhichTree.find("/")] +"/" + "TriggerResults"
    print "Trigger names taken from:",histname
    triggerHist=f.Get(histname)
    # print triggerHist,triggerHist.GetEntries()
    # print triggerHist.GetXaxis().GetNbins()
    trigDict={}

    for itrig in xrange( triggerHist.GetXaxis().GetNbins() ):
        # print "\t",triggerHist.GetXaxis().GetBinLabel(itrig+1)
        trigDict[triggerHist.GetXaxis().GetBinLabel(itrig+1)]=itrig

    if not trigDict.has_key(jobpar.WhichTrigger):
        print "\n\tRequested trigger name not found -- EXITING"
        sys.exit(1)


    leaves = tree.GetListOfBranches();
    # printLeaves(leaves)

    EVENT  = ROOT.EventInfo()
    tree.SetBranchAddress("EVENT", AddressOf(EVENT, "run") );

    NEntries = tree.GetEntries()
    print "Number of entries on Tree:",NEntries
    nevt=NEntries
    if jobpar.nevts>-1:
        nevt=jobpar.nevts
    print "Processing: ",nevt," events"

    decade=0
    for jentry in xrange( nevt ):

        progress = 10.0*jentry/(1.0*nevt);
        k = math.floor(progress); 
        if (k > decade):
            print 10*k," %",jentry, "\tRun:Event:lumi: ",EVENT.run, EVENT.event, EVENT.lumi
        decade = k  

        tree.GetEntry(jentry) 
        processEvent(jentry)  ### work is done here


    outf.Write();
    outf.Close();
