#!/usr/bin/env python

import sys,string,math,os
from array import array

from PhysicsTools.PythonAnalysis import *
from ROOT import *
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()

ROOT.gROOT.LoadMacro('higgsStruct_new.h+')

######################################

def CreateTree(filename):

    f = TFile(filename,"recreate");
    tree = TTree( 'treeFR', 'Step 2 friend' )

    tree.Branch( 'HVdPhiFR', HVdPhiFR, 'HVdPhiFR/F' )

    SetOwnership( tree, False )
    return tree

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print 'Usage: ' + sys.argv[0] + ' infile  outfile'
        sys.exit()

    infile=sys.argv[1]
    outfile=sys.argv[2]

    print "Input: ", infile
    print "Output: ",outfile

    f=TFile.Open(infile)
    tree = f.Get("tree");


    NEntries = tree.GetEntries()
    print "Number of entries on Tree:",NEntries
    nevt=NEntries


    HVdPhiFR = array( 'f', [ 0 ] )

    # fr_tree = CreateTree(outfile)
    fout = TFile(outfile,"recreate");
    fr_tree = TTree( 'treeFR', 'Step 2 friend' )
    fr_tree.Branch( 'HVdPhiFR', HVdPhiFR, 'HVdPhiFR/F' )



    V      = ROOT.VInfo()
    EVENT  = ROOT.EventInfo()
    tree.SetBranchAddress("V",     AddressOf(V, "mass") );
    tree.SetBranchAddress("EVENT", AddressOf(EVENT, "run") );

    decade=0
    # nevt=1000
    for jentry in xrange( nevt ):

        progress = 10.0*jentry/(1.0*nevt);
        k = math.floor(progress); 
        if (k > decade):
            print 10*k," %",jentry, "\tRun:Event:lumi:json: ",EVENT.run, EVENT.event, EVENT.lumi, EVENT.json
        decade = k  

        tree.GetEntry(jentry) 
        ## if tree.eventFlav == 5 and V.pt > 100.0 and V.mass>70. and V.mass<110. : ntree.Fill();

        HVdPhiFR[0]=tree.HVdPhi
        fr_tree.Fill();


    # ntree.Print();
    fr_tree.AutoSave();

    # newfile.cd()
    
    f.Close()
    print "Done"


