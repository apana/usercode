from ROOT import TF1, TFile, TH1F, TTree, TChain, SetOwnership

def getRootChain(infiles,treeName):

    import glob
    # print infiles
    InputRootFiles = glob.glob (infiles)

    nStep1=0
    tr = TChain(treeName)    
    for rootfile in InputRootFiles:
        print "Adding to chain: ",rootfile,rootfile.find("pnfs")
        if rootfile.find("pnfs")>-1:
            rootfile="dcache:" + rootfile

        print "XXX: ", rootfile

        f=TFile.Open(rootfile)
        h=f.Get('Count')  # get histogram with Step 1 event count
        # print h.GetEntries()
        nStep1 = nStep1 + h.GetEntries()
        f.Close()

        tr.AddFile(rootfile)



    SetOwnership( tr, False )

    print "Number of Step1 events: ", nStep1

    return tr,nStep1

def Book1D(cname,ctitle,nbins,xmin,xmax,doSumw2=True):
    h=TH1F(cname,ctitle,nbins,xmin,xmax)
    if doSumw2: h.Sumw2()
    return h

def OutputRootFile(outfile):

    o = TFile(outfile,"RECREATE");
    print "Output written to: ", outfile
    SetOwnership( o, False )   # tell python not to take ownership

    return o

def printLeaves(leaves):

    leafEnts=leaves.GetSize();
    print "\nNumber of leaves: ", leafEnts
    for i in range(0,leafEnts):
        leafName = leaves[i].GetName();
        print "\t",leafName
