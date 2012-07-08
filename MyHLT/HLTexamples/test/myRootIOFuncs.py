from ROOT import TF1, TFile, TH1F, TH3F, TTree, TChain, SetOwnership

def getRootChain(infiles,treeName):

    import glob
    print infiles
    InputRootFiles = glob.glob (infiles)

    tr = TChain(treeName)    
    for rootfile in InputRootFiles:
        print "Adding to chain: ",rootfile
        if rootfile.find("pnfs")>-1:
            rootfile="dcache:" + rootfile
        print rootfile
        tr.AddFile(rootfile)

    SetOwnership( tr, False )
    return tr

def Book1D(cname,ctitle,nbins,xmin,xmax,doSumw2=True):
    h=TH1F(cname,ctitle,nbins,xmin,xmax)
    if doSumw2: h.Sumw2()
    return h

def Book3D(cname,ctitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax,nbinsz,zmin,zmax,doSumw2=True):
    h=TH3F(cname,ctitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax,nbinsz,zmin,zmax)
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
