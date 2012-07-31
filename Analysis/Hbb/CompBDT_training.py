import ROOT
from ROOT import gROOT, gStyle, gSystem, TCanvas, TF1, TFile, TH1F
from ROOT import TColor, TLine, TLegend, TLatex, TObjArray,TPad, TBox, TPaveText
from ROOT import SetOwnership

from ROOT import gDirectory, gPad

import sys,string,math,os
 
ROOT.gROOT.ProcessLine(".x ~/rootmacros/myStyle.cc");

sys.path.append(os.path.join(os.environ.get("HOME"),'rootmacros'))
from myPyRootSettings import prepPlot

# BDTvar="H.mass"
# BDTvar="H.pt"
# BDTvar="V.pt"
BDTvar="hJ12_MaxCsv"
# BDTvar="HVdPhi"
# BDTvar="H.dEta"

savePlot=False
showRatio=True

def DrawText(xtxt,ytxt,leglabel,txtsize=0.045):


    t1 = TLatex();
    t1.SetNDC();

    
    t1.SetTextSize(txtsize);
    t1.DrawLatex(xtxt,ytxt,leglabel);

    return

def DrawLegend(xl1,xl2,yl1,yl2,leglabel=""):

    leg =TLegend(xl1,yl1,xl2,yl2);
    leg.SetFillColor(0);
    leg.SetLineColor(0);
    leg.SetShadowColor(0);
    leg.SetBorderSize(0);
    # leglabel="HLT p_{T} = " + ptbin +" GeV"
    # leglabel="HLT_Jet" + ptbin
    if leglabel != "":
        leg.SetHeader(leglabel);

    leg.SetTextSize(0.04);

    return leg

def DrawPave(x1,y1,x2,y2):

    pave = TPaveText();
    pave.SetX1NDC(x1);
    pave.SetY1NDC(y1);
    pave.SetX2NDC(x2);
    pave.SetY2NDC(y2);
    pave.SetFillColor(ROOT.kAzure-7);
    pave.SetShadowColor(ROOT.kWhite);
    pave.SetTextAlign(12);
    pave.SetTextColor(ROOT.kWhite);
    pave.SetTextFont(42);
    pave.SetTextSize(0.04);
    
    return pave

if __name__ == '__main__':

    gROOT.Reset()
    gROOT.SetStyle("Plain");  
    gStyle.SetOptTitle(0);
    gStyle.SetOptStat(0);


    RootDir="."
    HistFile="TMVA_dijet.root"

    HistFile=os.path.join(RootDir,HistFile)
    print "\n",HistFile
    print "\n\tDistribution:" , BDTvar


    f = TFile(HistFile)

    f.cd("Method_BDT/BDT")
    f.ls()

    hnameSig= "Method_BDT/BDT/" + BDTvar + "__Signal"
    hnameBckg="Method_BDT/BDT/" + BDTvar + "__Background"

    print hnameSig
    hSig=f.Get(hnameSig)
    hBckg=f.Get(hnameBckg)

    xlabel="xxx"
    ylabel="yyy"
    xlegend=0.15
    ylegend=0.75

    SigCol =ROOT.kBlue
    SigFillCol =ROOT.kAzure-8
    BckgCol=ROOT.kRed
    hSig.SetLineWidth(3)
    hSig.SetLineColor(SigCol)
    hSig.SetFillColor(SigFillCol)
    hBckg.SetLineWidth(3)
    hBckg.SetLineColor(BckgCol)
    hBckg.SetFillColor(BckgCol)
    hBckg.SetFillStyle(3004)


    int1=hSig.Integral()
    int2=hBckg.Integral()

    hSig.Scale(1./int1)
    hBckg.Scale(1./int2)

    hSig.GetYaxis().SetTitleOffset(1.5);
    # hNum.GetXaxis().SetTitle(xlabel);
    # hNum.GetYaxis().SetTitle(ylabel);


    c1 = prepPlot("c1","BDT training",550,40,540,550)
    c1.SetTopMargin(0.1)
    # c1.SetLogy(ilogy);    
    hSig.Draw()
    hBckg.Draw("same")

    pave=DrawPave(0.05,0.93,0.6,0.98)
    pave.AddText("Input Variable " + BDTvar);
    pave.Draw()

    xl1=xlegend; yl1=ylegend; xl2=xl1+.3; yl2=yl1+.1;
    leg=DrawLegend(xl1,xl2,yl1,yl2)
    leg.AddEntry(hSig,"Signal","f");
    leg.AddEntry(hBckg,"Background","f");
    leg.Draw()

    if savePlot==1:
        plot="TMVAtraining_"+BDTvar+"_"+compare+".gif"

        print plot
        c1.Print(plot)

    raw_input('\npress return to end the program...')
