import sys, os
import ROOT
from ROOT import gROOT, gStyle, TFile, TH1, TH1F

from zmumu_samples_new import SAMPLES

sys.path.append(os.path.join(os.environ.get("HOME"),'rootmacros'))
from myPyRootSettings import prepPlot, prep1by2Plot

# scaling luminosity
lumi = 5.0 # /fb

savePlot=False
plotRatio=True

mlow=91.
mhigh=139.

YMIN=0.001
YMAX=1000.

XMIN=0.
XMAX=300.

LOGY=False

# mlow=0.
# mhigh=300.

# histogram attributes -----------------------------------------------------------------------------
class HistAttr2:
    def __init__(self, fileName_, histName_, label_, N_, xSection_, color_, flag_,addLeg_=True):
        self.fileName = fileName_
        self.histName = histName_
        self.label = label_
        self.N = N_
        self.xSection = xSection_
        self.color = color_
        self.flag = flag_
        self.addLeg = addLeg_

class HistAttr:
    def __init__(self, sampleInfo_, histName_, leglabel_, color_,flag_,addLeg_=True):

        self.fileName = sampleInfo_[0]
        #self.N = sampleInfo_[1]
        self.N = self.GetStep1Entries()
        self.xSection = sampleInfo_[2]

        self.histName = histName_
        self.label = leglabel_
        self.color = color_
        self.flag = flag_
        self.addLeg = addLeg_
    def printIt(self):
        print "Filename: ",self.fileName
        print "\tnumber of events in Step 1: ",self.N
        print "\tCross section: ",self.xSection
        print "\tHistogram name: ",self.histName," color: ",self.color," flag: ",self.flag

    def GetStep1Entries(self):

        f=TFile.Open(self.fileName)
        theHist = f.Get("NStep1")
        n=int(theHist.GetBinContent(1))
        f.Close()

        return n
        
        


# plot histograms ----------------------------------------------------------------------------------
class PlotHist:
    
    datHistAttrs = []
    sgnHistAttrs = []
    bkgHistAttrs = []
    rootFiles = []
    datHists = []
    sgnHists = []
    bkgHists = []
    
    # initialization
    def __init__(self, name_, title_, rebinFactor_):
        self.name = name_
        self.title = title_
        self.rebinFactor = rebinFactor_
    
    # add histogram
    def AddHist(self, sampleID,histName_, leglabel_, color_,flag_,addLeg_=True):
        if SAMPLES.has_key(sampleID):
            histAttr=HistAttr(SAMPLES[sampleID],histName_, leglabel_, color_,flag_,addLeg_)
            #histAttr.printIt()

            if histAttr.flag == 'D':
                self.datHistAttrs.append(histAttr)
            elif histAttr.flag == 'S':
                self.sgnHistAttrs.append(histAttr)
            elif histAttr.flag in ['b', 'B']:
                self.bkgHistAttrs.append(histAttr)
        else:
            print "Unknown sample identifier"
            sys.exit(1)

    # setup canvas
    def SetCanvas(self):
        ROOT.gROOT.SetStyle('Plain')
        self.canvas = ROOT.TCanvas('canvas_' + self.name)
        self.canvas.SetLogy(LOGY)
        self.canvas.SetTopMargin(0.05)
        self.canvas.SetRightMargin(0.05)
        self.canvas.SetTicks()
    
    # setup legend
    def SetLegend(self):
        self.legend = ROOT.TLegend(0.62, 0.88, 0.82, 0.88)
        if plotRatio:
            self.legend.SetTextSize(0.04)
        else:
            self.legend.SetTextSize(0.03)
        self.legend.SetFillColor(ROOT.kWhite)
        self.legend.SetLineColor(ROOT.kWhite)
        self.legend.SetShadowColor(ROOT.kWhite);
    
    # increase legend size on Y-axis
    def IncLegdYsize(self):
        self.legend.SetY1(self.legend.GetY1() - 0.04)
    
    # setup latex
    def SetLatex(self):
        self.latex = ROOT.TLatex()
        self.latex.SetNDC()
        self.latex.SetTextSize(0.04)
    
    # setup histogram
    def SetHist(self, histAttr):
        self.rootFiles.append(ROOT.TFile(histAttr.fileName))
        hist = self.rootFiles[-1].Get(histAttr.histName)
        hist.SetName(histAttr.histName + '_' + histAttr.label)
        hist.SetTitle("")
        hist.Rebin(self.rebinFactor)
        if histAttr.flag != 'D':
            hist.Scale(lumi * histAttr.xSection / histAttr.N)
        if histAttr.flag == 'D':
            hist.SetLineColor(histAttr.color)
            hist.SetMarkerColor(histAttr.color)
            hist.SetMarkerStyle(20)
            hist.SetMarkerSize(1.2)
            self.datHists.append(hist)
        elif histAttr.flag == 'S':
            hist.SetLineWidth(2)
            hist.SetLineColor(histAttr.color)
            self.sgnHists.append(hist)
        elif histAttr.flag in ['b', 'B']:
            hist.SetLineWidth(1)
            if histAttr.flag == 'B':
                hist.SetLineColor(ROOT.kBlack)
            else:
                hist.SetLineColor(histAttr.color)
            hist.SetFillColor(histAttr.color)
            hist.SetFillStyle(1001)
            self.bkgHists.append(hist)
    
    def Ratio(self):
        #print len(self.datHists)
        hData=self.datHists[0]

        nbins=hData.GetNbinsX()
        x1=hData.GetBinLowEdge(1)
        x2=hData.GetBinLowEdge(nbins)+hData.GetBinWidth(nbins)

        hMC=self.histStack.GetStack().Last()

        self.hRat= hData.Clone()
        self.hRat.SetName("Ratio")
        self.hRat.GetXaxis().SetRangeUser(XMIN,XMAX)

        self.hRat.Divide(hData,hMC,1.,1.,"");
        self.hRat.SetMinimum(0.)
        self.hRat.SetMaximum(2.49)
        
        pfit = ROOT.TF1("P0", 'pol0',x1,x2)
        pfit.SetLineWidth (2)
        pfit.SetLineStyle (2)
        self.hRat.Fit(pfit,"0QR")

    # print results
        par = pfit.GetParameters()
        err = pfit.GetParErrors()
        chisqNDF=pfit.GetChisquare()/pfit.GetNDF()
        print 'fit results: const =', par[0], ', +- ', err[0]
        lsize=0.08
        
        self.hRat.GetYaxis().SetTitle("Data/MC");
        self.hRat.GetXaxis().SetTitle(Analysis);

        self.hRat.GetXaxis().SetLabelSize(lsize);
        self.hRat.GetYaxis().SetLabelSize(lsize);

        self.hRat.GetYaxis().SetTitleSize(.09);
        self.hRat.GetXaxis().SetTitleSize(.09);
        self.hRat.GetYaxis().SetTitleOffset(.5);
        self.hRat.GetXaxis().SetTitleOffset(1.1);

        
        self.hRat.GetYaxis().SetNdivisions(505);
        kNotDraw = 1<<9;
        self.hRat.GetFunction("P0").ResetBit(kNotDraw)
        self.hRat.Draw("")

        latex=ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.1)
        
        x0=0.15
        y0=0.85
        latex.DrawLatex(x0,y0, 'P0 = ' + '{0:.2f}'.format(par[0]) + ' #pm ' + '{0:.2f}'.format(err[0]))
        y0=y0-0.1
        latex.DrawLatex(x0,y0, '#chi^{2}_{#nu} = ' + '{0:.2f}'.format(chisqNDF))


    # draw histograms
    def Draw(self):
        #self.SetCanvas()
        self.SetLegend()
        self.SetLatex()
        self.histStack = ROOT.THStack('histStack_' + self.name, '')
        for histAttr in self.bkgHistAttrs:
            self.SetHist(histAttr)
            self.histStack.Add(self.bkgHists[-1])
            if histAttr.flag == 'B':
                if histAttr.addLeg:
                    self.IncLegdYsize()

                imin=self.bkgHists[-1].GetXaxis().FindBin(mlow)
                imax=self.bkgHists[-1].GetXaxis().FindBin(mhigh)
                error= ROOT.Double(0.)
                int= self.bkgHists[-1].IntegralAndError(imin,imax,error)
                label=histAttr.label

                if histAttr.addLeg:
                    self.legend.AddEntry(self.bkgHists[-1], label, 'f')
 
                self.hxxx=self.bkgHists[-1].Clone(self.name + '_xxx')


        self.histTemp = self.bkgHists[-1].Clone(self.name + '_temp')
        self.histTemp.Reset()
        self.histTemp.SetStats(0)
        self.histTemp.GetXaxis().SetTitle(self.title)
        bw=self.histTemp.GetBinWidth(0)
        self.histTemp.GetYaxis().SetTitle("Events"+'/ {0:.2f}'.format(bw))
        self.histTemp.GetYaxis().SetRangeUser(YMIN,YMAX)

        if (plotRatio):
            self.histTemp.GetXaxis().SetLabelOffset(100.);
            self.histTemp.GetYaxis().SetTitleOffset(1.2);

        self.histTemp.GetXaxis().SetRangeUser(XMIN,XMAX)
        self.histTemp.Draw()
        self.histStack.Draw('same hist')
        
        for histAttr in self.sgnHistAttrs:
            self.SetHist(histAttr)
            self.IncLegdYsize()

            imin=self.sgnHists[-1].GetXaxis().FindBin(mlow)
            imax=self.sgnHists[-1].GetXaxis().FindBin(mhigh)
            print imin,imax
            error= ROOT.Double(0.)
            int= self.sgnHists[-1].IntegralAndError(imin,imax,error)

            label=histAttr.label
            self.legend.AddEntry(self.sgnHists[-1], label, 'l')
            self.sgnHists[-1].Draw('same')


        for histAttr in self.datHistAttrs:
            self.SetHist(histAttr)
            self.IncLegdYsize()
            self.legend.AddEntry(self.datHists[-1], histAttr.label, 'p')
            self.datHists[-1].Draw('same x0 e p')
        self.histTemp.Draw('same axis')
        self.legend.Draw()

        ltx = ROOT.TLatex()
        ltx.SetNDC()
        ltx.SetTextSize(0.05)
        ltx.DrawLatex(0.15,0.89, 'CMS')
        
        yl=0.84
        xl=0.15
        self.latex.DrawLatex(xl,yl, '#sqrt{s} = 8 TeV')
        yl=yl-0.05
        self.latex.DrawLatex(xl,yl, 'L = ' + str(lumi)+ ' fb^{-1}')
        yl=yl-0.05
        self.latex.DrawLatex(xl,yl, 'Z(mu^{+}mu^{-})H(b#bar{b})')

        #self.canvas.Update()

        if savePlot == True:
            plot=histname[histname.find("/")+1:] +"_udscgContrl.gif"
            print plot
            self.canvas.Print(plot)

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: ' + sys.argv[0] + ' histName histTitle rebinFactor'
        sys.exit()


    gROOT.Reset()
    gROOT.SetStyle("Plain");    
    # gStyle.SetOptLogy(0);
    gStyle.SetPalette(1);
    gStyle.SetOptTitle(0);
    gStyle.SetOptStat(0);
    gStyle.SetPadTickX(1);

    histname=sys.argv[1]
    Analysis=sys.argv[2]
    CR=sys.argv[3]

    if histname.find("ttbarCR")>-1:
        if histname.find("Hmass_")>-1:
            YMAX=60
        elif histname.find("Hpt_")>-1:
            YMAX=140
        elif histname.find("csv_jet1")>-1:
            YMAX=300
        elif histname.find("csv_jet2")>-1:
            YMAX=160
        elif histname.find("Mll")>-1:
            YMAX=100
        elif histname.find("Zll_pT")>-1:
            YMAX=140
        elif histname.find("BDT_")>-1:
            YMAX=170
            XMIN=-1
    elif histname.find("bbCR")>-1:
        if histname.find("Hmass_")>-1:
            YMAX=50
        elif histname.find("Hpt_")>-1:
            YMAX=80
        elif histname.find("csv_jet1")>-1:
            YMAX=300
        elif histname.find("csv_jet2")>-1:
            YMAX=100
        elif histname.find("Mll_pT")>-1:
            YMAX=100
        elif histname.find("Zll_pT")>-1:
            YMAX=100
        elif histname.find("BDT_")>-1:
            YMAX=100
            XMIN=-1
    elif histname.find("udscgCR")>-1:
        if histname.find("Hmass_")>-1:
            YMAX=2e7
            LOGY=True
        elif histname.find("Hpt_")>-1:
            YMAX=1300
        elif histname.find("csv_jet1")>-1:
            YMAX=1.8e3
        elif histname.find("csv_jet2")>-1:
            YMAX=4.2e3
        elif histname.find("Mll_pT")>-1:
            YMAX=2500
        elif histname.find("Zll_pT")>-1:
            YMAX=1.8e3
        elif histname.find("BDT_")>-1:
            YMAX=2e7
            YMIN=1e-2
            XMIN=-1
            LOGY=True


    plotHist = PlotHist(sys.argv[1], sys.argv[2], int(sys.argv[3]))

    analysis="Zmumu"

    histName=sys.argv[1]

    # SingleTop
    plotHist.AddHist("ST-s",     histName, 'SingleTop', ROOT.kCyan, 'b')
    plotHist.AddHist("STbar-s",  histName, 'SingleTop', ROOT.kCyan, 'b')
    ### NA  plotHist.AddHist("ST-t",  histName, 'SingleTop', ROOT.kCyan, 'b')
    plotHist.AddHist("STbar-t",  histName, 'SingleTop', ROOT.kCyan, 'b')
    plotHist.AddHist("ST-tW",    histName, 'SingleTop', ROOT.kCyan, 'b')
    plotHist.AddHist("STbar-tW", histName, 'SingleTop', ROOT.kCyan, 'B')
    
    ## # TTbar
    plotHist.AddHist("TTbar", histName, 'TTbar', ROOT.kBlue, 'B')

    darkYellow=ROOT.TColor.GetColorDark(ROOT.kYellow);
    plotHist.AddHist("ZJets_udscg_ptL",histName,'Z + udscg',ROOT.kOrange-2, 'b')
    plotHist.AddHist("ZJets_udscg_ptH",histName,'Z + udscg',ROOT.kOrange-2, 'B')

    plotHist.AddHist("ZJets_bJets_ptL", histName, 'Z + bb', ROOT.kYellow, 'b')
    plotHist.AddHist("ZJets_bJets_ptH", histName, 'Z + bb', ROOT.kYellow, 'B')


    # DiBosons
    plotHist.AddHist("ZZ", histName, 'ZZ', ROOT.kGray+1, 'b')
    plotHist.AddHist("WW", histName, 'WW', ROOT.kGray+1, 'b')
    plotHist.AddHist("WZ", histName, 'VV (ZZ/WW/WZ)', ROOT.kGray+1, 'B')

    ## # Higgs
    plotHist.AddHist("HIGGS", histName, 'ZH (125)', ROOT.kRed, 'S')


    ## # Data
    plotHist.AddHist("DATA", histName, 'Data', ROOT.kBlack, 'D')
    
    if not plotRatio:
        c1=prepPlot("c1","Comp",600,20,600,600)
        c1.SetLogy(LOGY)
    else:
        c1,pad1,pad2 = prep1by2Plot("c1","Comp",600,20,600,600)
        print LOGY
        pad1.SetLogy(LOGY)
        pad1.Draw()
        pad2.Draw()
        pad1.cd()
    
    plotHist.Draw()

    if plotRatio:
        pad2.cd()
        plotHist.Ratio()

    c1.Update()

    raw_input('\npress return to end the program...')
