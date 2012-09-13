import ROOT, sys, os, re, string
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2D,TF1, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad
from ROOT import *

from array import array

from ROOT import TGraph
from ROOT import TGraphErrors
from ROOT import TMultiGraph

import numpy

#from tdrStyle import *
#setTDRStyle()

import math
from math import *

#RootFiles = TFile("/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/ZH_125_summer12_33b_120GeV.root")
RootFiles = TFile("/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/En7TeV/ZH_125_summer11_1_Niklas.root")
Tree = RootFiles.Get("tree")

ANALYSIS = "Dijet"
#ANALYSIS = "Subjet"
 

def make_plot(var, bin, low, high, ylabel, xlabel, save, setLog = False):

    _analysis = ANALYSIS

    if _analysis == "Dijet":
        cut = 'H.mass>0. &&  H.mass<200. && V.pt>100. && max(hJet_csv[0],hJet_csv[1])>0.244  && min(hJet_csv[0],hJet_csv[1])>0.244 && hJet_pt[0]>20. && hJet_pt[1]>20. && hJet_eta[0]<2.4 && hJet_eta[1]<2.4 && V.mass>75 && V.mass<105 && fabs(deltaPullAngle) < 2.4'
        cutReg = 'newHiggsMass>0. &&  newHiggsMass<200. && V.pt>100. && max(hJet_csv[0],hJet_csv[1])>0.244  && min(hJet_csv[0],hJet_csv[1])>0.244 &&  hJet_genPtReg0>20. && hJet_genPtReg1>20. && hJet_eta[0]<2.4 && hJet_eta[1]<2.4 && V.mass>75. && V.mass<105. && fabs(deltaPullAngle) < 2.4'

    if _analysis == "Subjet":
        cut = 'FatH.filteredmass>0. &&  FatH.filteredmass<200. && V.pt>250. && max(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.244  && min(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.244 && fathFilterJets_pt[0]>20. && fathFilterJets_pt[1]>20. && fathFilterJets_eta[0]<2.4 && fathFilterJets_eta[1]<2.4 && V.mass>75 && V.mass<105'
        cutReg = 'newfatHiggsMass>0. &&  newfatHiggsMass<200. && V.pt>250. && max(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.244  && min(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.244 &&  fathFilterJets_genPtReg0>20. && fathFilterJets_genPtReg1>20. && fathFilterJets_eta[0]<2.4 && fathFilterJets_eta[1]<2.4 && V.mass>75 && V.mass<105'

    #if _analysis == "Dijet":
    #    cut = 'H.mass>0. &&  H.mass<200. && V.pt>250. && max(hJet_csv[0],hJet_csv[1])>0.9  && min(hJet_csv[0],hJet_csv[1])>0.5 && hJet_pt[0]>40. && hJet_pt[1]>40.'
    #    cutReg = 'newHiggsMass>0. &&  newHiggsMass<200. && V.pt>250. && max(hJet_csv[0],hJet_csv[1])>0.9  && min(hJet_csv[0],hJet_csv[1])>0.5 &&  hJet_genPtReg0>40. && hJet_genPtReg1>40.'
    #
    #if _analysis == "Subjet":
    #    cut = 'FatH.filteredmass>0. &&  FatH.filteredmass<200. && V.pt>250. && max(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.9  && min(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.5 && fathFilterJets_pt[0]>40. && fathFilterJets_pt[1]>40.'
    #    cutReg = 'newfatHiggsMass>0. &&  newfatHiggsMass<200. && V.pt>250. && max(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.9  && min(fathFilterJets_csv[0],fathFilterJets_csv[1])>0.5 &&  fathFilterJets_genPtReg0>40. && fathFilterJets_genPtReg1>40.'
    #

    var1 = var[0]
    var2 = var[1]
    
    histo1 = TH1D(var1, var1, bin, low, high)
    histo2 = TH1D(var2, var2, bin, low, high)
    
    Tree.Draw(var1 + " >> " + var1, cut, 'goff')
    Tree.Draw(var2 + " >> " + var2, cutReg, 'goff')
    
    legend = TLegend(.15,.75,.4,.85)
    c = TCanvas("c","Higgs -> b#b", 1000, 800)
    if setLog:
        c.SetLogy()

    histo1.SetLineColor(1)
    histo1.SetLineWidth(2)
    histo1.SetMaximum(1.2*histo1.GetMaximum())
    histo1.GetXaxis().SetTitle(xlabel)    
    histo2.GetXaxis().SetTitle(xlabel)
    histo1.GetYaxis().CenterTitle()
    histo2.GetYaxis().CenterTitle()
    histo1.GetYaxis().SetTitle(ylabel)
    histo2.GetYaxis().SetTitle(ylabel)
    histo2.SetLineColor(2)
    histo2.SetLineWidth(2)
    histo2.Draw()
    histo1.Draw("SAMES")
    #histo1.Draw()
    #histo2.Draw("SAMES")

    myfit = TF1("myfit","[0]*exp(-0.5*(pow(((x-[1])/[2]),2))) + [3] + [4]*x + [5]*x*x", 90, 140)
    myfit.SetParNames('constant','mean','sigma','const','slope','curve')
    myfit.SetParameter(0, 500)
    myfit.SetParameter(1, 130)
    myfit.SetParameter(2, 16)
    myfit.SetParameter(3, 220)
    myfit.SetParameter(4, -.1)
    myfit.SetParameter(5, -.006)
    #myfit = TF1("myfit", "gaus", 100, 140) 
    histo1.Fit(myfit, "0R")  

    myfit2 = TF1("myfit2","[0]*exp(-0.5*(pow(((x-[1])/[2]),2))) + [3] + [4]*x + [5]*x*x", 90, 140)
    myfit2.SetParNames('constant','mean','sigma','const','slope','curve')
    myfit2.SetParameter(0, 400)
    myfit2.SetParameter(1, 120)
    myfit2.SetParameter(2, 17)
    myfit2.SetParameter(3, 220)
    myfit2.SetParameter(4, -.1)
    myfit2.SetParameter(5, -.006)
    #myfit2 = TF1("myfit2", "gaus", 100, 140) 
    histo2.Fit (myfit2, "0R")

    myfit.SetLineColor(1)
    myfit2.SetLineColor(2)
    myfit2.Draw("SAMES") 
    myfit.Draw("SAMES")

    legend.AddEntry( histo1, var1 , "l")
    legend.AddEntry( histo2, var2 , "l")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw() 

    Incr = (1 - ( (myfit2.GetParameter(2) / myfit2.GetParameter(1)) / (myfit.GetParameter(2) / myfit.GetParameter(1)) ))*100
    print "Increase = ", Incr

    latex=ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
        
    x0=0.12
    y0=0.65
    latex.DrawLatex(x0,y0, '#color[8]{ Improvement = ' + '{0:.0f}'.format(Incr) + '%}')  

    latex.SetTextSize(0.03)
    latex.DrawLatex(.65, 0.5, '#color[1]{ Mean = ' + '{0:.1f}'.format(myfit.GetParameter(1)) + '}')
    latex.DrawLatex(.65, 0.45, '#color[1]{ Sigma = ' + '{0:.2f}'.format(myfit.GetParameter(2)) + '}')
    latex.DrawLatex(.65, 0.4, '#color[2]{ MeanReg = ' + '{0:.1f}'.format(myfit2.GetParameter(1)) + '}')
    latex.DrawLatex(.65, 0.35, '#color[2]{ SigmaReg = ' + '{0:.2f}'.format(myfit2.GetParameter(2)) + '}')

    raw_input('')
    
    #c.SaveAs("Plots/" + save + ".pdf")


_analysis = ANALYSIS
#var = ['H.mass','newHiggsMass']; bin = 70; low = 0; high = 250; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_mass'; setLog = False
if _analysis == "Dijet":
    var = ['H.mass','newHiggsMass']; bin = 80; low = 40; high = 200; ylabel = 'Events / 2 GeV'; xlabel = 'mass'; save = 'Dijet_bJetRegression_mass'; setLog = False
if _analysis == "Subjet":
    var = ['FatH.filteredmass','newfatHiggsMass']; bin = 40; low = 40; high = 200; ylabel = 'Events / 4 GeV'; xlabel = 'Filtered Jet Hmass'; save = 'Subjet_bJetRegression_mass'; setLog = False
    
#var = ['H.pt','HCorr.pt']; bin = 80; low = 0; high = 450; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_pt'; setLog = False
#var = ['hJet_pt[1]','hJet_genPtReg1']; bin = 80; low = 0; high = 450; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_hjetpt0'; setLog = False
#var = ['hJet_pt[1]','hJet_ptCorr[1]']; bin = 50; low = 0; high = 250; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_hjetpt1'; setLog = False
make_plot(var, bin, low, high, ylabel, xlabel, save, setLog)
