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

RootFiles = TFile("/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/ZH_125_summer12_33b.root")
Tree = RootFiles.Get("tree")

def make_plot(var, bin, low, high, ylabel, xlabel, save, setLog = False):
    

    #cut = 'H.pt>120. && H.mass<250. && Vtype==2 &&  V.pt>120. &&  Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )==0  && max(hJet_csv[0],hJet_csv[1])>0.4  && min(hJet_csv[0],hJet_csv[1])>0.4 && hJet_pt[0]>30. && hJet_pt[1]>30.'

    #cut = 'H.mass>0. &&  H.mass<250. && H.pt>120. && ((Vtype==2 && vLepton_pt[0]>20.)) && V.pt>170. &&  Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )==0 && MET.et<1000. && V.pt<1000 && max(hJet_csv[0],hJet_csv[1])>0.4  && min(hJet_csv[0],hJet_csv[1])>0.4 && hJet_csv[0]>0. && hJet_csv[1]>0.&& hJet_pt[0]>30. && hJet_pt[1]>30. && abs(deltaPullAngle)<5.'

    #cutReg = 'newHiggsMass>0. &&  newHiggsMass<250. && H.pt>120. && ((Vtype==2 && vLepton_pt[0]>20.)) && V.pt>170. &&  Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCorrIso < 0.1) )==0 && MET.et<100000. && V.pt<100000 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.4  && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.4 && hJet_csv_nominal[0]>0. && hJet_csv_nominal[1]>0.&& hJet_genPtReg0>30. && hJet_genPtReg1>30. && hbhe && ((EVENT.run<193834 && (triggerFlags[22]>0 || triggerFlags[23]>0)) || (EVENT.run>=193834  && (triggerFlags[14]>0 ||triggerFlags[21]>0))) && EVENT.json>0'


    cut = 'H.mass>0. &&  H.mass<200. && V.pt>150. && max(hJet_csv[0],hJet_csv[1])>0.4  && min(hJet_csv[0],hJet_csv[1])>0.4 && hJet_pt[0]>30. && hJet_pt[1]>30.'
    cutReg = 'newHiggsMass>0. &&  newHiggsMass<200. && V.pt>150. && max(hJet_csv[0],hJet_csv[1])>0.4  && min(hJet_csv[0],hJet_csv[1])>0.4 &&  hJet_genPtReg0>30. && hJet_genPtReg1>30.'


    #cut = 'H.pt>120. && Vtype==2 &&  V.pt>120. &&  Sum$(aLepton_pt > 15 && abs(aLepton_eta) < 2.5 && (aLepton_pfCombRelIso < 0.15) )==0  && max(hJet_csv[0],hJet_csv[1])>0.898  && min(hJet_csv[0],hJet_csv[1])>0.4 && hJet_pt[0]>30. && hJet_pt[1]>30. && Sum$(aJet_pt>20 && abs(aJet_eta)<4.5)==0 && HVdPhi > 2.95'

    var1 = var[0]
    var2 = var[1]
    
    histo1 = TH1D(var1, var1, bin, low, high)
    histo2 = TH1D(var2, var2, bin, low, high)
    
    Tree.Draw(var1 + " >> " + var1, cut, 'goff')
    Tree.Draw(var2 + " >> " + var2, cutReg, 'goff')
    
    legend = TLegend(.75,.60,.90,.90)
    c = TCanvas("c","Higgs -> b#b", 1000, 800)
    if setLog:
        c.SetLogy()

    histo1.SetLineColor(1)
    histo1.SetLineWidth(2)
    histo1.SetMaximum(1.2*histo1.GetMaximum())
    histo1.GetXaxis().SetTitle(xlabel)
    histo1.GetYaxis().CenterTitle()
    histo1.GetYaxis().SetTitle(ylabel)
    histo1.Draw()
    histo2.SetLineColor(2)
    histo2.SetLineWidth(2)
    histo2.Draw("SAMES")


    myfit = TF1("myfit","gaus", 100, 150 )
    histo1 . Fit (myfit, "0")
    myfit.SetLineColor(1)
    myfit.Draw("SAMES")

    myfit2 = TF1("myfit2","gaus", 100, 150 )
    histo2 . Fit (myfit2, "0")
    myfit2.SetLineColor(2)
    myfit2.Draw("SAMES")

    print '% Increase = ', (1 - ( (myfit2.GetParameter(2) / myfit2.GetParameter(1)) / (myfit.GetParameter(2) / myfit.GetParameter(1)) ))*100

    legend . AddEntry( histo1, var1 , "l")
    legend . AddEntry( histo2, var2 , "l")

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw()     

    raw_input('')
    
    #c.SaveAs("Plots/" + save + ".pdf")


#var = ['H.mass','newHiggsMass']; bin = 70; low = 0; high = 250; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_mass'; setLog = False
var = ['H.mass','newHiggsMass']; bin = 40; low = 40; high = 200; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_mass'; setLog = False
#var = ['H.pt','HCorr.pt']; bin = 80; low = 0; high = 450; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_pt'; setLog = False
#var = ['hJet_pt[1]','hJet_genPtReg1']; bin = 80; low = 0; high = 450; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_hjetpt0'; setLog = False
#var = ['hJet_pt[1]','hJet_ptCorr[1]']; bin = 50; low = 0; high = 250; ylabel = 'Events / 5 GeV'; xlabel = 'mass'; save = 'bJetRegression_hjetpt1'; setLog = False
make_plot(var, bin, low, high, ylabel, xlabel, save, setLog)
