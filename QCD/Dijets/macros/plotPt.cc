#include "rootfuncs.h"

void plotPt(){

  mySetup(1); // 0 -- linear plots, 1 -- log plots  

  TString thresh="110";

  TString rootname="../hsts/DijetHistograms_JetET" + thresh + "_Gen_Eff_eta3_unwght_la.root";
  TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;

  TString hname="JetPtMax";
  TH1 *h1=GetHist(rootfile,hname); if (!h1) return;

  double ptmin=0,ptmax=1000.;
  if (thresh == "110"){
    ptmin=150; ptmax=250;
  }else if (thresh == "80"){
    ptmin=120; ptmax=200;
  }else if (thresh == "50"){
    ptmin=70; ptmax=140;
  }else if (thresh == "30"){
    ptmin=40; ptmax=100;
  }else if (thresh == "20"){
    ptmin=25; ptmax=80;
  }
  h1->GetXaxis()->SetRangeUser(ptmin,ptmax);
  h1->Draw();
  
}
