#include "rootfuncs.h"

void compSpectrum(){

  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetOptFit(111);

  TString thresh="150";

  //TString rootname1="../nhsts/DiJets_Jet" + thresh + "_Scone7_Gen_smrEng_adjMass_eta3_wght1.root";
  //TString rootname2="../nhsts/DiJets_Jet" + thresh + "_Scone7_Gen_adjMass_eta3_wght1.root";

  TString rootname1="../nhsts/DiJets_Jet" + thresh + "_Scone7_Gen_smrEng_eta3_wght1_noEtaCut_v2.root";
  TString rootname2="../nhsts/DiJets_Jet" + thresh + "_Scone7_Gen_eta3_wght1_noEtaCut_v2.root";

  TFile *rootfile1=OpenRootFile(rootname1); if (!rootfile1) return;
  TFile *rootfile2=OpenRootFile(rootname2); if (!rootfile2) return;

  //TString hname="JetPtMax";
  TString hname="JetXComp3";
  TH1 *h1=GetHist(rootfile1,hname); if (!h1) return;
  TH1 *h2=GetHist(rootfile2,hname); if (!h2) return;

  //Int_t rebin=1;
  //h1->Rebin(rebin);
  //h2->Rebin(rebin);

  divideByBinWidth(h1);
  divideByBinWidth(h2);

  double min=1.,max=20.;

  Int_t istart=h1->GetXaxis()->FindBin(min);
  Int_t ifin=h1->GetXaxis()->FindBin(max);

  cout << "Integrating bins from: " << istart << " - " << ifin << endl;

  float n1=h1->Integral(istart,ifin); h1->Scale(1./n1);
  float n2=h2->Integral(istart,ifin); h2->Scale(1./n2);

  scaleErrors(h2,0.);
  TH1F *rat = (TH1F*)h1->Clone(); rat->SetName("ratio");
  rat->Divide(h2,h1,1.,1.,"");

  h1->GetXaxis()->SetRangeUser(min,max);
  h2->SetLineColor(kRed);
  h1->Draw();
  h2->Draw("same");
  
 Int_t wtopx,wtopy; UInt_t ww, wh;
 c1->GetCanvasPar(wtopx,wtopy,ww,wh); // Gets Canvas Parameters
 TCanvas *c2 = new TCanvas("c2","Root Canvas 2",20, wtopy, ww, wh);
 gPad->SetLogy(0);

 Double_t ymin=0,ymax=1.5;
 rat->GetXaxis()->SetRangeUser(min,max);
 rat->GetYaxis()->SetRangeUser(ymin,ymax);
 //TF1 *p1fit = new TF1("P1","pol1",min,max);
 //rat->Fit(p1fit,"RQ");

 TF1 *p0fit = new TF1("P0","pol0",min,max);
 rat->Fit(p0fit,"RQ");

 rat->Draw();

}
