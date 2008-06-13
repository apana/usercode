#include "rootfuncs.h"

void comp_hists(){

 gROOT->Reset();
 mySetup(0); // 0 -- linear plots, 1 -- log plots  

 TString rootname = "../nhsts/DiJets_Jet150_Scone7_Gen_smrEng_eta3_wght1_hmass.root";

 TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;
 rootfile->GetListOfKeys()->Print();

 TString hname="JetXComp1";
 TH1F *hnum=(TH1F*)GetHist(rootfile,hname); if (!hnum) return;
 TString hname="JetXComp3";
 TH1F *hden=(TH1F*)GetHist(rootfile,hname); if (!hden) return;

 divideByBinWidth(hnum);
 divideByBinWidth(hden);

 TH1F *rat = (TH1F*)hden->Clone(); rat->SetName("eff");
 rat->Divide(hnum,hden,1.,1.,"B");


 Double_t minpt=0,maxpt=60;

 TH1F *h1 = hnum;
 TH1F *h2 = hden;

 h1->SetLineColor(kRed);
 h1->GetXaxis()->SetRangeUser(minpt,maxpt);
 h2->GetXaxis()->SetRangeUser(minpt,maxpt);
 h1->GetXaxis()->SetTitle("p_{T} (GeV)");
 h2->GetXaxis()->SetTitle("p_{T} (GeV)");

 h1->SetTitleSize(.05,"x");

 h1->Draw();
 h2->Draw("same");


 Int_t wtopx,wtopy; UInt_t ww, wh;
 c1->GetCanvasPar(wtopx,wtopy,ww,wh); // Gets Canvas Parameters 
 TCanvas *c2 = new TCanvas("c2","Root Canvas 2",20, wtopy, ww, wh);
 gPad->SetLogy(0);

 Double_t ymin=0,ymax=1.5;
 rat->GetXaxis()->SetRangeUser(minpt,maxpt);
 rat->GetYaxis()->SetRangeUser(ymin,ymax);

 rat->Draw();

 //gROOT->Reset();

}
