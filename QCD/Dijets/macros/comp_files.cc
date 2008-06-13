#include "rootfuncs.h"

void comp_files(){

 gROOT->Reset();
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 //gStyle->SetTitleOffset(0.8,"y");
 //gStyle->SetOptLogy(0);

 //gStyle->SetTitleSize(0.05,"y");
 //gStyle->SetTitleSize(0.05,"x");

 //gStyle->SetLabelSize(0.04,"y");
 //gStyle->SetLabelSize(0.045,"x");

 //gStyle->SetHistFillColor(kYellow);
 gStyle->SetOptLogy(1);
 gROOT->ForceStyle();


 TCanvas *c1 = new TCanvas();
 //TString rootname2 = "./hlt_hists_CSA07_QCD_167.root";
 //TString rootname1 = "./hlt_hists_sbQCD_167.root";

 TString rdir="../nhsts/";
 TString rootname2 = rdir + "DiJets_Jet80_Scone7_CaloMCCor_engShift+0.1_eta3_wght1_extrap.root";
 TString rootname1 = rdir + "DiJets_Jet80_Scone7_CaloMCCor_engShift+0.1_eta3_wght1.root";

 TString rdir="../nhsts/";
 TString rootname2 = "../DiJets_Jet20_Scone7_Gen_eta3_wght1_v2.root";
 TString rootname1 = "../nhsts/DiJets_Jet20_Scone7_Gen_eta3_wght1.root";


 TFile *rootfile1=OpenRootFile(rootname1); if (!rootfile1) return;
 rootfile1->GetListOfKeys()->Print();
 TFile *rootfile2=OpenRootFile(rootname2); if (!rootfile2) return;

 TString hname="JetMass";
 //TString hname="JetPtMax";
 TH1F *h1=(TH1F*)GetHist(rootfile1,hname); if (!h1) return;
 TH1F *h2=(TH1F*)GetHist(rootfile2,hname); if (!h2) return;

 divideByBinWidth(h1);
 divideByBinWidth(h2);

 int rebin=5;
 h1->Rebin(rebin);
 h2->Rebin(rebin);

 TH1F *rat = (TH1F*)h1->Clone(); rat->SetName("ratio");
 rat->Divide(h2,h1,1.,1.,"");

 Double_t minpt=0,maxpt=3000;

 h1->SetLineColor(kRed);
 h1->GetXaxis()->SetRangeUser(minpt,maxpt);
 h2->GetXaxis()->SetRangeUser(minpt,maxpt);
 h1->GetXaxis()->SetTitle("p_{T} (GeV)");
 h2->GetXaxis()->SetTitle("p_{T} (GeV)");

 h1->SetTitleSize(.05,"x");

 h2->Draw();
 h1->Draw("same");


 Int_t wtopx,wtopy; UInt_t ww, wh;
 c1->GetCanvasPar(wtopx,wtopy,ww,wh); // Gets Canvas Parameters 
 TCanvas *c2 = new TCanvas("c2","Root Canvas 2",20, wtopy, ww, wh);
 gPad->SetLogy(0);
 
 Double_t ymin=0,ymax=2.1;
 rat->GetXaxis()->SetRangeUser(minpt,maxpt);
 rat->GetYaxis()->SetRangeUser(ymin,ymax);

 rat->Draw();


 //gROOT->Reset();

}
