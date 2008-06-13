#include "rootfuncs.h"
#include "plotTriggerEff.h"


void plotChiDist(){

  gROOT->Reset();
  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetTitleOffset(1.2,"y");

  bool saveIt=true;
  bool do1pb = true;

  const int nhsts=6;
  TString trignames[nhsts]={"Jet20","Jet30","Jet50","Jet80","Jet110","Jet150"};
  int minMass[nhsts]={380,600, 900,1160,1460,1880};
  int maxMass[nhsts]={600,900,1160,1460,1880,12000};

  TH1F *hists[nhsts];
  TLine *lines[nhsts];
  TLatex *latex[nhsts];

  double ymin=0., ymax=.15;

  TString whichJets;
  //whichJets="_Gen_";
  whichJets="_CaloMCCor_";

  TString whichAlgo = "_Scone7";

  TString xtitle="#chi";
  TString ytitle="1/N dN/d#chi";

  Double_t xl1=.54, yl1=0.59, xl2=xl1+.22, yl2=yl1+.25;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);


  float min=1.,max=20.;

  for (int itrig=0;itrig<nhsts;itrig++){

    //TString rootname="../hsts/DijetHistograms_" + trignames[itrig] + whichJets + "Eff_eta3_wght_yb1.5_ptcut.root";
    //TString rootname="../nhsts/DiJets_" + trignames[itrig] + whichAlgo + whichJets + "eta3_wght1.root";
    TString rootname="../nhsts/DiJets_" + trignames[itrig] + whichAlgo + whichJets + "eta3_wght1";
    if (do1pb) rootname = rootname + "_1pb-1";
    rootname = rootname + ".root";
    TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;

    TString hname="JetXComp3";
    if (itrig>2){
      hname="JetXComp4";
    }

    hists[itrig]=get1DHist(rootfile,"",hname); if (!hists[itrig]) return;

    Int_t istart=hists[itrig]->GetXaxis()->FindBin(min);
    Int_t ifin=hists[itrig]->GetXaxis()->FindBin(max);

    float n=hists[itrig]->Integral(istart,ifin);
    hists[itrig]->Scale(1./n);
    cout << "N= " << n << endl;

    divideByBinWidth(hists[itrig]);

    hists[itrig]->GetXaxis()->SetTitle(xtitle);
    hists[itrig]->GetYaxis()->SetTitle(ytitle);

    latex[itrig] = new TLatex(); latex[itrig]->SetNDC(); latex[itrig]->SetTextAlign(22);  latex[itrig]->SetTextSize(0.05);
  }

  Int_t wtopx=700,wtopy=20; UInt_t ww=540, wh=550;
  int istart=0, ifin=nhsts;
  if (ifin>nhsts) ifin=nhsts;
  for (int itrig=istart;itrig<ifin;itrig++){

    ostringstream canvas("");
    canvas << itrig << endl;
    string cname= "c" + canvas.str();

    TCanvas *c1 = new TCanvas(cname.c_str(),"Root Canvas",wtopx, wtopy, ww, wh);
    wtopy=wtopy+20;

    TH1F *h = hists[itrig]; TLatex *t= latex[itrig];
    h->SetMarkerStyle(20);
    h->SetLineColor(kBlue);
    h->SetMarkerColor(kBlue);

    h->GetXaxis()->SetRangeUser(1,20);
    h->GetYaxis()->SetRangeUser(ymin,ymax);
    
    h->Draw();
    //l->Draw();
    ostringstream ch_min("");
    ostringstream ch_max("");
    ch_min << minMass[itrig];
    ch_max << maxMass[itrig];
    //cout << ch_min.str() << endl;
    //cout << ch_max.str() << endl;
    TString range = ch_min.str() + " < Mass(GeV/c) < " + ch_max.str();
    if ( itrig == nhsts-1 ) range = "Mass(GeV/c) > " + ch_min.str();
    t->DrawLatex(.5,.8,range);
    
    if (saveIt){
      TString psfile="ChiDistribution_" + ch_min.str() + "_M_" + ch_max.str();
      //TString psfile="ChiDistribution";

      if (do1pb) psfile = psfile + "_1pb-1";      
      psfile=psfile + ".eps";
      cout << "Writing file: " << psfile << endl;
      c1->Print(psfile);
    }
  }
}
