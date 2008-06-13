#include "rootfuncs.h"

void plotPtSpectrum(){

  mySetup(1); // 0 -- linear plots, 1 -- log plots  

  const int ntrigs=6;
  TString trignames[ntrigs]={"JetET20","JetET30","JetET50","JetET80","JetET110","JetET150"};
  TH1 *hists[ntrigs];

  for (int itrig=0;itrig<ntrigs;itrig++){

    TString rootname="../hsts/DijetHistograms_" + trignames[itrig] + "_CaloMCCor_Eff_eta3_norm_la.root";
    TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;

    TString hname="JetPtMax";
    hists[itrig]=GetHist(rootfile,hname); if (!hists[itrig]) return;

    TH1 *h = hists[itrig];
    TString chopt="e,same";
    if (itrig ==0){
      double ptmin=0,ptmax=500.;
      h->GetXaxis()->SetRangeUser(ptmin,ptmax);
      chopt="hist";
    }
    h->Draw(chopt);
  }

  return;

  double ptmin=0,ptmax=400.;
  h1->GetXaxis()->SetRangeUser(ptmin,ptmax);
  h1->Draw("e,hist");
  h2->Draw("e,same");
  
  h1->Sumw2(); h2->Sumw2(); // store sum of squares of weights (if not already done)
  TH1F *h12 = h1->Clone(); h12->SetName("Sum"); // Clone one of the histograms
  h12->Add(h1,h2,1.,1.);

  //h12->Draw("same");
  //h12->Draw("hist");
}
