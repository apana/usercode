#include "rootfuncs.h"
#include "plotTriggerEff.h"


void plotTriggerEff(){

  gROOT->Reset();
  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  //gStyle->SetMarkerSize(.8);
  gStyle->SetOptFit(1);
  bool saveIt=true;

  const int neffs=5;
  TString trignames[neffs]={"JetET30","JetET50","JetET80","JetET110","JetET150"};
  Int_t trigcols[neffs]={kRed,38,kDarkGreen,kMagenta,kBlue};
  Int_t trigmrks[neffs]={20,20,20,20,20};
  double mintrg[neffs]={ 320, 500, 750,1000,1240};
  double maxtrg[neffs]={1000,1200,1450,1700,2300};
  double xarr[neffs]={600,900,1160,1460,1770};
  TH1F *hists[neffs];
  TLine *lines[neffs];
  TLatex *latex[neffs];
  TLine *arrow[neffs];

  double ymin=0., ymax=1.49;
  //if (combineTrgs) ymax=1e8;

  TString whichJets;
  //whichJets="_Gen_";
  whichJets="_CaloMCCor_";

  TString xtitle="Dijet Mass (GeV)";
  TString ytitle="Efficiency";

  Double_t xl1=.54, yl1=0.59, xl2=xl1+.22, yl2=yl1+.25;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);

  TString rootname="../eff/TriggerEfficiency_JetMass.root";
  TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;

  for (int ieff=0;ieff<neffs;ieff++){
    ostringstream ch_i("");
    ch_i << ieff+1;
    TString hname="Efficiency_sample" + ch_i.str() + "_Y0";

    hists[ieff]=get1DHist(rootfile,"",hname); if (!hists[ieff]) return;
    hists[ieff]->GetXaxis()->SetRangeUser(mintrg[ieff],maxtrg[ieff]);
    hists[ieff]->GetYaxis()->SetRangeUser(ymin,ymax);
    hists[ieff]->GetXaxis()->SetTitle(xtitle);
    hists[ieff]->GetYaxis()->SetTitle(ytitle);

    lines[ieff] = new TLine(mintrg[ieff],1.,maxtrg[ieff],1.); lines[ieff]->SetLineStyle(2);
    //cuts[ieff] = new TLine(mintrg[ieff],1.,maxtrg[ieff],1.); lines[ieff]->SetLineStyle(2);
    latex[ieff] = new TLatex(); latex[ieff]->SetNDC(); latex[ieff]->SetTextAlign(12);  latex[ieff]->SetTextSize(0.05);

    //Float_t yarr1=1.3, yarr2=1.02, arrsiz=0.03;
    Float_t yarr1=0.7, yarr2=0.95, arrsiz=0.03;
    if (ieff == 1){yarr1=1.3; yarr2=1.02;}
    arrow[ieff] = new TArrow(xarr[ieff],yarr1,xarr[ieff],yarr2,arrsiz,"->");
    arrow[ieff]->SetLineStyle(1); arrow[ieff]->SetLineWidth(2); arrow[ieff]->SetLineColor(kDarkGreen);
  }

  int n=neffs;
  //n=1;
  for (int ieff=0;ieff<n;ieff++){
    //    int ieff=4;
    TH1F *h = hists[ieff]; TLine *l = lines[ieff];  TLatex *t= latex[ieff]; TLine *arr= arrow[ieff];
    h->SetMarkerStyle(20);
    h->SetLineColor(kBlue);
    //TF1* f= h->GetFunction("FitEfficiency_sample0_Y0");
    //f->SetLineColor(kRed);
    TList *mylist= h->GetListOfFunctions();
    
    cout << mylist->First() << endl;
    TF1* f= mylist->First();
    f->SetLineColor(kRed);
    f->SetLineWidth(3);
    h->SetMarkerColor(kBlue);
    
    h->Draw();
    l->Draw();
    t->DrawLatex(.15,.85,trignames[ieff]);
    arr->Draw();
    
    if (saveIt){
      TString psfile="TriggerEff_" + trignames[ieff];      
      psfile=psfile + ".eps";
      cout << "Writing file: " << psfile << endl;
      c1->Print(psfile);
    }
  }
}
