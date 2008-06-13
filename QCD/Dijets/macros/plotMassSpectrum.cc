#include "rootfuncs.h"
// #include "d0xs.h"

Double_t fitfunc(Double_t *x, Double_t *par)
{

  //double value=1;
  //float N=1.e13;
  //float a=-4;
  //float a=-5;
  //float b=1.;

  double pt=x[0];

  float N=par[0];
  float a=par[1];
  float b=par[2];

  double roots=par[3];


  double value=N*pow(pt,a)*pow((1-2*pt/roots),b);

  //double value=N*pow(pt,a);
  //double value=pow((1-2*pt/roots),b);
  return value;
}

void plotMassSpectrum(){

  gROOT->Reset();
  mySetup(1); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetMarkerSize(.8);
  gStyle->SetLineWidth(1);

  bool do1pb = true;
  bool combineTrgs=true;
  bool compD0=true;
  bool saveIt=true;
  saveIt=false;

  double xmin=0, xmax=3000;
  double ymin=2e1, ymax=3e7;

  if (do1pb){ ymin=2e-3; ymax=9e3;}

  if (compD0){ xmin=200, xmax=1000; ymin=2e-1;}

  //double xmin=2000, xmax=4000;
  //double ymin=2e-2, ymax=1.5e1;

  //if (combineTrgs) ymax=1e8;
  const int ntrigs=6;
  //TString trignames[ntrigs]={"JetET20","JetET30","JetET50","JetET80","JetET110","JetET150"};
  TString trignames[ntrigs]={"Jet20","Jet30","Jet50","Jet80","Jet110","Jet150"};
  Int_t trigcols[ntrigs]={kRed,38,kDarkGreen,kMagenta,kBlue,48};
  Int_t trigmrks[ntrigs]={20,20,20,20,20,20};
  double mintrg[ntrigs]={  0,150,250,350,450,550};
  double maxtrg[ntrigs]={660,1000,1500,2000,2500,3000};

  Int_t combcol=kDarkGreen;
  Int_t combmrk=20;

  //TString trignames[ntrigs]={"JetET20","JetET30"}; xmin=400; xmax=800;
  //TString trignames[ntrigs]={"JetET30","JetET50"};  xmin=650; xmax=1000;
  //TString trignames[ntrigs]={"JetET50","JetET80"}; xmin=750; xmax=1400;
  //TString trignames[ntrigs]={"JetET80","JetET110"}; xmin=1000; xmax=2000;
  //TString trignames[ntrigs]={"JetET110","JetET150"}; xmin=1200; xmax=2000;
  TH1 *hists[ntrigs];

  int n=ntrigs;
  //n=5;
  //if (do1pb) n=4;
  TString whichAlgo = "_Scone7";
  TString whichJets;
  whichJets="_Gen_";
  // whichJets="_CaloMCCor_";

  TString xtitle="Dijet Mass (GeV)";
  TString ytitle="dN/dm";

  Double_t xl1=.54, yl1=0.59, xl2=xl1+.22, yl2=yl1+.25;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);

  cout << ymin << " " << ymax << endl;
  h2 = new TH2F("h2"," ",100,xmin,xmax,100,ymin,ymax);
  h2->GetXaxis()->SetTitle(xtitle);
  h2->GetYaxis()->SetTitle(ytitle);
  h2->Draw();
  for (int itrig=0;itrig<n;itrig++){

    //TString rootname="../hsts/DijetHistograms_" + trignames[itrig] + whichJets + "Eff_eta3_wght_yb1.5_ptcut.root";
    TString rootname="../nhsts/DiJets_" + trignames[itrig] + whichAlgo + whichJets + "eta3_wght1";
    if (do1pb) rootname = rootname + "_1pb-1";
    rootname = rootname + ".root";
    TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;

    //TString hname="JetPtMax";
    TString hname="JetMass1";
    if (combineTrgs){
      hname="JetMassComp1";
    }
    //hists[itrig]=GetHist(rootfile,hname); if (!hists[itrig]) return;
    hists[itrig]=get1DHist(rootfile,"",hname); if (!hists[itrig]) return;

    TH1F *h = hists[itrig];
    TH1F *hcomb;

    divideByBinWidth(h);
    Int_t rebin=1;
    h->Rebin(rebin);
    TString chopt="e,same";
    if (itrig ==0){
      h->GetXaxis()->SetRangeUser(xmin,xmax);
      h->GetYaxis()->SetRangeUser(ymin,ymax);
      chopt="hist";
      TH1F *htemp = h->Clone(); 
      htemp->SetName(hname="_combined"); 
      hcomb=htemp;
      hcomb->Reset();
    }
    hcomb->Add(hcomb,h,1.,1.);
    if (! combineTrgs) {
      double xmint=mintrg[itrig]; double xmaxt=maxtrg[itrig];
      h->GetXaxis()->SetRangeUser(xmint,xmaxt);
      h->SetMarkerStyle(trigmrks[itrig]);
      h->SetLineColor(trigcols[itrig]);
      h->SetMarkerColor(trigcols[itrig]);
      leg->AddEntry(h,trignames[itrig],"pl");
      h->Draw("e,same");
    }
  }

  if (combineTrgs) {
    //hcomb->GetYaxis()->SetRangeUser(ymin,ymax);
    xmin=300.;

    double minf=300.0, maxf=2500.;
    int npar=4;
    TF1 *xsfit = new TF1("xsfit",fitfunc,minf,maxf,npar); 
    xsfit->SetParameter(0,2e18);
    xsfit->SetParameter(1,-6);
    xsfit->SetParameter(2,-0.005);
    xsfit->SetParameter(3,14000.);
    xsfit->FixParameter(3,14000.);
    xsfit->SetLineWidth(2);

    hcomb->SetMarkerStyle(combmrk);
    hcomb->SetLineColor(combcol);
    hcomb->SetMarkerColor(combcol);
    hcomb->GetXaxis()->SetRangeUser(xmin,xmax);
    //hcomb->Draw("same");

    cout << "Now fitting" << endl;
    hcomb->Fit("xsfit","R");
    cout << "Done fitting" << endl;

  }else{
    leg->SetTextSize(0.035);
    leg->Draw();
  }

  if (compD0){

    const int nd0=15;
    Double_t xd0[nd0] = {209.1, 229.2, 253.3, 283.4, 309.3, 333.6, 367.6, 407.8, 447.9, 488.0, 528., 572., 638.9, 739.2, 873.2};
    Double_t yd0[nd0] = {3.78e-2, 2.10e-2, 1.16e-2, 6.18e-3, 3.55e-3, 2.12e-3, 1.18e-3, 5.84e-4, 2.89e-4, 1.64e-4, 
			 8.74e-5, 4.49e-5, 1.73e-5, 4.58e-6, 2.39e-7};
    
    double normD0=2.4e6;
    for (Int_t i=0; i<nd0; i++) { // normalize to CMS spectrum
      yd0[i]=yd0[i]*normD0;
    }

    TGraph *gr1 = new TGraph (nd0, xd0, yd0);
    gr1->SetMarkerStyle(21);

    double minf=350.0, maxf=800.;
    int npar=4;
    TF1 *xsfitd0 = new TF1("xsfitd0",fitfunc,minf,maxf,npar); 
    xsfitd0->SetParameter(0,3e18);
    xsfitd0->SetParameter(1,-6);
    xsfitd0->SetParameter(2,-0.05);
    xsfitd0->SetParameter(3,14000.);
    xsfitd0->FixParameter(3,14000.);

    xsfitd0->SetLineWidth(2);
    cout << " now d0 " << endl;
    gr1->Fit("xsfitd0","R");

    gr1->Draw("same,p");

  }
  if (saveIt){
    TString psfile="DiJetMass";
    if (combineTrgs) 
      psfile = psfile + "_compositeTrg";
    else
      psfile = psfile + "_allTrg";

    if (do1pb) psfile = psfile + "_1pb-1";
    psfile=psfile + ".eps";
    cout << "Writing file: " << psfile << endl;
    c1->Print(psfile);
  }
}
