#include "rootfuncs.h"

double invMass(double *xx, double *pars);

double invMass(double *xx, double *pars)
{

  double x=xx[0];
  double et=pars[0];

  double value=2.*pow(et,2)*(cosh(log(x))+1);
  double fitval=value;

  return sqrt(fitval);
}

void plotMvsChi(){

  gROOT->Reset();
  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetOptStat(1);

  TString thresh="30";

  TString rootname="../hsts/DijetHistograms_JetET" + thresh + "_Gen_Eff_eta2_tst.root";
  //TString rootname="../hsts/DijetHistograms_JetET" + thresh + "_Gen_Eff.root";
  TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;

  TString hname="MassVsChi";
  TH1 *h1=GetHist(rootfile,hname); if (!h1) return;

  TF1 *f1 = new TF1("f1",invMass,1.,30.,1);

  double etcut=0.;
  if (thresh == "110"){
    etcut=310;
  }else if (thresh == "80"){
    etcut=250;
  }else if (thresh == "50"){
    etcut=200;
  }else if (thresh == "30"){
    etcut=125;
  }else if (thresh == "20"){
    etcut=60;
  }
  f1->SetParameter(0,etcut);

  f1->GetXaxis()->SetTitle("#chi");
  f1->GetYaxis()->SetTitle("Invariant Mass (GeV)");
  f1->GetYaxis()->SetRangeUser(60,1400);

  f1->SetLineColor(kRed);
  f1->SetLineStyle(2);

  h1->Draw();
  f1->Draw("same");

}
