#include "rootfuncs.h"

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

  const double roots=2000.;


  //double value=N*pow(pt,a)*pow((1-2*pt/roots),b);
  double value=N*pow(pt,a);
  //double value=pow((1-2*pt/roots),b);
  return value;
}

// d0 dijet cross section PRL v82, n12 p2460
void d0xs(){

  gROOT->Reset();
  mySetup(1); // 0 -- linear plots, 1 -- log plots  
  const int nd0=15;
  Double_t xd0[nd0] = {209.1, 229.2, 253.3, 283.4, 309.3, 333.6, 367.6, 407.8, 447.9, 488.0, 528., 572., 638.9, 739.2, 873.2};
  Double_t yd0[nd0] = {3.78e-2, 2.10e-2, 1.16e-2, 6.18e-3, 3.55e-3, 2.12e-3, 1.18e-3, 5.84e-4, 2.89e-4, 1.64e-4, 
		       8.74e-5, 4.49e-5, 1.73e-5, 4.58e-6, 2.39e-7};

  Double_t exd0[nd0], eyd0[nd0];
  for (Int_t i=0; i<nd0; i++) {
    cout << xd0[i] << " " << yd0[i]<< endl;
    // fake the errors for now
    exd0[i]=0.;
    eyd0[i]=0.2*yd0[i];
  }

  double minf=200.0, maxf=400.;
  int npar=2;
  TF1 *xsfit = new TF1("xsfit",fitfunc,200.,300.,2); 
  xsfit->SetParameter(0,1e13);
  xsfit->SetParameter(1,-4);

  TGraphErrors *gr1 = new TGraphErrors (nd0, xd0, yd0,exd0,eyd0);
  gr1->SetMarkerStyle(21);
  gr1->Draw("AP");
  
  //histo->Fit("f1","R");
  gr1->Fit("xsfit");
}
