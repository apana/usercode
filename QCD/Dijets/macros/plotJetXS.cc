#include "rootfuncs.h"


double xs(double *xx, double *pars)
{

  double value=1;
  float N=1.e13;
  //float a=-4;
  float a=-5;
  float b=1.;

  const double roots=14000.;

  double pt=xx[0];
  double value=N*pow(pt,a)*pow((1-2*pt/roots),b);
  //double value=N*pow(pt,a);
  //double value=pow((1-2*pt/roots),b);
  return value;
}

void plotJetXS(){

  mySetup(1); // 0 -- linear plots, 1 -- log plots  

  double ptmin=50, ptmax=3500;
  double xsmin=1e-6; xsmax=1e4;


  TF1 *f1 = new TF1("f1",xs,ptmin,ptmax,1.);
  TH2F *h2 = new TH2F("h2"," ",100,ptmin,ptmax,100,xsmin,xsmax);

  h2->Draw();
  f1->GetYaxis()->SetRangeUser(1e-6,1e4);
  f1->Draw("same");

}
