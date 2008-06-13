#include "rootfuncs.h"


double resol(double *xx, double *pars)
{

  double a,b,c;
  // data-driven barrel
  //a=4.7;
  //b=1.3;
  //c=0.044;

  // mc truth barrel
  //a=4.;
  // b=1.4;
  //c=0.046;

  // mc truth endcap
  a=4.5;
  b=1.2;
  c=0.057;

  // mc truth forward
  //a=4.;
  //b=0;
  //c=0.11;

  double pt=xx[0];

  double aa=pow(a/pt,2);
  double bb=pow(b/sqrt(pt),2);
  double cc=pow(c,2);

  double value=sqrt(aa+bb+cc);

  return value;
}

void plotResol(){

  mySetup(0); // 0 -- linear plots, 1 -- log plots  

  h2 = new TH2F("h2"," ",100,0,799,100,0,.45);
  TF1 *f1 = new TF1("f1",resol,20,799,1.);
  f1->GetXaxis()->SetTitle("p_{T}");
  f1->GetYaxis()->SetTitle("#sigma/p_{T}");

  h2->Draw();
  f1->Draw("same");

}
