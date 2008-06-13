#include "rootfuncs.h"


double chi(double *xx, double *pars)
{

  double etastar=xx[0];

  double value=exp(2*etastar);

  return value;
}

void plotChivsEtaStar(){

  mySetup(0); // 0 -- linear plots, 1 -- log plots  

  TF1 *f1 = new TF1("f1",chi,0.3,1.5,1.);
  f1->GetXaxis()->SetTitle("#eta^{*}");
  f1->GetYaxis()->SetTitle("#chi");
  f1->Draw();

}
