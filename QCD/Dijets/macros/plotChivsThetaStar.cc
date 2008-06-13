#include "rootfuncs.h"


double chi(double *xx, double *pars)
{

  double theta=xx[0];
  double costheta=cos(theta);
  double value=(1+costheta)/(1-costheta);

  return value;
}

void plotChivsThetaStar(){

  mySetup(0); // 0 -- linear plots, 1 -- log plots  

  double thetamin=0.3, thetamax=1.57;
  thetamax=.6;
  TF1 *f1 = new TF1("f1",chi,thetamin,thetamax,1.);
  
  f1->Draw();

}
