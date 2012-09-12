#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <cmath>

#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLorentzVector.h"

using namespace std;


// forward declaration  
double deltaPhi(double phi1, double phi2); 
double evalEt(double pt, double eta, double phi, double e);
double evalMt(double pt, double eta, double phi, double e);
double ptHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1, double e1, double oldpt0, double oldpt1);
double massHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1, double e1, double oldpt0, double oldpt1);
double resolutionBias(double eta);
double METdeltaPhi(double metphi, double jet1phi, double jet2phi);
double deltaR2(double eta1, double phi1, double eta2, double phi2);
double deltaR(double eta1, double phi1, double eta2, double phi2); 
bool deltaRCone(double eta1, double phi1, double eta2, double phi2); 



double deltaPhi(double phi1, double phi2) { 
  double PI = 3.14159265;
  double result = phi1 - phi2;
  while (result > PI) result -= 2*PI;
  while (result <= -PI) result += 2*PI;
  return result;
}

double evalEt( double pt, double eta, double phi, double e){
  TLorentzVector j;
  j.SetPtEtaPhiE(pt,eta,phi,e);
  return j.Et(); 
}

double evalMt( double pt, double eta, double phi, double e){
  TLorentzVector j;
  j.SetPtEtaPhiE(pt,eta,phi,e);
  return j.Mt(); 
}
 
double ptHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1 , double e1, double oldpt0, double oldpt1) {
  TLorentzVector j0, j1 , H;

  j0.SetPtEtaPhiE(pt0,eta0,phi0,e0*pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1,e1*pt1/oldpt1);
  H = j0 +j1;
  return H.Pt(); 
}

double massHiggs(double pt0, double eta0, double phi0, double e0, double pt1,double  eta1, double phi1, double e1, double oldpt0, double oldpt1) {
  TLorentzVector j0, j1, H;
  j0.SetPtEtaPhiE(pt0,eta0,phi0, e0 * pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1, e1  * pt1/oldpt1);
  H = j0 +j1;
  return fabs(H.M());
}

double resolutionBias(double eta) {
  if (eta< 1.1) return 0.05;
  if (eta< 2.5) return 0.10;
  if (eta< 5) return 0.30;
  return 0;
}

double METdeltaPhi(double metphi, double jet1phi, double jet2phi) {
  double METJ1dPhi;
  double METJ2dPhi;
  double METJdPhi;
  METJ1dPhi=fabs(deltaPhi(metphi, jet1phi));
  METJ2dPhi=fabs(deltaPhi(metphi, jet2phi));
  if (METJ1dPhi<METJ2dPhi) METJdPhi=METJ1dPhi;
  else METJdPhi=METJ2dPhi;
  return METJdPhi;
}


double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return deta*deta + dphi*dphi;
}


double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return std::sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}


bool deltaRCone(double eta1, double phi1, double eta2, double phi2) {
  
  cout << "Delta R = "<< std::sqrt(deltaR2 (eta1, phi1, eta2, phi2)) << endl;
  return std::sqrt(deltaR2 (eta1, phi1, eta2, phi2)) < 0.4;
  
}



double MAX(double val1, double val2) {
  return TMath::Max(val1,val2);
}

double MIN(double val1, double val2) {
  return TMath::Min(val1,val2);
}
