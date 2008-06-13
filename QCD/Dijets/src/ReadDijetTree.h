#ifndef READDIJET_H
#define READDIJET_H

#include <map>
#include <string>
#include <vector>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include <TLorentzVector.h>

using namespace std;

//bool weightXS=true;
//bool weightXS=false;

const int NmassBins = 57;
const Float_t massBins[NmassBins+1] = {10,20,30,40,50,60,70,80,90,100,115,130,145,160,180,200,225,250,280,
				  310,340,370,410,450,500,550,610,670,730,790,870,960,1060,1170,1290,
				  1420,1570,1730,1900,2100,2300,2530,2790,3070,3380,3720,4100, 4510,
				  4960,5460,6010,6610,7280,8010,8810,10000,11000,12000};

const int NmassBins1 = 50;
const Float_t massBins1[NmassBins1+1] = {0,
					 100,120,140,160,180,200,220,240,260,280,
					 300,320,340,360,380,400,420,440,460,480,
					 500,520,540,560,580,600,620,640,660,680,
					 700,750,800,850,900,950,1000,1050,1100,1150,
					 1200,1300,1400,1500,1600,1800,2000,2250,2500,3000};


const int NchiBins1 = 20;
const Float_t chiBins1[NchiBins1+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
const int NchiBins2 = 15;
const Float_t chiBins2[NchiBins2+1] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20};
const int NchiBins3 = 11;
const Float_t chiBins3[NchiBins3+1] = {0,1,2,4,6,8,10,12,14,16,18,20};
const int NchiBins4 = 7;
const Float_t chiBins4[NchiBins4+1] = {0,1,3,6,9,12,15,20};


TH1F *hResB1,*hResB2,*hResB3,*hResB4;
TH1F *hResE1,*hResE2,*hResE3,*hResE4;
TH1F *hResF1,*hResF2,*hResF3,*hResF4;

void loadVar();
double jetResolution(const double, const double);
void fillJetResolutionPlots(const double, const double, const double, const double);
double scaleJetEnergies(const double);
vector <TLorentzVector> sortJets(const vector <TLorentzVector>);
// dijet cross section functions
const int npar=3;
double parsCMS[npar]={2.6e+18,-5.8,-.0581};
double parsD0[npar]={5.42e+18,-5.8,20.};
double xs(const double xx, const double *pars);

string jetAlgo;
string triggerStream;
string whichJets;
string outputDir, hComment;

int nevt,engShift;
bool weightXS, smrGenEng, adjustSpectrum;
double sclFactor, do1pb;

TRandom Trand(129344);

// tree variables
const int njmax=10;
float px[njmax],py[njmax],pz[njmax],e[njmax];
float dphi,m;
int njets;

// histogram pointers
TH1F *hNjets;

// jet and dijet variables
//TLorentzVector jetP4[njmax], dijetP4;
TLorentzVector dijetP4;
//vector<TLorentzVector> vjetP4; 

double mass_,rapidity_,etaStar_,etaBoost_,cosThetaStar_,chi_,dphi_;
double yStar_,yBoost_;
double  pt1, pt2, rap1, rap2, eta1, eta2, phi1, phi2;

const double CHI_CUT = 20.; 
//const double ETA_CUT = 10.0; 
//const double ETA_CUT = 3.0; 
//const double ETA_CUT = 2.0; 

double ETA_CUT;
double ETASTAR_CUT;
double ETABOOST_CUT;

// const double COSTHETASTAR_CUT = 10.; 
const double MASS_MIN = 40;
const double MASS_MAX = 12000;

//
double TRIGGER_ET_CUT;
double TRIGGER_PRESCALE;
double TRIGGER_XS_WEIGHT;
double TRIGGER_NEVT;
double TRIGGER_MINMASS;
double TRIGGER_MAXMASS;

/*
const double Jet20ptCut = 1;
const double Jet30ptCut = 1;
const double Jet50ptCut = 1;
const double Jet80ptCut = 1;
const double Jet110ptCut = 1;
const double Jet150ptCut = 1;
*/

const double Jet20ptCut  =  80;
const double Jet30ptCut  = 125;
const double Jet50ptCut  = 190;
const double Jet80ptCut  = 250;
const double Jet110ptCut = 310;
const double Jet150ptCut = 400;

const double Jet20prescl_all  =  6370;
const double Jet30prescl_all  =  1213;
const double Jet50prescl_all  =   152;
const double Jet80prescl_all  =  26.6;
const double Jet110prescl_all =  4.99;
const double Jet150prescl_all =  1.00;

const double Jet20prescl_1pb  =  29.65;
const double Jet30prescl_1pb  =   5.16;
const double Jet50prescl_1pb  =   1.00;
const double Jet80prescl_1pb  =   1.00;
const double Jet110prescl_1pb =   1.00;
const double Jet150prescl_1pb =   1.00;

const double Jet20xswght=101600000.;
const double Jet30xswght= 21550000.;
const double Jet50xswght=  2484000.;
const double Jet80xswght=   323700.;
const double Jet110xswght=   88730.;
const double Jet150xswght=   17120.;

const double Jet20Nevt_all= 3926600;
const double Jet30Nevt_all= 4131600;
const double Jet50Nevt_all= 4010400;
const double Jet80Nevt_all= 2891200;
const double Jet110Nevt_all=3980000;
const double Jet150Nevt_all=3885600;

const double Jet20Nevt_1pb= 3926600;
const double Jet30Nevt_1pb= 4131600;
const double Jet50Nevt_1pb= 2506500;
const double Jet80Nevt_1pb=  324850;
const double Jet110Nevt_1pb=  88839;
const double Jet150Nevt_1pb=  17121;

const double Jet20AllSamples= 3131233;
const double Jet30AllSamples= 4342174;
const double Jet50AllSamples= 4022000;
const double Jet80AllSamples= 2554579;
const double Jet110AllSamples=3848194;
const double Jet150AllSamples=5213220;

const double Jet20MassRng[2]={380,600};
const double Jet30MassRng[2]={600,900};
const double Jet50MassRng[2]={900,1160};
const double Jet80MassRng[2]={1160,1460};
const double Jet110MassRng[2]={1460,1880};
const double Jet150MassRng[2]={1880,12000};
//const double Jet150MassRng[2]={3000,12000};


vector<std::string> hltnames;
std::map<string,double> hltTrigEtCut, hltPrescale,hltXSweight,hltAllSmpls,hltNevents, 
  hltTrigMinMass,hltTrigMaxMass;
std::map<string,double>::iterator hltTrig_iter;
typedef std::map<string,double>::value_type valType;

void fillTriggerInfo();
void fillTriggerInfo(){

  hltnames.push_back("Jet20");
  hltnames.push_back("Jet30");
  hltnames.push_back("Jet50");
  hltnames.push_back("Jet80");
  hltnames.push_back("Jet110");
  hltnames.push_back("Jet150");

  string trigName;
  double Jet20prescl=Jet20prescl_all, Jet30prescl=Jet30prescl_all, Jet50prescl=Jet50prescl_all, 
    Jet80prescl=Jet80prescl_all, Jet110prescl= Jet110prescl_all, Jet150prescl= Jet150prescl_all;
  if (do1pb){
    Jet20prescl = Jet20prescl_1pb;
    Jet30prescl = Jet30prescl_1pb;
    Jet50prescl = Jet50prescl_1pb;
    Jet80prescl = Jet80prescl_1pb; 
    Jet110prescl= Jet110prescl_1pb;
    Jet150prescl= Jet150prescl_1pb;
  }

  double Jet20Nevt=Jet20Nevt_all, Jet30Nevt=Jet30Nevt_all, Jet50Nevt=Jet50Nevt_all, 
    Jet80Nevt=Jet80Nevt_all, Jet110Nevt= Jet110Nevt_all, Jet150Nevt= Jet150Nevt_all;
  if (do1pb){
    Jet20Nevt = Jet20Nevt_1pb;
    Jet30Nevt = Jet30Nevt_1pb;
    Jet50Nevt = Jet50Nevt_1pb;
    Jet80Nevt = Jet80Nevt_1pb; 
    Jet110Nevt =Jet110Nevt_1pb; 
    Jet150Nevt= Jet150Nevt_1pb;
  }

  trigName=hltnames.at(0);
  hltTrigEtCut.insert(valType(trigName,Jet20ptCut));  
  hltTrigMinMass.insert(valType(trigName,Jet20MassRng[0]));  
  hltTrigMaxMass.insert(valType(trigName,Jet20MassRng[1]));  
  hltPrescale.insert(valType(trigName,Jet20prescl));  
  hltXSweight.insert(valType(trigName,Jet20xswght));
  hltNevents.insert(valType(trigName,Jet20Nevt));
  hltAllSmpls.insert(valType(trigName,Jet20AllSamples));

  trigName=hltnames.at(1);
  hltTrigEtCut.insert(valType(trigName,Jet30ptCut));  
  hltTrigMinMass.insert(valType(trigName,Jet30MassRng[0]));  
  hltTrigMaxMass.insert(valType(trigName,Jet30MassRng[1]));  
  hltPrescale.insert(valType(trigName,Jet30prescl));
  hltXSweight.insert(valType(trigName,Jet30xswght));  
  hltNevents.insert(valType(trigName,Jet30Nevt));
  hltAllSmpls.insert(valType(trigName,Jet30AllSamples));

  trigName=hltnames.at(2);
  hltTrigEtCut.insert(valType(trigName,Jet50ptCut));  
  hltTrigMinMass.insert(valType(trigName,Jet50MassRng[0]));  
  hltTrigMaxMass.insert(valType(trigName,Jet50MassRng[1]));  
  hltPrescale.insert(valType(trigName,Jet50prescl));  
  hltXSweight.insert(valType(trigName,Jet50xswght));
  hltNevents.insert(valType(trigName,Jet50Nevt));
  hltAllSmpls.insert(valType(trigName,Jet50AllSamples));

  trigName=hltnames.at(3);
  hltTrigEtCut.insert(valType(trigName,Jet80ptCut));  
  hltTrigMinMass.insert(valType(trigName,Jet80MassRng[0]));  
  hltTrigMaxMass.insert(valType(trigName,Jet80MassRng[1]));  
  hltPrescale.insert(valType(trigName,Jet80prescl));  
  hltXSweight.insert(valType(trigName,Jet80xswght));
  hltNevents.insert(valType(trigName,Jet80Nevt));
  hltAllSmpls.insert(valType(trigName,Jet80AllSamples));
  trigName=hltnames.at(4);

  hltTrigEtCut.insert(valType(trigName,Jet110ptCut)); 
  hltTrigMinMass.insert(valType(trigName,Jet110MassRng[0]));  
  hltTrigMaxMass.insert(valType(trigName,Jet110MassRng[1]));  
  hltPrescale.insert(valType(trigName,Jet110prescl));   
  hltXSweight.insert(valType(trigName,Jet110xswght));
  hltNevents.insert(valType(trigName,Jet110Nevt));
  hltAllSmpls.insert(valType(trigName,Jet110AllSamples));

  trigName=hltnames.at(5);
  hltTrigEtCut.insert(valType(trigName,Jet150ptCut));
  hltTrigMinMass.insert(valType(trigName,Jet150MassRng[0]));  
  hltTrigMaxMass.insert(valType(trigName,Jet150MassRng[1]));  
  hltPrescale.insert(valType(trigName,Jet150prescl));
  hltXSweight.insert(valType(trigName,Jet150xswght));
  hltNevents.insert(valType(trigName,Jet150Nevt));
  hltAllSmpls.insert(valType(trigName,Jet150AllSamples));
}

double getTriggerEtCut(const string&);
double getTriggerEtCut(const string& tStream){

  double etCut=1e10;

  std::map<string,double>::iterator iter=hltTrigEtCut.find(tStream);
  if (iter==hltTrigEtCut.end())
    std::cout << "%getTriggerETCut -- Could not find stream with name: " << tStream << std::endl;
  else
    etCut=iter->second;

  return etCut;
}

double getTriggerPrescl(const string&);
double getTriggerPrescl(const string& tStream){

  double preScale=0.;

  std::map<string,double>::iterator iter=hltPrescale.find(tStream);
  if (iter==hltPrescale.end())
    std::cout << "%getTriggerPrescl -- Could not find stream with name: " << tStream << std::endl;
  else
    preScale=iter->second;

  return preScale;
}

double getTriggerXSWeight(const string&);
double getTriggerXSWeight(const string& tStream){

  double xsWeight=0.;

  std::map<string,double>::iterator iter=hltXSweight.find(tStream);
  if (iter==hltXSweight.end())
    std::cout << "%getTriggerXSweight -- Could not find stream with name: " << tStream << std::endl;
  else
    xsWeight=iter->second;

  return xsWeight;
}

double getTriggerAllSamples(const string&);
double getTriggerAllSamples(const string& tStream){

  double N=0.;

  std::map<string,double>::iterator iter=hltAllSmpls.find(tStream);
  if (iter==hltAllSmpls.end())
    std::cout << "%getTriggerAllSmples -- Could not find stream with name: " << tStream << std::endl;
  else
    N=iter->second;

  return N;
}

double getTriggerNevt(const string&);
double getTriggerNevt(const string& tStream){

  double N=0.;

  std::map<string,double>::iterator iter=hltNevents.find(tStream);
  if (iter==hltNevents.end())
    std::cout << "%getTriggerAllSmples -- Could not find stream with name: " << tStream << std::endl;
  else
    N=iter->second;

  return N;
}

double getTriggerMinMass(const string&);
double getTriggerMinMass(const string& tStream){

  double mass=-999.;

  std::map<string,double>::iterator iter=hltTrigMinMass.find(tStream);
  if (iter==hltTrigMinMass.end())
    std::cout << "%getTriggerMinMass -- Could not find stream with name: " << tStream << std::endl;
  else
    mass=iter->second;

  return mass;
}

double getTriggerMaxMass(const string&);
double getTriggerMaxMass(const string& tStream){

  double mass=-999.;

  std::map<string,double>::iterator iter=hltTrigMaxMass.find(tStream);
  if (iter==hltTrigMaxMass.end())
    std::cout << "%getTriggerMaxMass -- Could not find stream with name: " << tStream << std::endl;
  else
    mass=iter->second;

  return mass;
}

TH1F* bookTH1F(const TString& hname, const TString& htitle, const int nbins, const Double_t xmin, const Double_t xmax);
TH1F* bookTH1F(const TString& hname, const TString& htitle, const int nbins, const Double_t xmin, const Double_t xmax){
  TH1F *h = new TH1F (hname,htitle,nbins,xmin,xmax);
  h->Sumw2();
  return h;
}

TH1F* bookTH1F(const TString& hname, const TString& htitle, const int nbins, const Float_t xbins[]);
TH1F* bookTH1F(const TString& hname, const TString& htitle, const int nbins, const Float_t xbins[]){
  TH1F *h = new TH1F (hname,htitle,nbins,xbins);
  h->Sumw2();
  return h;
}

#endif

void readInputCards(const string&);
void readInputCards(const string& inpc){

  cout << "Input cards read from: " << inpc << endl;
  ifstream inputCards(inpc.c_str());

  char temp[100];
  inputCards.getline(temp, 100);  jetAlgo=temp;
  inputCards.getline(temp, 100);  outputDir=temp;
  inputCards.getline(temp, 100);  triggerStream=temp;
  inputCards.getline(temp, 100);  nevt=strtol(temp,0,10);
  inputCards.getline(temp, 100);  whichJets=temp;
  inputCards.getline(temp, 100);  ETA_CUT=strtof(temp,0);
  inputCards.getline(temp, 100);  weightXS=strtol(temp,0,10);
  inputCards.getline(temp, 100);  engShift=strtol(temp,0,10);
  inputCards.getline(temp, 100);  smrGenEng=strtol(temp,0,10);
  inputCards.getline(temp, 100);  do1pb=strtol(temp,0,10);
  inputCards.getline(temp, 100);  adjustSpectrum=strtol(temp,0,10);

  if (whichJets != "Gen") smrGenEng=false;

  /*
  cout << "Jet Algorithm: " << jetAlgo << "\n";
  cout << "Trigger Stream: " << triggerStream << "\n";
  cout << "nevt: " << nevt << "\n";
  cout << "whichJets: " << whichJets << "\n";
  cout << "ETA_CUT: " << ETA_CUT << "\n";
  if (weightXS)
    cout << "Weighting Cross Sections" << endl;
  else
    cout << "Not weighting Cross Sections" << endl;
  cout << "" << endl;
  */

  //float fl;
  //fl=strtof(temp);
  //int i;
  //someVarName >> i;


}

string createOutputFileName();
string createOutputFileName(){

  ostringstream ch_eta(""), ch_xsWeight(""), ch_engShift("");
  ch_eta << ETA_CUT;
  ch_xsWeight << weightXS;
  ch_engShift << engShift;

  hComment="xxx";
  const char *hcomm = gSystem->Getenv("HCOMMENT");
  if (hcomm) hComment = hcomm;

  string basename="DiJets_";
  string outfile=outputDir + "/" + basename;
  outfile= outfile + triggerStream + "_";
  outfile= outfile + jetAlgo + "_";
  outfile= outfile + whichJets + "_";

  if (abs(engShift) > 0){
    outfile= outfile + "engShift" + ch_engShift.str() + "_";
  }

  if (smrGenEng)
    outfile= outfile + "smrEng_";

  if (adjustSpectrum)
    outfile= outfile + "adjMass_";

  outfile= outfile + "eta" + ch_eta.str() + "_";
  outfile= outfile + "wght" + ch_xsWeight.str();

  if (do1pb) outfile= outfile + "_1pb-1";

  if (hComment != "xxx")
    outfile= outfile + "_" + hComment;
  outfile= outfile + "_v2.root";
  return outfile;

}

double Phi_0_2pi(double x) {
  while (x >= 2*M_PI) x -= 2*M_PI;
  while (x <     0.)  x += 2*M_PI;
  return x;
}

class PtGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.Pt() > j.Pt());
  }
};
