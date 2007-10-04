#ifndef __CINT__
#include <iostream>
#include <sstream>
#include <cmath>
#include "Riostream.h"
#include <vector>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TChain.h"
#endif
#include <map>

//  constants
const double convPbtoCm2 = 1.E-36;
const int nhltmax=500;

const float mineta=-5.4;
const float maxeta=5.4;
//const float mineta=-2.;
//const float maxeta=2.;


int nhlt;
vector<string> hltnames;
std::map<string,int> hltTriggers;
std::map<string,int>::iterator hltTrig_iter;
typedef std::map<string,int>::value_type valType;


// ### tree variables ####
float metpt,metphi,metsum;
float genmetpt;

float htpt,htsum;
int nl1cenjet,nl1forjet,nl1taujet;

float mcpthat;

int njetcal,njetgen;
const int kMaxJetCal = 10000;
float caljetpt[kMaxJetCal],caljeteta[kMaxJetCal];
float genjetpt[kMaxJetCal],genjeteta[kMaxJetCal];
  
const int kMaxL1Jet = 4;
float l1cenjetpt[kMaxL1Jet],l1forjetpt[kMaxL1Jet],l1taujetpt[kMaxL1Jet];
float l1met,l1metphi,l1mettot,l1mettothad;

// Level-1 bits
int i_l1bit15,i_l1bit30,i_l1bit50,i_l1bit70,i_l1bit100,i_l1bit150,i_l1bit200;
int i_etm20,i_etm30,i_etm40,i_etm50,i_etm60;
int i_l1htt100,i_l1htt200,i_l1htt250,i_l1htt300,i_l1htt400,i_l1htt500;

// HLT Trigger names
string hlt1jet,hlt1jetPre1,hlt1jetPre2,hlt1jetPre3,hlt1jetPre4;
string hlt1met,hlt1metPre1,hlt1metPre2,hlt1metPre3;
string hlt2jet,hlt3jet,hlt4jet;

int n_dijet200;
int n_1jet;
int n_dijet200_1jet;

//  ### end of tree variables ###

// Histogram pointers //
TString hname,htitle;
TObjArray *Hlist = new TObjArray();

TH1F* h_l1;

TH1F* h_gen[nhltmax];
TH1F* h_cal[nhltmax];

TH1F *jetmet_jet0, *jetmet_jet50, *jetmet_jet80, *jetmet_jet140, *jetmet_jet180;
TH1F *jetmet_l1_jet0, *jetmet_l1_jet50, *jetmet_l1_jet80, *jetmet_l1_jet140, *jetmet_l1_jet180;

TH1F *met_Untriggered, *met_trg, *met_trgPre1, *met_trgPre2, *met_trgPre3, *met_etm20, *met_etm30, *met_etm40, *met_etm50, *met_etm60;
TH1F *met_rat;
TH1F *summet_Untriggered;
TH1F *sumht_Untriggered;
TH1F *sumhtjets_Untriggered;

TH1F *h_genmet, *met_Untriggered_gencut;

TH1F *h_l1met,*h_l1sumet,*h_l1sumht, *h_l1sumht_l1htt200;
TH1F *h_l1_metDphi_1,*h_l1_metDphi_2,*h_l1_metDphi_3,*h_l1_metDphi_4,*h_l1_metDphi_5;
TH1F *h_l1_metRat_1,*h_l1_metRat_2,*h_l1_metRat_3,*h_l1_metRat_4,*h_l1_metRat_5;
TH1F *h_l1met_1,*h_l1met_2,*h_l1met_3;

TH2F *htvspt, *htvsavept, *htvsl1ht;
TH2F *h_l1vsGen, *h_l1metvsMet;

TH1F *h_caljetpt;
TH1F *h_genjetpt;

TH1F *mx_caljetpt_Untriggered;  
TH1F *mx_caljetpt_l15, *mx_caljetpt_l30, *mx_caljetpt_l50, *mx_caljetpt_l70, *mx_caljetpt_l100, *mx_caljetpt_l150, *mx_caljetpt_l200;
TH1F *mx_caljetpt_m15, *mx_caljetpt_m30, *mx_caljetpt_m50, *mx_caljetpt_m70, *mx_caljetpt_m100, *mx_caljetpt_m150, *mx_caljetpt_m200;

TH1F *mx_genjetpt_Untriggered;  
TH1F *mx_genjetpt_l15, *mx_genjetpt_l30, *mx_genjetpt_l50, *mx_genjetpt_l70, *mx_genjetpt_l100, *mx_genjetpt_l150, *mx_genjetpt_l200;
TH1F *mx_genjetpt_m15, *mx_genjetpt_m30, *mx_genjetpt_m50, *mx_genjetpt_m70, *mx_genjetpt_m100, *mx_genjetpt_m150, *mx_genjetpt_m200;

TH1F *mx_genjetpt_hlt30,*mx_genjetpt_hlt60,*mx_genjetpt_hlt110,*mx_genjetpt_hlt150,*mx_genjetpt_hlt200,*mx_genjetpt_hlt250;
TH1F *mx_genjetpt_hlt65,*mx_genjetpt_hlt70;

TH1F *nl1jetc, *nl1jetf, *nl1jett;
TH1F *mx_l1cenpt, *mx_l1forpt, *mx_l1taupt;

TH1F *hlt_rates, *met_rates;

TH1F *h_dijetpt, *h_dijetpt_dijet30, *h_dijetpt_dijet60, *h_dijetpt_dijet110, *h_dijetpt_dijet150, *h_dijetpt_dijet200;
TH1F *h_dijetpt_l15, *h_dijetpt_l30, *h_dijetpt_l70, *h_dijetpt_l100, *h_dijetpt_l150;
///////////
// Files 
///////////

const Int_t Nfil = 21; // Total number of files, of different x-section, to be read
vector<TString>ProcFil;

TChain* TabChain[Nfil];

vector<int>mcSample_bin;
double xsec[Nfil], skimEff[Nfil];


// ########## Functions ####################

void loadHLTNames(){

  hlt1jet     ="HLT1jet";  
  hlt1jetPre1 ="HLT1jetPE1"; hlt1jetPre2 ="HLT1jetPE3";  hlt1jetPre3 ="HLT1jetPE5";  hlt1jetPre4 ="CandHLT1jetPE7";
  hltnames.push_back(hlt1jet);
  hltnames.push_back(hlt1jetPre1);
  hltnames.push_back(hlt1jetPre2);
  hltnames.push_back(hlt1jetPre3);
  hltnames.push_back(hlt1jetPre4);

  hlt2jet="HLT2jet";  hlt3jet="HLT3jet";  hlt4jet="HLT4jet";
  hltnames.push_back(hlt2jet);
  hltnames.push_back(hlt3jet);
  hltnames.push_back(hlt4jet);

  hlt1met="HLT1MET";
  hlt1metPre1="CandHLT1METPre1";  hlt1metPre2="CandHLT1METPre2";  hlt1metPre3="CandHLT1METPre3";
  hltnames.push_back(hlt1met);
  hltnames.push_back(hlt1metPre1);
  hltnames.push_back(hlt1metPre2);
  hltnames.push_back(hlt1metPre3);

  hltnames.push_back("HLT1jet1MET");
  hltnames.push_back("HLT2jet1MET");
  hltnames.push_back("HLT3jet1MET");
  hltnames.push_back("HLT4jet1MET");

  hltnames.push_back("HLT2jetAco");
  hltnames.push_back("HLT1jet1METAco");

  hltnames.push_back("HLT1MET1HT");
  hltnames.push_back("HLT2jetvbfMET");

  hltnames.push_back("HLTS2jet1METNV");
  hltnames.push_back("HLTS2jet1METAco");
  //hltnames.push_back("HLT1HT");

  hltnames.push_back("CandHLT2jetAve30");
  hltnames.push_back("CandHLT2jetAve60");
  hltnames.push_back("CandHLT2jetAve110");
  hltnames.push_back("CandHLT2jetAve150");
  hltnames.push_back("CandHLT2jetAve200");

  nhlt=hltnames.size();

  for (int i=0; i<nhlt; i++ ){
    string trigName=hltnames.at(i);
    hltTriggers.insert(valType(trigName,0));
  }
}

TH1F* Book1dHist(const char* name, const char* title, Int_t nbins, Double_t xmin, Double_t xmax){
  TH1F *h= new TH1F(name,title,nbins,xmin,xmax);
  h->Sumw2();
  Hlist->Add(h);
  return h;
}

TH2F* Book2dHist(const char* name, const char* title, Int_t nbinsx, Double_t xmin, Double_t xmax,
		 Int_t nbinsy, Double_t ymin, Double_t ymax){
  TH2F *h= new TH2F(name,title,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  h->Sumw2();
  Hlist->Add(h);
  return h;
}

void bookHistograms(){


  const Int_t nmetbins=50;
  const Double_t minmet=0.;
  const Double_t maxmet=150.;

  h_l1=Book1dHist("h_l1","l1pt  45 < gen pt < 55 ",20,0.,2.5);
  h_l1vsGen=Book2dHist("h_l1vsGen","Leading l1pt  vs Leading genjet",50,0.,250,50,0.,250);
  h_l1metvsMet=Book2dHist("h_l1metvsMet","L1 Missing ET  vs Offline Missing ET",50,0.,250,50,0.,250);

  h_l1_metDphi_1=Book1dHist("h_l1_metDphi_1"," RecoMet #phi - L1Met #phi,  0 < MET < 10 ",50,-M_PI,M_PI);
  h_l1_metDphi_2=Book1dHist("h_l1_metDphi_2"," RecoMet #phi - L1Met #phi,  10 < MET < 20 ",50,-M_PI,M_PI);
  h_l1_metDphi_3=Book1dHist("h_l1_metDphi_3"," RecoMet #phi - L1Met #phi,  20 < MET < 30 ",50,-M_PI,M_PI);
  h_l1_metDphi_4=Book1dHist("h_l1_metDphi_4"," RecoMet #phi - L1Met #phi,  30 < MET < 50 ",50,-M_PI,M_PI);
  h_l1_metDphi_5=Book1dHist("h_l1_metDphi_5"," RecoMet #phi - L1Met #phi,  MET > 50 ",50,-M_PI,M_PI);

  h_l1_metRat_1=Book1dHist("h_l1_metRat_1"," L1MET/RecoMet --  0 < MET < 10 ",50,0.,5.);
  h_l1_metRat_2=Book1dHist("h_l1_metRat_2"," L1MET/RecoMet --  10 < MET < 20 ",50,0.,5.);
  h_l1_metRat_3=Book1dHist("h_l1_metRat_3"," L1MET/RecoMet --  20 < MET < 30 ",50,0.,5.);
  h_l1_metRat_4=Book1dHist("h_l1_metRat_4"," L1MET/RecoMet --  30 < MET < 50 ",50,0.,5.);
  h_l1_metRat_5=Book1dHist("h_l1_metRat_5"," L1MET/RecoMet --  MET > 50 ",50,0.,5.);

  met_Untriggered = Book1dHist("met_Untriggered","Met",nmetbins,minmet,maxmet);
  met_Untriggered_gencut = Book1dHist("met_Untriggered_gencut","Met -- gencut",nmetbins,minmet,maxmet);
  met_trg = Book1dHist("met_trg","Met -- HLT",nmetbins,minmet,maxmet);
  met_trgPre1 = Book1dHist("met_trgPre1","Met -- HLT",nmetbins,minmet,maxmet);
  met_trgPre2 = Book1dHist("met_trgPre2","Met -- HLT",nmetbins,minmet,maxmet);
  met_trgPre3 = Book1dHist("met_trgPre3","Met -- HLT",nmetbins,minmet,maxmet);
  met_etm20 = Book1dHist("met_etm20"," Met -- L1 etm>20 ",nmetbins,minmet,maxmet);
  met_etm30 = Book1dHist("met_etm30"," Met -- L1 etm>30 ",nmetbins,minmet,maxmet);
  met_etm40 = Book1dHist("met_etm40"," Met -- L1 etm>40 ",nmetbins,minmet,maxmet);
  met_etm50 = Book1dHist("met_etm50"," Met -- L1 etm>50 ",nmetbins,minmet,maxmet);
  met_etm60 = Book1dHist("met_etm60"," Met -- L1 etm>60 ",nmetbins,minmet,maxmet);

  met_rat = Book1dHist("met_rat"," Met pt / MetSum ",100,0.,1.);

  h_genmet = Book1dHist("h_genmet","Generated Met",nmetbins,minmet,maxmet);

  h_l1met = Book1dHist("h_l1met","L1 Met ",nmetbins,minmet,maxmet);
  h_l1met_1 = Book1dHist("h_l1met_1","L1 Met, Met #phi - L1MET #phi < 1.",nmetbins,minmet,maxmet);
  h_l1met_2 = Book1dHist("h_l1met_2","L1 Met, MaxGenJet in Barrel",nmetbins,minmet,maxmet);
  h_l1met_3 = Book1dHist("h_l1met_3","L1 Met, GenMET < 5.",nmetbins,minmet,maxmet);

  const Int_t nsumbins=50;
  const Double_t minsumet=0.;
  const Double_t maxsumet=1000.;
  
  summet_Untriggered = Book1dHist("summet_Untriggered","Sumet -- Untriggered",nsumbins,minsumet,maxsumet);
  sumht_Untriggered  = Book1dHist("sumht_Untriggered","Sumht -- Untriggered",nsumbins,minsumet,maxsumet);
  sumhtjets_Untriggered = Book1dHist("sumhtjets_Untriggered","Sumht from Jets -- Untriggered",nsumbins,minsumet,maxsumet);

  h_l1sumet = Book1dHist("l1sumEt","L1 SumEt ",nsumbins,minsumet,maxsumet);
  h_l1sumht = Book1dHist("l1sumHt","L1 SumHt ",nsumbins,minsumet,maxsumet);
  h_l1sumht_l1htt200 = Book1dHist("l1sumHt_l1htt200","L1 SumHt l1_htt200>0",nsumbins,minsumet,maxsumet);

  const Int_t nsumbinsx=25, nsumbinsy=50;
  const Double_t minsumetx=0.,minsumety=0.;
  const Double_t maxsumetx=500.,maxsumety=1000.;
  
  htvspt  = Book2dHist("sumht_vs_leadingpt","Sumht vs p_{T} of Leading Jet",nsumbinsx,minsumetx,maxsumetx,nsumbinsy,minsumety,maxsumety);

  htvsl1ht = Book2dHist("sumht_vs_l1sumht","Sumht vs L1 Sumht",nsumbins,minsumet,maxsumet,nsumbins,minsumet,maxsumet);

  htvsavept  = Book2dHist("sumht_vs_avept","Sumht vs Averaged p_{T} of Leading Jets",nsumbinsx,minsumetx,maxsumetx,nsumbinsy,minsumety,maxsumety);

  jetmet_jet0   = Book1dHist("jetmet_jet0","1Jet1Met -- Jet Thresh >0 GeV",nmetbins,minmet,maxmet);
  jetmet_jet50  = Book1dHist("jetmet_jet50","1Jet1Met -- Jet Thresh >50 GeV",nmetbins,minmet,maxmet);
  jetmet_jet80  = Book1dHist("jetmet_jet80","1Jet1Met -- Jet Thresh >80 GeV",nmetbins,minmet,maxmet);
  jetmet_jet140 = Book1dHist("jetmet_jet140","1Jet1Met -- Jet Thresh >150 GeV",nmetbins,minmet,maxmet);
  jetmet_jet180 = Book1dHist("jetmet_jet180","1Jet1Met -- Jet Thresh >180 GeV",nmetbins,minmet,maxmet);

  jetmet_l1_jet0   = Book1dHist("jetmet_l1_jet0","1Jet1Met -- Jet Thresh >0 GeV  --  L1 etm>30",nmetbins,minmet,maxmet);
  jetmet_l1_jet50  = Book1dHist("jetmet_l1_jet50","1Jet1Met -- Jet Thresh >50 GeV --  L1 etm>30",nmetbins,minmet,maxmet);
  jetmet_l1_jet80  = Book1dHist("jetmet_l1_jet80","1Jet1Met -- Jet Thresh >80 GeV --  L1 etm>30",nmetbins,minmet,maxmet);
  jetmet_l1_jet140 = Book1dHist("jetmet_l1_jet140","1Jet1Met -- Jet Thresh >150 GeV --  L1 etm>30",nmetbins,minmet,maxmet);
  jetmet_l1_jet180 = Book1dHist("jetmet_l1_jet180","1Jet1Met -- Jet Thresh >180 GeV --  L1 etm>30",nmetbins,minmet,maxmet);

  const Int_t njetbins=100;
  const Double_t minjet=0.;
  const Double_t maxjet=500.;

  h_caljetpt = Book1dHist("h_caljetpt","Reconstructed Jet p_{T}",njetbins,minjet,maxjet);
  h_genjetpt = Book1dHist("h_genjetpt","Generated Jet p_{T}",njetbins,minjet,maxjet);

  mx_caljetpt_Untriggered= Book1dHist("mx_caljetpt_Untriggered","Maximum Reconstructed Jet P_{T}",njetbins,minjet,maxjet);  
  mx_caljetpt_l15= Book1dHist("mx_caljetpt_l15","Maximum Reconstructed Jet P_{T} -- A_SingleJet15",njetbins,minjet,maxjet);  
  mx_caljetpt_l30= Book1dHist("mx_caljetpt_l30","Maximum Reconstructed Jet P_{T} -- A_SingleJet30",njetbins,minjet,maxjet);  
  mx_caljetpt_l50= Book1dHist("mx_caljetpt_l50","Maximum Reconstructed Jet P_{T} -- A_SingleJet50",njetbins,minjet,maxjet);  
  mx_caljetpt_l70= Book1dHist("mx_caljetpt_l70","Maximum Reconstructed Jet P_{T} -- A_SingleJet70",njetbins,minjet,maxjet);  
  mx_caljetpt_l100= Book1dHist("mx_caljetpt_l100","Maximum Reconstructed Jet P_{T} -- A_SingleJet100",njetbins,minjet,maxjet);  
  mx_caljetpt_l150= Book1dHist("mx_caljetpt_l150","Maximum Reconstructed Jet P_{T} -- A_SingleJet150",njetbins,minjet,maxjet);  
  mx_caljetpt_l200= Book1dHist("mx_caljetpt_l200","Maximum Reconstructed Jet P_{T} -- A_SingleJet200",njetbins,minjet,maxjet);  

  mx_caljetpt_m15= Book1dHist("mx_caljetpt_m15","Maximum Reconstructed Jet P_{T} -- A_SingleJet15",njetbins,minjet,maxjet);  
  mx_caljetpt_m30= Book1dHist("mx_caljetpt_m30","Maximum Reconstructed Jet P_{T} -- A_SingleJet30",njetbins,minjet,maxjet);  
  mx_caljetpt_m50= Book1dHist("mx_caljetpt_m50","Maximum Reconstructed Jet P_{T} -- A_SingleJet50",njetbins,minjet,maxjet);  
  mx_caljetpt_m70= Book1dHist("mx_caljetpt_m70","Maximum Reconstructed Jet P_{T} -- A_SingleJet70",njetbins,minjet,maxjet);  
  mx_caljetpt_m100= Book1dHist("mx_caljetpt_m100","Maximum Reconstructed Jet P_{T} -- A_SingleJet100",njetbins,minjet,maxjet);  
  mx_caljetpt_m150= Book1dHist("mx_caljetpt_m150","Maximum Reconstructed Jet P_{T} -- A_SingleJet150",njetbins,minjet,maxjet);  
  mx_caljetpt_m200= Book1dHist("mx_caljetpt_m200","Maximum Reconstructed Jet P_{T} -- A_SingleJet200",njetbins,minjet,maxjet);  

  mx_genjetpt_Untriggered= Book1dHist("mx_genjetpt_Untriggered","Maximum Generated Jet P_{T}",njetbins,minjet,maxjet);  
  mx_genjetpt_l15= Book1dHist("mx_genjetpt_l15","Maximum Generated Jet P_{T} -- A_SingleJet15",njetbins,minjet,maxjet);  
  mx_genjetpt_l30= Book1dHist("mx_genjetpt_l30","Maximum Generated Jet P_{T} -- A_SingleJet30",njetbins,minjet,maxjet);  
  mx_genjetpt_l50= Book1dHist("mx_genjetpt_l50","Maximum Generated Jet P_{T} -- A_SingleJet50",njetbins,minjet,maxjet);  
  mx_genjetpt_l70= Book1dHist("mx_genjetpt_l70","Maximum Generated Jet P_{T} -- A_SingleJet70",njetbins,minjet,maxjet);  
  mx_genjetpt_l100= Book1dHist("mx_genjetpt_l100","Maximum Generated Jet P_{T} -- A_SingleJet100",njetbins,minjet,maxjet);  
  mx_genjetpt_l150= Book1dHist("mx_genjetpt_l150","Maximum Generated Jet P_{T} -- A_SingleJet150",njetbins,minjet,maxjet);  
  mx_genjetpt_l200= Book1dHist("mx_genjetpt_l200","Maximum Generated Jet P_{T} -- A_SingleJet200",njetbins,minjet,maxjet);  

  mx_genjetpt_m15= Book1dHist("mx_genjetpt_m15","Maximum Generated Jet P_{T} -- A_SingleJet15",njetbins,minjet,maxjet);  
  mx_genjetpt_m30= Book1dHist("mx_genjetpt_m30","Maximum Generated Jet P_{T} -- A_SingleJet30",njetbins,minjet,maxjet);  
  mx_genjetpt_m50= Book1dHist("mx_genjetpt_m50","Maximum Generated Jet P_{T} -- A_SingleJet50",njetbins,minjet,maxjet);  
  mx_genjetpt_m70= Book1dHist("mx_genjetpt_m70","Maximum Generated Jet P_{T} -- A_SingleJet70",njetbins,minjet,maxjet);  
  mx_genjetpt_m100= Book1dHist("mx_genjetpt_m100","Maximum Generated Jet P_{T} -- A_SingleJet100",njetbins,minjet,maxjet);  
  mx_genjetpt_m150= Book1dHist("mx_genjetpt_m150","Maximum Generated Jet P_{T} -- A_SingleJet150",njetbins,minjet,maxjet);  
  mx_genjetpt_m200= Book1dHist("mx_genjetpt_m200","Maximum Generated Jet P_{T} -- A_SingleJet200",njetbins,minjet,maxjet);  

  for (int i=0; i<nhlt; i++ ){
    TString hltname=hltnames[i];
    hname="mx_caljetpt_" + hltname;
    htitle="Maximum Reconstructed Jet P_{T} -- " + hltname;
    h_cal[i]=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

    hname="mx_genjetpt_" + hltname;
    htitle="Maximum Generated Jet P_{T} -- " + hltname;
    h_gen[i]=Book1dHist(hname,htitle,njetbins,minjet,maxjet);
  }

  hname="mx_genjetpt_hlt30";    htitle="Maximum Generated Jet P_{T} -- hlt30";
  mx_genjetpt_hlt30=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt60";    htitle="Maximum Generated Jet P_{T} -- hlt60";
  mx_genjetpt_hlt60=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt65";    htitle="Maximum Generated Jet P_{T} -- hlt65";
  mx_genjetpt_hlt65=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt70";    htitle="Maximum Generated Jet P_{T} -- hlt70";
  mx_genjetpt_hlt70=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt110";    htitle="Maximum Generated Jet P_{T} -- hlt110";
  mx_genjetpt_hlt110=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt150";    htitle="Maximum Generated Jet P_{T} -- hlt150";
  mx_genjetpt_hlt150=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt200";    htitle="Maximum Generated Jet P_{T} -- hlt200";
  mx_genjetpt_hlt200=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  hname="mx_genjetpt_hlt250";    htitle="Maximum Generated Jet P_{T} -- hlt250";
  mx_genjetpt_hlt250=Book1dHist(hname,htitle,njetbins,minjet,maxjet);

  // dijet stuff
  const Int_t ndijetbins=50;
  const Double_t mindijet=0.;
  const Double_t maxdijet=500.;

  hname="h_dijetpt";    htitle="(pt1 + pt2) / 2";
  h_dijetpt=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_dijet30";    htitle="(pt1 + pt2) / 2 -- DiJet 30 Trigger ";
  h_dijetpt_dijet30=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_dijet60";    htitle="(pt1 + pt2) / 2 -- DiJet 60 Trigger ";
  h_dijetpt_dijet60=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_dijet110";    htitle="(pt1 + pt2) / 2 -- DiJet 110 Trigger ";
  h_dijetpt_dijet110=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_dijet150";    htitle="(pt1 + pt2) / 2 -- DiJet 150 Trigger ";
  h_dijetpt_dijet150=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_dijet200";    htitle="(pt1 + pt2) / 2 -- DiJet 200 Trigger ";
  h_dijetpt_dijet200=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);

  hname="h_dijetpt_l15";    htitle="(pt1 + pt2) / 2 -- L1_SingleJet15 ";
  h_dijetpt_l15=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_l30";    htitle="(pt1 + pt2) / 2 -- L1_SingleJet30 ";
  h_dijetpt_l30=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_l70";    htitle="(pt1 + pt2) / 2 -- L1_SingleJet70 ";
  h_dijetpt_l70=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_l100";    htitle="(pt1 + pt2) / 2 -- L1_SingleJet100 ";
  h_dijetpt_l100=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);
  hname="h_dijetpt_l150";    htitle="(pt1 + pt2) / 2 -- L1_SingleJet150 ";
  h_dijetpt_l150=Book1dHist(hname,htitle,ndijetbins,mindijet,maxdijet);


  // level1
  nl1jetc = Book1dHist("h_nl1jetc","Number of L1 Central Jets",5,-0.5,4.5);
  nl1jetf = Book1dHist("h_nl1jetf","Number of L1 Forward Jets",5,-0.5,4.5);
  nl1jett = Book1dHist("h_nl1jett","Number of L1 Tau Jets",5,-0.5,4.5);

  mx_l1cenpt= Book1dHist("mx_l1cenpt","Maximum L1 Central Jet P_{T}",njetbins,minjet,maxjet);
  mx_l1forpt= Book1dHist("mx_l1forpt","Maximum L1 Forward Jet P_{T}",njetbins,minjet,maxjet);
  mx_l1taupt= Book1dHist("mx_l1taupt","Maximum L1 Tau Jet P_{T}",njetbins,minjet,maxjet);

  // hlt overall rates
  hlt_rates=Book1dHist("hlt_rates","HLT Rates",nhlt,0.,nhlt);
  for (int i=0; i<nhlt; i++ ){
    hlt_rates->GetXaxis()->SetBinLabel(i+1,hltnames[i].c_str());
  }
  met_rates=Book1dHist("met_rates","MET Rates",9,0.,9.);
  met_rates->GetXaxis()->SetBinLabel(1,"A_ETM20");
  met_rates->GetXaxis()->SetBinLabel(2,"A_ETM30");
  met_rates->GetXaxis()->SetBinLabel(3,"A_ETM40");
  met_rates->GetXaxis()->SetBinLabel(4,"A_ETM50");
  met_rates->GetXaxis()->SetBinLabel(5,"HLT1MET");
  met_rates->GetXaxis()->SetBinLabel(6,"A_ETM30 PT>50");
  met_rates->GetXaxis()->SetBinLabel(7,"A_ETM30 PT>60");
  met_rates->GetXaxis()->SetBinLabel(8,"A_ETM30 PT>70");
  met_rates->GetXaxis()->SetBinLabel(9,"A_ETM30 PT>75");
}

void loadChain(int istart, int ifin){

  for (int ip = istart; ip != ifin; ++ip){

    TString filename=ProcFil.at(ip);
    cout << "Adding to chain: " << filename <<endl;
    TabChain[ip] = new TChain("HltTree");
    TabChain[ip]->Add(filename);

    TabChain[ip]->SetBranchStatus("*",0);
    TabChain[ip]->SetBranchStatus("recoMetCal",1);
    TabChain[ip]->SetBranchStatus("recoMetCalPhi",1);
    TabChain[ip]->SetBranchStatus("recoMetCalSum",1);

    TabChain[ip]->SetBranchStatus("recoMetGen",1);

    TabChain[ip]->SetBranchStatus("recoHTCal",1);
    TabChain[ip]->SetBranchStatus("recoHTCalSum",1);

    TabChain[ip]->SetBranchStatus("recoJetCalPt",1);
    TabChain[ip]->SetBranchStatus("recoJetCalEta",1);
    TabChain[ip]->SetBranchStatus("NrecoJetCal",1);
    
    TabChain[ip]->SetBranchStatus("recoJetGenPt",1);
    TabChain[ip]->SetBranchStatus("recoJetGenEta",1);
    TabChain[ip]->SetBranchStatus("NrecoJetGen",1);
        
    TabChain[ip]->SetBranchStatus("NL1ForJet",1);
    TabChain[ip]->SetBranchStatus("NL1CenJet",1);
    TabChain[ip]->SetBranchStatus("NL1Tau",1);
    
    TabChain[ip]->SetBranchStatus("L1CenJetEt",1);
    TabChain[ip]->SetBranchStatus("L1CenJetEta",1);
    TabChain[ip]->SetBranchStatus("L1CenJetPhi",1);
    
    TabChain[ip]->SetBranchStatus("L1ForJetEt",1);
    TabChain[ip]->SetBranchStatus("L1ForJetEta",1);
    TabChain[ip]->SetBranchStatus("L1ForJetPhi",1);
    
    TabChain[ip]->SetBranchStatus("L1TauEt",1);
    TabChain[ip]->SetBranchStatus("L1TauEta",1);
    TabChain[ip]->SetBranchStatus("L1TauPhi",1);
    
    TabChain[ip]->SetBranchStatus("L1Met",1);
    TabChain[ip]->SetBranchStatus("L1MetPhi",1);
    TabChain[ip]->SetBranchStatus("L1MetTot",1);
    TabChain[ip]->SetBranchStatus("L1MetHad",1);
    
    //Trigger bit branches

    TabChain[ip]->SetBranchStatus("L1_SingleJet15",1);  
    TabChain[ip]->SetBranchStatus("L1_SingleJet30",1);  
    TabChain[ip]->SetBranchStatus("L1_SingleJet50",1);  
    TabChain[ip]->SetBranchStatus("L1_SingleJet70",1);  
    TabChain[ip]->SetBranchStatus("L1_SingleJet100",1);  
    TabChain[ip]->SetBranchStatus("L1_SingleJet150",1);  
    TabChain[ip]->SetBranchStatus("L1_SingleJet200",1);  
    TabChain[ip]->SetBranchStatus("L1_ETM20",1);  
    TabChain[ip]->SetBranchStatus("L1_ETM30",1);  
    TabChain[ip]->SetBranchStatus("L1_ETM40",1);  
    TabChain[ip]->SetBranchStatus("L1_ETM50",1);  
    TabChain[ip]->SetBranchStatus("L1_ETM60",1);  

    //Monte Carlo branches
    TabChain[ip]->SetBranchStatus("MCPtHat",1);  

    //HLT Trigger bit branches
    int nhlt=hltnames.size();
    for (int i=0; i<nhlt; i++ ){
      TString hltname=hltnames[i];
      //cout << "Adding HLT trigger: " << hltname << endl;
      TabChain[ip]->SetBranchStatus(hltname,1);  
    }
  }
}

void setBranchAdds(int ip){

  TabChain[ip]->SetBranchAddress("recoMetCal",&metpt);
  TabChain[ip]->SetBranchAddress("recoMetCalPhi",&metphi);
  TabChain[ip]->SetBranchAddress("recoMetCalSum",&metsum);

  TabChain[ip]->SetBranchAddress("recoMetGen",&genmetpt);

  TabChain[ip]->SetBranchAddress("recoHTCal",&htpt);
  TabChain[ip]->SetBranchAddress("recoHTCalSum",&htsum);

  TabChain[ip]->SetBranchAddress("NrecoJetCal",&njetcal);
  TabChain[ip]->SetBranchAddress("recoJetCalPt",caljetpt);
  TabChain[ip]->SetBranchAddress("recoJetCalEta",caljeteta);
  
  TabChain[ip]->SetBranchAddress("NrecoJetGen",&njetgen);
  TabChain[ip]->SetBranchAddress("recoJetGenPt",genjetpt);
  TabChain[ip]->SetBranchAddress("recoJetGenEta",genjeteta);
  
  TabChain[ip]->SetBranchAddress("NL1CenJet",&nl1cenjet);
  TabChain[ip]->SetBranchAddress("L1CenJetEt",l1cenjetpt);
  
  TabChain[ip]->SetBranchAddress("NL1ForJet",&nl1forjet);
  TabChain[ip]->SetBranchAddress("L1ForJetEt",l1forjetpt);
  
  TabChain[ip]->SetBranchAddress("NL1Tau",&nl1taujet);
  TabChain[ip]->SetBranchAddress("L1TauEt",l1taujetpt);

  TabChain[ip]->SetBranchAddress("L1Met",&l1met);
  TabChain[ip]->SetBranchAddress("L1MetPhi",&l1metphi);
  TabChain[ip]->SetBranchAddress("L1MetTot",&l1mettot);
  TabChain[ip]->SetBranchAddress("L1MetHad",&l1mettothad);
  

  TabChain[ip]->SetBranchAddress("L1_SingleJet15",&i_l1bit15);
  TabChain[ip]->SetBranchAddress("L1_SingleJet30",&i_l1bit30);
  TabChain[ip]->SetBranchAddress("L1_SingleJet50",&i_l1bit50);
  TabChain[ip]->SetBranchAddress("L1_SingleJet70",&i_l1bit70);
  TabChain[ip]->SetBranchAddress("L1_SingleJet100",&i_l1bit100);
  TabChain[ip]->SetBranchAddress("L1_SingleJet150",&i_l1bit150);
  TabChain[ip]->SetBranchAddress("L1_SingleJet200",&i_l1bit200);
  
  TabChain[ip]->SetBranchAddress("L1_ETM20",&i_etm20);  
  TabChain[ip]->SetBranchAddress("L1_ETM30",&i_etm30);  
  TabChain[ip]->SetBranchAddress("L1_ETM40",&i_etm40);  
  TabChain[ip]->SetBranchAddress("L1_ETM50",&i_etm50);  
  TabChain[ip]->SetBranchAddress("L1_ETM60",&i_etm60);  
  
  TabChain[ip]->SetBranchAddress("L1_HTT100",&i_l1htt100);
  TabChain[ip]->SetBranchAddress("L1_HTT200",&i_l1htt200);
  TabChain[ip]->SetBranchAddress("L1_HTT250",&i_l1htt250);
  TabChain[ip]->SetBranchAddress("L1_HTT300",&i_l1htt300);
  TabChain[ip]->SetBranchAddress("L1_HTT400",&i_l1htt400);
  TabChain[ip]->SetBranchAddress("L1_HTT500",&i_l1htt500);

  TabChain[ip]->SetBranchAddress("MCPtHat",&mcpthat);


  // Get HLT Trigger Info;
  for (int i=0; i<nhlt; i++ ){
    string hltname=hltnames[i];
    TabChain[ip]->SetBranchAddress(hltname.c_str(),&hltTriggers[hltname]);  
  }
}

double getXSWeight(const float pthat,const int mcSample){

  double wt=0.;
  if (pthat<0.001)return wt;

  if (mcSample == 0){ // Ordinary QCD binned samples
      if (pthat < 15.) {
	cout << "%getXSWeight--Error pthat < 15. and mcSample=0 -- Setting cross section weight to 0" << endl;
	wt=0.0;
      }
      else if (pthat < 20. ) wt=1.46E9;
      else if (pthat < 30. ) wt=6.32E8;
      else if (pthat < 50. ) wt=1.63E8;
      else if (pthat < 80. ) wt=2.16E7;
      else if (pthat < 120. ) wt=3.08E6;
      else if (pthat < 170.) wt=4.94E5;
      else if (pthat < 230.) wt=1.01E5;
      else if (pthat < 300.) wt=2.45E4;
      else if (pthat < 380.) wt=6.24E3;
      else if (pthat < 470.) wt=1.78E3;
      else if (pthat < 600.) wt=6.83E2;
      else if (pthat < 800.) wt=2.04E2; 
      else if (pthat < 1000.) wt=3.51E1;
      else if (pthat < 1400.) wt=1.09E1;
      else if (pthat < 1800.) wt=1.06;
      else if (pthat < 2200.) wt=1.45E-1;
      else if (pthat < 2600.) wt=2.38E-2;
      else if (pthat < 3000.) wt=4.29E-3;
      else if (pthat < 3500.) wt=8.44E-4;
      else wt=1.08E-4;
  } else if (mcSample == 1){ // QCD 0-15 sample -- eliminated events with pthat > 15 GeV/c
      wt=5.52E10;
    if (pthat > 15.) {
      cout << "%getXSWeight--QCD 0-15 sample with pthat > 15. -- Setting cross section weight to 0" << endl;
      wt=0.0;
    }
  } else {
    cout << "%getXSWeight--Undefined MC Sample type -- Setting cross section weight to 0" << endl;
    wt=0.0;
  }
  return wt;
}

double getSkimEff(const float pthat,const int mcSample){

  double wt=0.;
  if (pthat<0.001)return wt;

  if (mcSample <= 1){ // QCD binned samples
    if (pthat < 15.) wt=0.00283936;
    else if (pthat < 20. ) wt=0.0214774;
    else if (pthat < 30. ) wt=0.0587037;
    else if (pthat < 50. ) wt=0.219188;
    else if (pthat < 80. ) wt=0.640025;
    else if (pthat < 120. ) wt=0.959605;
    else if (pthat < 170.) wt=0.998667;
    else if (pthat < 230.) wt=1.;
    else if (pthat < 300.) wt=1.;
    else if (pthat < 380.) wt=1.;
    else if (pthat < 470.) wt=1.;
    else if (pthat < 600.) wt=1.;
    else if (pthat < 800.) wt=1.; 
    else if (pthat < 1000.) wt=1.;
    else if (pthat < 1400.) wt=1.;
    else if (pthat < 1800.) wt=1.;
    else if (pthat < 2200.) wt=1.;
    else if (pthat < 2600.) wt=1.;
    else if (pthat < 3000.) wt=1.;
    else if (pthat < 3500.) wt=1.;
    else wt=1.;
  } else {
    cout << "%getSkimEff--Undefined MC Sample type -- Setting SkimEff to 0" << endl;
    wt=0.0;
  }
  return wt;
}

void jetplusmetRates(double metpt,double maxjetpt,int itrig,double wt){
  //cout << "Here" << endl;

  bool jet50=false,jet80=false,jet140=false,jet180=false;

  if (maxjetpt>50) jet50=true;
  if (maxjetpt>80) jet80=true;
  if (maxjetpt>140) jet140=true;
  if (maxjetpt>180) jet180=true;

  jetmet_jet0->Fill(metpt,wt);
  if (jet50) jetmet_jet50->Fill(metpt,wt);
  if (jet80) jetmet_jet80->Fill(metpt,wt);
  if (jet140) jetmet_jet140->Fill(metpt,wt);
  if (jet180) jetmet_jet180->Fill(metpt,wt);

  if (itrig>0 ){ // no level-1 cuts
    jetmet_l1_jet0->Fill(metpt,wt);
    if (jet50) jetmet_l1_jet50->Fill(metpt,wt);
    if (jet80) jetmet_l1_jet80->Fill(metpt,wt);
    if (jet140) jetmet_l1_jet140->Fill(metpt,wt);
    if (jet180) jetmet_l1_jet180->Fill(metpt,wt);
  }
}

int readFilelist(std::string filelist){
  int istat=0; //assume success

  std::string CommentLine="#"; // treat lines that begin with "#" as comment lines
  ifstream in;

  string input_line,word,filename;

  cout << "%readFilelist--Reading: " << filelist << endl;
  in.open(filelist.c_str());
  if (!in){
    cerr << "%Could not open file: " << filelist << endl;
    cerr << "%Terminating program" << endl;
    return -1;
  }

  //  while (1) {
  while (getline(in,input_line)){

    if (!in.good()) break;

    uint p1 = input_line.find (CommentLine,0);
    if ( p1 == std::string::npos){
      istringstream stream(input_line);
      std::vector<string> elements;

      while (stream >> word) {
	elements.push_back(word);
      }

      if (elements.size() >= 2){
	ProcFil.push_back(elements[0]);
	mcSample_bin.push_back(atoi(elements[1].c_str()));
      }else if (elements.size() == 0){
	cout << "%readFilelist--Blank line encountered -- Treating as EOF" << endl;
	return 0;
      }else {
	cout << "%readFilelist--Error parsing input filelist" << endl;
	return -1;
      }
    }
  }
  return istat;

}


float deltaPhi (float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}
