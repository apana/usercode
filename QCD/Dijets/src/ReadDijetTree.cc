#include <iostream>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <TFile.h>
#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TTree.h>
#include <TSystem.h>
#include "ReadDijetTree.h"

//---------------------------------------------------------------------
//                  HISTORY
//  V0 -- original from Kostas Kousouris
//  V1 -- Read trees with leading 2 jet four-vector information included
//  V2 -- Read trees storing top 10 jets in event
//---------------------------------------------------------------------

using namespace std;

//int main()
int main(int argc, char *argv[])
{
  string inputFile="in.root";
  if (argc>1) {
    inputFile=argv[1];
    cout << "Input file:  "<< inputFile << "\n";
  } else {
    cout << "Please supply input filename" << endl;
    return 1;
  }

  // read input control cards
  const char *inpc = gSystem->Getenv("INPC");
  if (inpc) { 
    readInputCards(inpc);
  } else {
    cout << "INPC environment variable not set!!" << endl;
    return 1;
  }
  //TFile *inf = new TFile(inputFile.c_str(),"R");
  TFile *inf = TFile::Open(inputFile.c_str(),"R");

  string outputFile=createOutputFileName();
  cout << "Output histograms written to: " << outputFile << "\n" << endl;
  TFile *outf = new TFile(outputFile.c_str(),"RECREATE");

  hNjets        = bookTH1F("Njets","Number of Saved Jets in event",11, -0.5,10.5); 

  TH1F *hMass         = bookTH1F("JetMass","JetMass",300, 0,3000); 
  TH1F *hMassTrg      = bookTH1F("JetMassTrg","JetMassTrg -- All trigger cuts, Mass Cut",300, 0,3000);
  TH1F *hMassComp     = bookTH1F("JetMassComp","JetMassComp -- All trigger cuts, Mass Cut",300, 0,3000);

  TH1F *hMass1         = bookTH1F("JetMass1","JetMass",NmassBins1,massBins1); 
  TH1F *hMassTrg1      = bookTH1F("JetMassTrg1","JetMassTrg -- All trigger cuts, Mass Cut",NmassBins1,massBins1);
  TH1F *hMassComp1     = bookTH1F("JetMassComp1","JetMassComp -- All trigger cuts, Mass Cut",NmassBins1,massBins1);

  //TH1F *hMass2         = bookTH1F("JetMass2","JetMass",NmassBins,massBins); 
  //TH1F *hMassTrg2      = bookTH1F("JetMassTrg2","JetMassTrg -- All trigger cuts, Mass Cut",NmassBins2,massBins2);
  //TH1F *hMassComp2     = bookTH1F("JetMassComp2","JetMassComp -- All trigger cuts, Mass Cut",NmassBins2,massBins2);

  TH1F *hEta1         = bookTH1F("Jet1Eta","Jet1Eta",100,-5,5);
  TH1F *hEta2         = bookTH1F("Jet2Eta","Jet2Eta",100,-5,5);
  TH1F *hEtaBoost     = bookTH1F("JetEtaBoost","JetEtaBoost",100,-5,5);
  TH1F *hEtaStar      = bookTH1F("JetEtaStar","JetEtaStar",100,-5,5);
  TH1F *hRap1         = bookTH1F("Jet1Rap","Jet1Rap",100,-5,5);
  TH1F *hRap2         = bookTH1F("Jet2Rap","Jet2Rap",100,-5,5);

  TH1F *hPt1          = bookTH1F("Jet1Pt","Jet1Pt",1200,0,6000);
  TH1F *hPt2          = bookTH1F("Jet2Pt","Jet2Pt",1200,0,6000);

  TH1F *hPtMax        = bookTH1F("JetPtMax","JetPtMax",1200,0,6000);
  TH1F *hYMax1        = bookTH1F("JetYMax1","JetYMax pt > 60 ",100,-5,5);
  TH1F *hYMax2        = bookTH1F("JetYMax2","JetYMax pt > 200 ",100,-5,5);

  TH1F *hPhi1         = bookTH1F("Jet1Phi","Jet1Phi",100,0,6.28318);
  TH1F *hPhi2         = bookTH1F("Jet2Phi","Jet2Phi",100,0,6.28318);
  TH1F *hdPhi         = bookTH1F("JetdPhi","JetdPhi",100,0,6.28318);
  TH1F *hCosThetaStar = bookTH1F("JetCosThetaStar","JetCosThetaStar",100,-1,1);

  TH1F *hChi          = bookTH1F("JetX","JetX",100,0,50);
  //TH1F *hChiTrg       = bookTH1F("JetXTrg","JetX -- All trigger cuts, Mass Cut",20,0,20);
  //TH1F *hChiComp      = bookTH1F("JetXComp","JetX -- All trigger cuts, Mass Cut",20,0,20);
  TH1F *hChiTrg       = bookTH1F("JetXTrg","JetX -- All trigger cuts, Mass Cut",NchiBins1,chiBins1);
  TH1F *hChiComp1     = bookTH1F("JetXComp1","JetX -- All trigger cuts, Mass Cut",NchiBins1,chiBins1);
  TH1F *hChiComp2     = bookTH1F("JetXComp2","JetX -- All trigger cuts, Mass Cut",NchiBins2,chiBins2);
  TH1F *hChiComp3     = bookTH1F("JetXComp3","JetX -- All trigger cuts, Mass Cut",NchiBins3,chiBins3);
  TH1F *hChiComp4     = bookTH1F("JetXComp4","JetX -- All trigger cuts, Mass Cut",NchiBins4,chiBins4);


  TH2F *hPhi12        = new TH2F("Jet12Phi","Jet12Phi",100,0,6.28318,100,0,6.28318);
  TH2F *hEta12        = new TH2F("Jet12Eta","Jet12Eta",100,-5,5,100,-5,5);
  TH2F *hCMEta12      = new TH2F("Jet12CMEta","JetCMEta",100,-5,5,100,-5,5);   
  //ccla
  TH2F *hMassVsChi    = new TH2F("MassVsChi","MassVsChi",30,0,30,140,0,1400);

  if (smrGenEng){
    Int_t nres=50;
    Double_t xmin=0.,xmax=2.0;
    hResB1          = bookTH1F("ResB1","Smeared GenJet p_{T} / GenJet p_{T} --  80 < p_{T} < 100 ",nres,xmin,xmax);
    hResB2          = bookTH1F("ResB2","Smeared GenJet p_{T} / GenJet p_{T} -- 200 < p_{T} < 250 ",nres,xmin,xmax);
    hResB3          = bookTH1F("ResB3","Smeared GenJet p_{T} / GenJet p_{T} -- 400 < p_{T} < 450 ",nres,xmin,xmax);
    hResB4          = bookTH1F("ResB4","Smeared GenJet p_{T} / GenJet p_{T} -- 900 < p_{T} < 1100",nres,xmin,xmax);

    hResE1          = bookTH1F("ResE1","Smeared GenJet p_{T} / GenJet p_{T} --  80 < p_{T} < 100 ",nres,xmin,xmax);
    hResE2          = bookTH1F("ResE2","Smeared GenJet p_{T} / GenJet p_{T} -- 200 < p_{T} < 250 ",nres,xmin,xmax);
    hResE3          = bookTH1F("ResE3","Smeared GenJet p_{T} / GenJet p_{T} -- 400 < p_{T} < 450 ",nres,xmin,xmax);
    hResE4          = bookTH1F("ResE4","Smeared GenJet p_{T} / GenJet p_{T} -- 900 < p_{T} < 1100",nres,xmin,xmax);

    hResF1          = bookTH1F("ResF1","Smeared GenJet p_{T} / GenJet p_{T} --  80 < p_{T} < 100 ",nres,xmin,xmax);
    hResF2          = bookTH1F("ResF2","Smeared GenJet p_{T} / GenJet p_{T} -- 200 < p_{T} < 250 ",nres,xmin,xmax);
    hResF3          = bookTH1F("ResF3","Smeared GenJet p_{T} / GenJet p_{T} -- 400 < p_{T} < 450 ",nres,xmin,xmax);
    hResF4          = bookTH1F("ResF4","Smeared GenJet p_{T} / GenJet p_{T} -- 900 < p_{T} < 1100",nres,xmin,xmax);
  }
  bool cut_eta,cut_cosThetaStar,cut_mass;

  int i;
  TTree *tr = (TTree*)inf->Get("DijetTree");

  TBranch *brNjets = (TBranch*)tr->GetBranch("njets");
  brNjets->SetAddress(&njets);
  TBranch *brPx = (TBranch*)tr->GetBranch("px");
  brPx->SetAddress(&px);
  TBranch *brPy = (TBranch*)tr->GetBranch("py");
  brPy->SetAddress(&py);
  TBranch *brPz = (TBranch*)tr->GetBranch("pz");
  brPz->SetAddress(&pz);
  TBranch *brE = (TBranch*)tr->GetBranch("e");
  brE->SetAddress(&e);
  TBranch *brMass = (TBranch*)tr->GetBranch("mass");
  brMass->SetAddress(&m);
  
  //setup the trigger cuts
  fillTriggerInfo();
  TRIGGER_ET_CUT=getTriggerEtCut(triggerStream);
  TRIGGER_PRESCALE=getTriggerPrescl(triggerStream);
  TRIGGER_XS_WEIGHT=getTriggerXSWeight(triggerStream);
  TRIGGER_NEVT=getTriggerNevt(triggerStream);
  TRIGGER_MINMASS=getTriggerMinMass(triggerStream);
  TRIGGER_MAXMASS=getTriggerMaxMass(triggerStream);

  if (fabs(ETA_CUT-3.0)<0.01){
    ETASTAR_CUT=1.5;
    ETABOOST_CUT=1.5;
  }else if (fabs(ETA_CUT-2.0)<0.01){
    ETASTAR_CUT=1.5;
    ETABOOST_CUT=0.5;
  }else{
    cout << "Illegal value of ETA_CUT -- Ending program execution" << endl;
    return 1;
  }
    
  //cout << "Number of HLT triggers: " << hltnames.size() << endl;
  cout << "TRIGGER_ET_CUT: " << TRIGGER_ET_CUT << "\n";
  cout << "TRIGGER_MINMASS: " << TRIGGER_MINMASS << "\n" ;
  cout << "TRIGGER_MAXMASS: " << TRIGGER_MAXMASS << "\n" ;
  cout << "TRIGGER_PRESCALE: " << TRIGGER_PRESCALE << "\n" ;
  cout << "TRIGGER_NEVT:     " << TRIGGER_NEVT << "\n" ;
  cout << "TRIGGER_XS_WEIGHT: " << TRIGGER_XS_WEIGHT << "\n";

  cout << "ETA 1 and ETA 2 MAX: " << ETA_CUT << "\n";
  cout << "ETABOOST_CUT: " << ETABOOST_CUT << "\n";
  cout << "ETASTAR_CUT: " << ETASTAR_CUT << "\n" << endl;

  if (abs(engShift) > 0){
    cout << "Shifting Jet energies using option: " << engShift << "\n";
  }

  if (do1pb){
    cout << "Running 1pb-1 luminosity" << "\n";
  }

  //  if (nevt<0) nevt=tr->GetEntries();
  if (nevt<0) nevt=int(TRIGGER_NEVT);
  cout << "\nNumber of events to process: " << nevt << endl;
  if (nevt != tr->GetEntries()){
    cout << "\t!Warning: Number of events to process does not match number of events on tree" << endl;
    cout << "Number of events on tree: " << tr->GetEntries() << endl;
    if (nevt > tr->GetEntries()) nevt=tr->GetEntries();
  }
  cout << " " << endl;

  double wght=1.,xswght=1.,twght=1.;
  if (weightXS){
    twght=TRIGGER_PRESCALE;
    xswght=TRIGGER_XS_WEIGHT/nevt;
  }
  for(i=0;i<nevt;i++)
    {
      if (i % 1000000 == 0)
        cout<< "Processed: "<<i<< " events so far" << endl;
      tr->GetEntry(i);
     
      wght=twght;
      // calculate jet and dijet kinematic variables
      loadVar();

      cut_eta = (fabs(eta1)<ETA_CUT && fabs(eta2)<ETA_CUT);
      cut_mass = (mass_>TRIGGER_MINMASS && mass_<TRIGGER_MAXMASS);

      double ptmax = pt1>pt2 ? pt1 : pt2;
      double ymax = pt1>pt2 ? rap1 : rap2;

      if (adjustSpectrum){
	double xsD0 =xs(mass_,parsD0);
	double xsCMS=xs(mass_,parsCMS);
	wght=wght*xsD0/xsCMS;
	//cout << "mass: " << mass_ << " fact: " << xsD0/xsCMS << endl;
      }
      if (fabs(yBoost_)<ETABOOST_CUT && 
	  fabs(yStar_)<ETASTAR_CUT &&
	  ptmax > TRIGGER_ET_CUT &&
	  chi_<CHI_CUT
	  ){
	//if (fabs(etaBoost)<1.5 && fabs(etaStar)<1.5 && chi<20){


	hMass->Fill(mass_,wght);
	hMass1->Fill(mass_,wght);
	hEtaStar->Fill(etaStar_,wght);
	hEtaBoost->Fill(etaBoost_,wght);

 	hMassTrg->Fill(mass_,wght);
 	hMassTrg1->Fill(mass_,wght);
	hChiTrg->Fill(chi_,wght);

	if (cut_mass){
	  hMassComp->Fill(mass_,wght);
	  hMassComp1->Fill(mass_,wght);

	  hChiComp1->Fill(chi_,wght);
	  hChiComp2->Fill(chi_,wght);
	  hChiComp3->Fill(chi_,wght);
	  hChiComp4->Fill(chi_,wght);
	}

	hEta1->Fill(eta1,wght);
	hEta2->Fill(eta2,wght);	
	hRap1->Fill(rap1,wght);
	hRap2->Fill(rap2,wght);
	hPt1->Fill(pt1,wght);
	hPt2->Fill(pt2,wght);
      }

      hPhi1->Fill(phi1,wght);
      hPhi2->Fill(phi2,wght);
      hdPhi->Fill(fabs(dphi_),wght);
      hChi->Fill(chi_,wght);
      hCosThetaStar->Fill(cosThetaStar_,wght);
      hPhi12->Fill(phi1,phi2,wght);
      hEta12->Fill(eta1,eta2,wght);
      hCMEta12->Fill(etaBoost_,etaStar_,wght);  


      hPtMax->Fill(ptmax,xswght);
      if (ptmax > 60){
	hYMax1->Fill(ymax,xswght);
      }
      if (ptmax > 200){
	hYMax2->Fill(ymax,xswght);
      }
      
      if (ptmax > TRIGGER_ET_CUT){
	if (abs(pt1-pt2)/pt1 < 0.1)
	  hMassVsChi->Fill(chi_,mass_,wght);  
      }
      
    }
  outf->Write();
  outf->Close();
}

void loadVar(){

  if (njets > 10 ) {
    cout << "%!Problem with jet size -- njets= " << njets << endl;
    njets=10;
  }

  if (njets<2) return;
  hNjets->Fill(njets,1.);

  vector<TLorentzVector> vjetP4; 
  for (int ijet=0; ijet<njets; ijet++){

    sclFactor=1.;
    if (abs(engShift) > 0){
      double pt;
      pt=sqrt(pow(px[ijet],2)+pow(py[ijet],2));
      sclFactor=scaleJetEnergies(pt);
    }

    /*
    jetP4[ijet].SetPx(px[ijet]*sclFactor);
    jetP4[ijet].SetPy(py[ijet]*sclFactor);
    jetP4[ijet].SetPz(pz[ijet]*sclFactor);
    jetP4[ijet].SetE(  e[ijet]*sclFactor);
    */
    TLorentzVector jetP4;
    jetP4.SetPx(px[ijet]*sclFactor);
    jetP4.SetPy(py[ijet]*sclFactor);
    jetP4.SetPz(pz[ijet]*sclFactor);
    jetP4.SetE(  e[ijet]*sclFactor);

    if (smrGenEng){
      double pt =jetP4.Pt();
      double rap=jetP4.Rapidity();
      double sigma=jetResolution(pt,rap);
      double sclFactor=Trand.Gaus(1,sigma);

      jetP4.SetPx(px[ijet]*sclFactor);
      jetP4.SetPy(py[ijet]*sclFactor);
      jetP4.SetPz(pz[ijet]*sclFactor);
      jetP4.SetE(  e[ijet]*sclFactor);

      fillJetResolutionPlots(jetP4.Pt(),jetP4.Rapidity(),pt,rap);
    }
    vjetP4.push_back(jetP4);
  }// end loop over input jets

  vjetP4=sortJets(vjetP4);
  //for (int ijet=0; ijet<njets; ijet++){
  //  cout << ijet << " " << vjetP4.at(ijet).Pt() << endl;
  //}
    
  pt1=vjetP4.at(0).Pt();
  pt2=vjetP4.at(1).Pt();
  rap1=vjetP4.at(0).Rapidity();
  rap2=vjetP4.at(1).Rapidity();
  eta1=vjetP4.at(0).Eta();
  eta2=vjetP4.at(1).Eta();  
  phi1=Phi_0_2pi(vjetP4.at(0).Phi());
  phi2=Phi_0_2pi(vjetP4.at(1).Phi());
  
  dijetP4=vjetP4.at(0)+vjetP4.at(1);
  
  mass_= sqrt(dijetP4*dijetP4);
  rapidity_=dijetP4.Rapidity();
  etaBoost_ = 0.5*(eta1+eta2);
  etaStar_ = 0.5*(eta1-eta2); 
  yBoost_ = 0.5*(rap1+rap2);
  yStar_ = 0.5*fabs(rap1-rap2); 
  dphi_ = phi1-phi2; 
  cosThetaStar_ = tanh(yStar_); 
  chi_ = exp(2*yStar_);

  //if (abs(mass_-m)>0.001 && fabs(1.-sclFactor)<0.001) 
  //  cout << "Mass diff: " << mass_-m << " " << mass_ << " " << m << endl;

  //cout << abs(phi2 - phi[1]) << endl;


}

double jetResolution(const double pt,const double y){
  //
  // returns fractional jet resol sigma/pt
  // from Agata's MC Truth resolution plots
  //
  double res=0.;
  double a,b,c;
  if (abs(y)<1.3){ // barrel
    a=4.0;
    b=1.4;
    c=0.046;    
  } else if (abs(y)< 3){ // endcap
    a=0.;
    b=1.4;
    c=0.051; 
  } else { // forward
    a=4.0;
    b=0;
    c=0.11; 
  }

  double aa=pow(a/pt,2);
  double bb=pow(b/sqrt(pt),2);
  double cc=pow(c,2);

  res=sqrt(aa+bb+cc);
  //res=res*3.;
  return res;

}

void fillJetResolutionPlots(const double pt, const double rap, const double ptgen, const double rapgen){

  if (abs(rapgen)<1.3){
    if (ptgen > 80 && ptgen < 100) {
      hResB1->Fill(pt/ptgen,1.);
    } else if (ptgen > 200 && ptgen < 250) {
      hResB2->Fill(pt/ptgen,1.);
    } else if (ptgen > 400 && ptgen < 450) {
      hResB3->Fill(pt/ptgen,1.);
    } else if (ptgen > 900 && ptgen < 1100) {
      hResB4->Fill(pt/ptgen,1.);
    }
  } else if (abs(rapgen)<3){
    if (ptgen > 80 && ptgen < 100) {
      hResE1->Fill(pt/ptgen,1.);
    } else if (ptgen > 200 && ptgen < 250) {
      hResE2->Fill(pt/ptgen,1.);
    } else if (ptgen > 400 && ptgen < 450) {
      hResE3->Fill(pt/ptgen,1.);
    } else if (ptgen > 900 && ptgen < 1100) {
      hResE4->Fill(pt/ptgen,1.);
    }
  } else{
    if (ptgen > 80 && ptgen < 100) {
      hResF1->Fill(pt/ptgen,1.);
    } else if (ptgen > 200 && ptgen < 250) {
      hResF2->Fill(pt/ptgen,1.);
    } else if (ptgen > 400 && ptgen < 450) {
      hResF3->Fill(pt/ptgen,1.);
    } else if (ptgen > 900 && ptgen < 1100) {
      hResF4->Fill(pt/ptgen,1.);
    }
  }
}

double xs(const double xx, const double *pars)
{

  /*
  double value=1;
  float N=1.e13;
  //float a=-4;
  float a=-5;
  float b=1.;
  */

  double N=pars[0];
  double a=pars[1];
  double b=pars[2];

  const double roots=14000.;

  double pt=xx;
  double value=N*pow(pt,a)*pow((1-2*pt/roots),b);

  return value;
}

double scaleJetEnergies(const double pt){

  double fact=1;
  double eScl=0.;
  if (abs(engShift) == 1){ // 10% energy scale uncertainly
    eScl=0.1;
    fact = 1.+eScl;
    if (engShift < 0) fact = 1.-eScl;
  }else if(abs(engShift) == 2){ // linear dependency 10% at 100GeV, 20% at 2 TeV
    eScl=0.1;
    if (pt>100.) eScl=5.2632e-05*pt+0.094737;
    fact = 1.+eScl;
    if (engShift < 0) fact = 1.-eScl;
  }
  return fact;
}

vector <TLorentzVector> sortJets(const vector<TLorentzVector> p4){ 
  vector <TLorentzVector> p4sorted;

  p4sorted=p4;
  std::sort(p4sorted.begin(),p4sorted.end(),PtGreater());
  //cout << "Size of jet vector: " << p4sorted.size() << endl;

  //int njets=p4sorted.size();
  //for (int ijet=0; ijet<njets; ijet++){
  //  cout << ijet << " " << p4.at(ijet).Pt() << " " << p4sorted.at(ijet).Pt() << endl;
  // }
  return p4sorted;
}
