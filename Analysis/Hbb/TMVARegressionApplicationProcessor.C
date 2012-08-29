/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVARegressionApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#endif

using namespace TMVA;

TString ANALYSIS = "Dijet";
//TString ANALYSIS = "Subjet";

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


double ptfatHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1 , double e1, double pt2, double eta2, double phi2, double e2, double oldpt0, double oldpt1, double oldpt2) {
  TLorentzVector j0, j1 , j2, H;

  j0.SetPtEtaPhiE(pt0,eta0,phi0,e0*pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1,e1*pt1/oldpt1);
  j2.SetPtEtaPhiE(pt2,eta2,phi2,e2*pt2/oldpt2);
  H = j0 +j1 +j2;
  return H.Pt(); 
}

double massHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1 , double e1, double oldpt0, double oldpt1) {
  TLorentzVector j0, j1 , H;

  j0.SetPtEtaPhiE(pt0,eta0,phi0,e0*pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1,e1*pt1/oldpt1);
  H = j0 +j1;
  return fabs(H.M());
}

double massfatHiggs(double pt0, double eta0, double phi0, double e0, double pt1, double eta1, double phi1 , double e1, double pt2, double eta2, double phi2, double e2, double oldpt0, double oldpt1, double oldpt2) {
  TLorentzVector j0, j1 , j2, H;

  j0.SetPtEtaPhiE(pt0,eta0,phi0,e0*pt0/oldpt0);
  j1.SetPtEtaPhiE(pt1,eta1,phi1,e1*pt1/oldpt1);
  j2.SetPtEtaPhiE(pt2,eta2,phi2,e2*pt2/oldpt2);
  H = j0 +j1 +j2;
  return fabs(H.M()); 
}

float resolutionBias(float eta){
  if(eta< 1.1) return 0.05;
  if(eta< 2.5) return 0.10;
  if(eta< 5) return 0.30;
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


double evalJERBias( double ptreco, double ptgen, double eta){
  double cor =1;   
  if ((fabs(ptreco - ptgen)/ ptreco)<0.5) { //Limit the effect to the core 
    cor = (ptreco +resolutionBias(fabs(eta)) *(ptreco-ptgen))/ptreco;   
  }
  return ptreco*cor;
}


void Process( TString  *fname , bool doTheFirst, TString myMethodList = "" ) 
{

  TString _analysis = ANALYSIS;

  // 0 or 1
  int isFirst =doTheFirst;
  int isSecond = !doTheFirst;
  if ( (isFirst * isSecond) ==1) { std::cout << "error, you should either run on the first or the second daugthers...."<<std::endl;}

  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();

  std::cout << std::endl;
  std::cout << "==> Start TMVARegressionApplication" << std::endl;

  // --- Create the Reader object

  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  
  float hjetpt[2],  hjetptraw[2], hjetptrawold[2],  hjetptleadtrack[2], hjetphi[2],  hjetgenpt[2], hjeteta[2], hjetet[2], hjete[2],  hjetmt[2],  hjetchf[2], hjetnhf[2], hjetcef[2], hjetnef[2], hjetnconstituents[2], hjetnch[2], hjetvtxmass[2], hjetvtxpt[2],  hjetvtx3dl[2], hjetvtx3del[2], hjetnconstituents[2], hjetJECUnc[2], fathFilterJetspt[2],  fathFilterJetsptraw[2], fathFilterJetsptrawold[2],  fathFilterJetsptleadtrack[2], fathFilterJetsphi[2],  fathFilterJetsgenpt[2], fathFilterJetseta[2], fathFilterJetset[2], fathFilterJetse[2],  fathFilterJetsmt[2],  fathFilterJetschf[2], fathFilterJetsnhf[2], fathFilterJetscef[2], fathFilterJetsnef[2], fathFilterJetsnconstituents[2], fathFilterJetsnch[2], fathFilterJetsvtxmass[2], fathFilterJetsvtxpt[2],  fathFilterJetsvtx3dl[2], fathFilterJetsvtx3del[2];  
  float rho25, MET, METet, METdPhi;
  if (isFirst){
    if (_analysis == "Dijet") {
      reader->AddVariable("hJet_pt", &hjetpt[0]);
      reader->AddVariable("hJet_eta", &hjeteta[0]);
      reader->AddVariable("hJet_phi", &hjetphi[0]);
      reader->AddVariable("hJet_e", &hjete[0]);
      reader->AddVariable("hJet_ptRaw*((hJet_ptRaw+resolutionBias(fabs(hJet_eta))*(hJet_ptRaw-hJet_genPt))/hJet_ptRaw)", &hjetptraw[0]);
      //reader->AddVariable("newptRaw:=evalJERBias(hJet_ptRaw, hJet_genPt, hJet_eta)", &hjetptraw[0]);
      reader->AddVariable("hJet_Mt:=evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetmt[0]);
      reader->AddVariable("hJet_et:=evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetet[0]);
      reader->AddVariable("hJet_ptLeadTrack", &hjetptleadtrack[0]);
      reader->AddVariable("hJet_vtxPt", &hjetvtxpt[0]);
      reader->AddVariable("hJet_vtx3dL", &hjetvtx3dl[0]);
      reader->AddVariable("hJet_vtx3deL", &hjetvtx3del[0]);
      reader->AddVariable("hJet_vtxMass", &hjetvtxmass[0]);
      reader->AddVariable("hJet_chf", &hjetchf[0]);
      reader->AddVariable("hJet_nch", &hjetnch[0]);
      reader->AddVariable("hJet_nconstituents", &hjetnconstituents[0]);
      reader->AddVariable("hJet_JECUnc", &hjetJECUnc[0]);
      reader->AddVariable("rho25", &rho25);
      reader->AddVariable("MET.et", &METet);
      reader->AddVariable("METdPhi:=METdeltaPhi(MET.phi, hJet_phi[0], hJet_phi[1])", &METdPhi);
    }
    if (_analysis == "Subjet") {
      reader->AddVariable("fathFilterJets_pt", &fathFilterJetspt[0]);
      reader->AddVariable("fathFilterJets_eta", &fathFilterJetseta[0]);
      reader->AddVariable("fathFilterJets_phi", &fathFilterJetsphi[0]);
      reader->AddVariable("fathFilterJets_e", &fathFilterJetse[0]);
      reader->AddVariable("fathFilterJets_ptRaw*((fathFilterJets_ptRaw+resolutionBias(fabs(fathFilterJets_eta))*(fathFilterJets_ptRaw-fathFilterJets_genPt))/fathFilterJets_ptRaw)", &fathFilterJetsptraw[0]);
      reader->AddVariable("fathFilterJets_Mt:=evalMt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)", &fathFilterJetsmt[0]);
      reader->AddVariable("fathFilterJets_et:=evalEt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)", &fathFilterJetset[0]);
      reader->AddVariable("fathFilterJets_ptLeadTrack", &fathFilterJetsptleadtrack[0]);
      reader->AddVariable("fathFilterJets_vtxPt", &fathFilterJetsvtxpt[0]);
      reader->AddVariable("fathFilterJets_vtx3dL", &fathFilterJetsvtx3dl[0]);
      reader->AddVariable("fathFilterJets_vtx3deL", &fathFilterJetsvtx3del[0]);
      reader->AddVariable("fathFilterJets_vtxMass", &fathFilterJetsvtxmass[0]);
      reader->AddVariable("fathFilterJets_chf", &fathFilterJetschf[0]);
      reader->AddVariable("rho25", &rho25);
      reader->AddVariable("MET.et", &METet);
      reader->AddVariable("METdPhi:=METdeltaPhi(MET.phi, fathFilterJets_phi[0], fathFilterJets_phi[1])", &METdPhi);
    }
  }
  else{
    if (_analysis == "Dijet") {
      reader->AddVariable("hJet_pt", &hjetpt[1]);
      reader->AddVariable("hJet_eta", &hjeteta[1]);
      reader->AddVariable("hJet_phi", &hjetphi[1]);
      reader->AddVariable("hJet_e", &hjete[1]);
      reader->AddVariable("hJet_ptRaw*((hJet_ptRaw+resolutionBias(fabs(hJet_eta))*(hJet_ptRaw-hJet_genPt))/hJet_ptRaw)", &hjetptraw[1]);
      //reader->AddVariable("newptRaw:=evalJERBias(hJet_ptRaw, hJet_genPt, hJet_eta)", &hjetptraw[1]);
      reader->AddVariable("hJet_Mt:=evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetmt[1]);
      reader->AddVariable("hJet_et:=evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)", &hjetet[1]);
      reader->AddVariable("hJet_ptLeadTrack", &hjetptleadtrack[1]);
      reader->AddVariable("hJet_vtxPt", &hjetvtxpt[1]);
      reader->AddVariable("hJet_vtx3dL", &hjetvtx3dl[1]);
      reader->AddVariable("hJet_vtx3deL", &hjetvtx3del[1]);
      reader->AddVariable("hJet_vtxMass", &hjetvtxmass[1]);
      reader->AddVariable("hJet_chf", &hjetchf[1]);
      reader->AddVariable("hJet_nch", &hjetnch[1]);
      reader->AddVariable("hJet_nconstituents", &hjetnconstituents[1]);
      reader->AddVariable("hJet_JECUnc", &hjetJECUnc[1]);
      reader->AddVariable("rho25", &rho25);
      reader->AddVariable("MET.et", &METet);
      reader->AddVariable("METdPhi:=METdeltaPhi(MET.phi, hJet_phi[0], hJet_phi[1])", &METdPhi);
    }
    if (_analysis == "Subjet") {
      reader->AddVariable("fathFilterJets_pt", &fathFilterJetspt[1]);
      reader->AddVariable("fathFilterJets_eta", &fathFilterJetseta[1]);
      reader->AddVariable("fathFilterJets_phi", &fathFilterJetsphi[1]);
      reader->AddVariable("fathFilterJets_e", &fathFilterJetse[1]);
      reader->AddVariable("fathFilterJets_ptRaw*((fathFilterJets_ptRaw+resolutionBias(fabs(fathFilterJets_eta))*(fathFilterJets_ptRaw-fathFilterJets_genPt))/fathFilterJets_ptRaw)", &fathFilterJetsptraw[1]);
      reader->AddVariable("fathFilterJets_Mt:=evalMt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)", &fathFilterJetsmt[1]);
      reader->AddVariable("fathFilterJets_et:=evalEt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)", &fathFilterJetset[1]);
      reader->AddVariable("fathFilterJets_ptLeadTrack", &fathFilterJetsptleadtrack[1]);
      reader->AddVariable("fathFilterJets_vtxPt", &fathFilterJetsvtxpt[1]);
      reader->AddVariable("fathFilterJets_vtx3dL", &fathFilterJetsvtx3dl[1]);
      reader->AddVariable("fathFilterJets_vtx3deL", &fathFilterJetsvtx3del[1]);
      reader->AddVariable("fathFilterJets_vtxMass", &fathFilterJetsvtxmass[1]);
      reader->AddVariable("fathFilterJets_chf", &fathFilterJetschf[1]);
      reader->AddVariable("rho25", &rho25);
      reader->AddVariable("MET.et", &METet);
      reader->AddVariable("METdPhi:=METdeltaPhi(MET.phi, fathFilterJets_phi[0], fathFilterJets_phi[1])", &METdPhi);
    }
  }

  // --- Book the MVA methods

  TString dir    = "weights_Reg_8TeV";
  TString prefix = "TMVARegression_BDT";

  TString methodName = "BDT method";
  TString weightfile = dir + "_" + _analysis + "/" + prefix +  ".weights.xml";
  std::cout << "--- TMVARegressionApp        : Using weight file: " << weightfile << std::endl;
  reader->BookMVA( methodName, weightfile ); 


  TFile *input(0);
  input = TFile::Open( *fname ); // check if file in local directory exists 
  if (!input) {
    std::cout << "ERROR: could not open data file" << std::endl;
    exit(1);
  }
  std::cout << "--- TMVARegressionApp        : Using input file: " << input->GetName() << std::endl;

  // --- Event loop

  TTree*  theTree= (TTree*)input->Get("tree");
  std::cout << "--- Select signal sample" << std::endl;

  float hJet_genPtReg0 , hJet_genPtReg1, fathFilterJets_genPtReg0, fathFilterJets_genPtReg1; 
   
  if (_analysis == "Dijet") {
    if (isSecond) theTree->SetBranchAddress( "hJet_genPtReg0", &hJet_genPtReg0 );
    theTree->SetBranchAddress( "hJet_pt", &hjetpt );	  
    theTree->SetBranchAddress( "rho25", &rho25 );	  
    theTree->SetBranchAddress( "hJet_ptRaw", &hjetptrawold );	   
    //theTree->SetBranchAddress( "hJet_ptRaw*((hJet_ptRaw+resolutionBias(fabs(hJet_eta))*(hJet_ptRaw-hJet_genPt))/hJet_ptRaw)", &hjetptraw );
    theTree->SetBranchAddress( "hJet_ptLeadTrack", &hjetptleadtrack );	   
    theTree->SetBranchAddress( "hJet_phi", &hjetphi );	   
    theTree->SetBranchAddress( "hJet_genPt", &hjetgenpt );	   
    theTree->SetBranchAddress( "hJet_eta", &hjeteta );	   
    theTree->SetBranchAddress( "hJet_e", &hjete);		   
    theTree->SetBranchAddress( "hJet_chf", &hjetchf );	   
    theTree->SetBranchAddress( "hJet_nch", &hjetnch );	   
    theTree->SetBranchAddress( "hJet_nconstituents", &hjetnconstituents );	   
    theTree->SetBranchAddress( "hJet_nhf", &hjetnhf );	   
    theTree->SetBranchAddress( "hJet_cef", &hjetcef );	   
    theTree->SetBranchAddress( "hJet_chf", &hjetchf );	   
    theTree->SetBranchAddress( "hJet_JECUnc", &hjetJECUnc );	   
    theTree->SetBranchAddress( "hJet_vtxMass", &hjetvtxmass );
    theTree->SetBranchAddress( "hJet_vtxPt", &hjetvtxpt );
    theTree->SetBranchAddress( "hJet_vtx3dL", &hjetvtx3dl );  
    theTree->SetBranchAddress( "hJet_vtx3deL", &hjetvtx3del );  
    theTree->SetBranchAddress( "MET", &MET);
  }
  if (_analysis == "Subjet") {
    if (isSecond) theTree->SetBranchAddress( "fathFilterJets_genPtReg0", &fathFilterJets_genPtReg0 );
    theTree->SetBranchAddress( "fathFilterJets_pt", &fathFilterJetspt );	  
    theTree->SetBranchAddress( "rho25", &rho25 );	  
    theTree->SetBranchAddress( "fathFilterJets_ptRaw", &fathFilterJetsptrawold );	   
    //theTree->SetBranchAddress( "fathFilterJets_ptRaw*((fathFilterJets_ptRaw+resolutionBias(fabs(fathFilterJets_eta))*(fathFilterJets_ptRaw-fathFilterJets_genPt))/fathFilterJets_ptRaw)", &fathFilterJetsptraw );
    theTree->SetBranchAddress( "fathFilterJets_ptLeadTrack", &fathFilterJetsptleadtrack );	   
    theTree->SetBranchAddress( "fathFilterJets_phi", &fathFilterJetsphi );	   
    theTree->SetBranchAddress( "fathFilterJets_genPt", &fathFilterJetsgenpt );	   
    theTree->SetBranchAddress( "fathFilterJets_eta", &fathFilterJetseta );	   
    theTree->SetBranchAddress( "fathFilterJets_e", &fathFilterJetse);		   
    theTree->SetBranchAddress( "fathFilterJets_chf", &fathFilterJetschf );	   
    theTree->SetBranchAddress( "fathFilterJets_cef", &fathFilterJetscef );	   
    theTree->SetBranchAddress( "fathFilterJets_chf", &fathFilterJetschf );	   
    theTree->SetBranchAddress( "fathFilterJets_vtxMass", &fathFilterJetsvtxmass );
    theTree->SetBranchAddress( "fathFilterJets_vtxPt", &fathFilterJetsvtxpt );
    theTree->SetBranchAddress( "fathFilterJets_vtx3dL", &fathFilterJetsvtx3dl );  
    theTree->SetBranchAddress( "fathFilterJets_vtx3deL", &fathFilterJetsvtx3del );  
    theTree->SetBranchAddress( "MET", &MET);
  }
  
  TH1F *  count = (TH1F*)input->Get("Count");
  TH1F *  countWithPU = (TH1F*)input->Get("CountWithPU");
  TH1F *  countWithPU2011B = (TH1F*)input->Get("CountWithPU2011B");
  TH3F *  input3DPU = (TH3F*)input->Get("Input3DPU");

  std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
  TStopwatch sw;
  sw.Start();

  float newHiggsPt, newHiggsMass, newfatHiggsPt, newfatHiggsMass;
  float hJet_ptCorr[2], fathFilterJets_ptCorr[2];

  TFile *nfile2 = new TFile(fname->Data() ,"recreate");
  TTree *ntree = theTree->CopyTree("");

  std::cout << "Entries " << count->GetEntries() << std::endl;
  count->Clone()->Write();
  input3DPU->Clone()->Write();
  countWithPU->Clone()->Write();
  countWithPU2011B->Clone()->Write();

  delete count;
  delete countWithPU;
  delete countWithPU2011B;
  delete input3DPU;
    
  TBranch * b_hJet_genPtReg0 ;
  TBranch * b_hJet_genPtReg1 ;
  TBranch * b_hJet_ptCorr ;
  TBranch * b_newHiggsPt;
  TBranch * b_newHiggsMass;
  TBranch * b_fathFilterJets_genPtReg0 ;
  TBranch * b_fathFilterJets_genPtReg1 ;
  TBranch * b_fathFilterJets_ptCorr ;
  TBranch * b_newfatHiggsPt;
  TBranch * b_newfatHiggsMass;
  if (_analysis == "Dijet") {
    if (isFirst) {
      b_hJet_genPtReg0 = ntree->Branch("hJet_genPtReg0"   ,  &hJet_genPtReg0            ,  "hJet_genPtReg0/F"       );        
      b_hJet_ptCorr = ntree->Branch("hJet_ptCorr"   ,  hJet_ptCorr            ,  "hJet_ptCorr[2]/F"       );        
    } 
    else{
      b_hJet_genPtReg1 = ntree->Branch("hJet_genPtReg1"   ,  &hJet_genPtReg1            ,  "hJet_genPtReg1/F"       );          
      b_newHiggsPt = ntree->Branch("newHiggsPt"   ,  &newHiggsPt            ,  "newHiggsPt/F"       );          
      b_newHiggsMass = ntree->Branch("newHiggsMass"   ,  &newHiggsMass            ,  "newHiggsMass/F"       );          
    }
  //    TBranch * b_hjetpt = ntree->Branch("hjetpt" ,  &hjetpt            ,  "hjetpt[2]/F");        
  // TBranch * b_hjeteta = ntree->Branch("hjeteta" ,  &hjeteta            ,  "hjeteta[2]/F");  
  }  
  if (_analysis == "Subjet") {
    if (isFirst) {
      b_fathFilterJets_genPtReg0 = ntree->Branch("fathFilterJets_genPtReg0"   ,  &fathFilterJets_genPtReg0            ,  "fathFilterJets_genPtReg0/F"       );        
      b_fathFilterJets_ptCorr = ntree->Branch("fathFilterJets_ptCorr"   ,  fathFilterJets_ptCorr            ,  "fathFilterJets_ptCorr[2]/F"       );        
    } 
    else{
      b_fathFilterJets_genPtReg1 = ntree->Branch("fathFilterJets_genPtReg1"   ,  &fathFilterJets_genPtReg1            ,  "fathFilterJets_genPtReg1/F"       );          
      b_newfatHiggsPt = ntree->Branch("newfatHiggsPt"   ,  &newfatHiggsPt            ,  "newfatHiggsPt/F"       );          
      b_newfatHiggsMass = ntree->Branch("newfatHiggsMass"   ,  &newfatHiggsMass            ,  "newfatHiggsMass/F"       );          
    } 
  } 


  for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
  //for (Long64_t ievt=0; ievt<10000;ievt++) {

    //theTree->GetEntry(ievt);
    ntree->GetEntry(ievt);

    if (_analysis == "Dijet") {
      hjetet[0] = evalEt(hjetpt[0], hjeteta[0], hjetphi[0], hjete[0]);
      hjetet[1] = evalEt(hjetpt[1], hjeteta[1], hjetphi[1], hjete[1]);
      hjetmt[0] = evalMt(hjetpt[0], hjeteta[0], hjetphi[0], hjete[0]);
      hjetmt[1] = evalMt(hjetpt[1], hjeteta[1], hjetphi[1], hjete[1]);

      if ( hjetgenpt[0] < 0 ) {
	hjetptraw[0] = hjetptrawold[0];
      }
      else{
	hjetptraw[0] = hjetptrawold[0]*((hjetptrawold[0]+resolutionBias(fabs(hjeteta[0]))*(hjetptrawold[0]-hjetgenpt[0]))/hjetptrawold[0]);
	//hjetptraw[0] = evalJERBias( hjetptrawold[0], hjetgenpt[0], hjeteta[0]);
      
      }
      if ( hjetgenpt[0] < 0 ) {
	hjetptraw[1] = hjetptrawold[1];
      }
      else{
	hjetptraw[1] = hjetptrawold[1]*((hjetptrawold[1]+resolutionBias(fabs(hjeteta[1]))*(hjetptrawold[1]-hjetgenpt[1]))/hjetptrawold[1]);
	//hjetptraw[1] = evalJERBias( hjetptrawold[1], hjetgenpt[1], hjeteta[1]);
      }

      if (isFirst) {
	hJet_genPtReg0 = (reader->EvaluateRegression( "BDT method" )[0])     ; 
	b_hJet_genPtReg0->Fill();

      }else {
	hJet_genPtReg1 = (reader->EvaluateRegression( "BDT method" )[0])     ; 
	b_hJet_genPtReg1->Fill();
	hJet_ptCorr[1] = hJet_genPtReg1;

	newHiggsPt =     ptHiggs(hJet_genPtReg0, hjeteta[0], hjetphi[0], hjete[0], hJet_genPtReg1, hjeteta[1], hjetphi[1], hjete[1], hjetpt[0], hjetpt[1]);
	newHiggsMass = massHiggs(hJet_genPtReg0, hjeteta[0], hjetphi[0], hjete[0], hJet_genPtReg1, hjeteta[1], hjetphi[1], hjete[1], hjetpt[0], hjetpt[1]);
      
	b_newHiggsPt->Fill();
	b_newHiggsMass->Fill();
      }
    }

    if (_analysis == "Subjet") {
      fathFilterJetset[0] = evalEt(fathFilterJetspt[0], fathFilterJetseta[0], fathFilterJetsphi[0], fathFilterJetse[0]);
      fathFilterJetset[1] = evalEt(fathFilterJetspt[1], fathFilterJetseta[1], fathFilterJetsphi[1], fathFilterJetse[1]);
      fathFilterJetsmt[0] = evalMt(fathFilterJetspt[0], fathFilterJetseta[0], fathFilterJetsphi[0], fathFilterJetse[0]);
      fathFilterJetsmt[1] = evalMt(fathFilterJetspt[1], fathFilterJetseta[1], fathFilterJetsphi[1], fathFilterJetse[1]);
      
      if ( fathFilterJetsgenpt[0] < 0 ) {
	fathFilterJetsptraw[0] = fathFilterJetsptrawold[0];
      }
      else{
	fathFilterJetsptraw[0] = fathFilterJetsptrawold[0]*((fathFilterJetsptrawold[0]+resolutionBias(fabs(fathFilterJetseta[0]))*(fathFilterJetsptrawold[0]-fathFilterJetsgenpt[0]))/fathFilterJetsptrawold[0]);
	//fathFilterJetsptraw[0] = evalJERBias( fathFilterJetsptrawold[0], fathFilterJetsgenpt[0], fathFilterJetseta[0]);
      
      }
      if ( fathFilterJetsgenpt[0] < 0 ) {
	fathFilterJetsptraw[1] = fathFilterJetsptrawold[1];
      }
      else{
	fathFilterJetsptraw[1] = fathFilterJetsptrawold[1]*((fathFilterJetsptrawold[1]+resolutionBias(fabs(fathFilterJetseta[1]))*(fathFilterJetsptrawold[1]-fathFilterJetsgenpt[1]))/fathFilterJetsptrawold[1]);
	//fathFilterJetsptraw[1] = evalJERBias( fathFilterJetsptrawold[1], fathFilterJetsgenpt[1], fathFilterJetseta[1]);
      }

      if (isFirst) {
	fathFilterJets_genPtReg0 = (reader->EvaluateRegression( "BDT method" )[0])     ; 
	b_fathFilterJets_genPtReg0->Fill();
	fathFilterJets_ptCorr[0] = fathFilterJets_genPtReg0;

      }else {
	fathFilterJets_genPtReg1 = (reader->EvaluateRegression( "BDT method" )[0])     ; 
	b_fathFilterJets_genPtReg1->Fill();
	fathFilterJets_ptCorr[1] = fathFilterJets_genPtReg1;

	newfatHiggsPt =     ptHiggs(fathFilterJets_genPtReg0, fathFilterJetseta[0], fathFilterJetsphi[0], fathFilterJetse[0], fathFilterJets_genPtReg1, fathFilterJetseta[1], fathFilterJetsphi[1], fathFilterJetse[1], fathFilterJetspt[0], fathFilterJetspt[1]);
	newfatHiggsMass = massHiggs(fathFilterJets_genPtReg0, fathFilterJetseta[0], fathFilterJetsphi[0], fathFilterJetse[0], fathFilterJets_genPtReg1, fathFilterJetseta[1], fathFilterJetsphi[1], fathFilterJetse[1], fathFilterJetspt[0], fathFilterJetspt[1]);
      
	b_newfatHiggsPt->Fill();
	b_newfatHiggsMass->Fill();
      }
    }

    
    if (ievt%10000 == 0) {
      std::cout << "--- ... Processing event: " << ievt << std::endl;
      if (_analysis == "Dijet") {
	if (isFirst) 
	  {std::cout << "ievt: " << ievt <<" regres output, hjet pt 0 , hJet gen pt 0  " << hJet_genPtReg0 << " , " << hjetpt[0] << " , " << hjetgenpt[0] << "\n";}
	if (isSecond) {
	  std::cout << "ievt: " << ievt <<" regres output, hjet pt 0 , hJet gen pt 0  " << hJet_genPtReg0 << " , " << hjetpt[0] << " , " << hjetgenpt[0] << "\n";
	  std::cout << "regres output, hjet pt 1 , hJet gen pt 1  " << hJet_genPtReg1 << " , " << hjetpt[1] << " , " << hjetgenpt[1] << "\n";
	  std::cout << "new Higgs mass, pt  " << newHiggsMass << " , " << newHiggsPt  << "\n";
	}
      }
      if (_analysis == "Subjet") {
	if (isFirst) 
	  {std::cout << "ievt: " << ievt <<" regres output, fathFilterJet pt 0 , fathFilterJets gen pt 0  " << fathFilterJets_genPtReg0 << " , " << fathFilterJetspt[0] << " , " << fathFilterJetsgenpt[0] << "\n";}
	if (isSecond) {
	  std::cout << "ievt: " << ievt <<" regres output, fathFilterJets pt 0 , fathFilterJets gen pt 0  " << fathFilterJets_genPtReg0 << " , " << fathFilterJetspt[0] << " , " << fathFilterJetsgenpt[0] << "\n";
	  std::cout << "regres output, fathFilterJets pt 1 , fathFilterJets gen pt 1  " << fathFilterJets_genPtReg1 << " , " << fathFilterJetspt[1] << " , " << fathFilterJetsgenpt[1] << "\n";
	  std::cout << "new fat Higgs mass, pt  " << newfatHiggsMass << " , " << newfatHiggsPt  << "\n";
	}
      }
    }
      
  }

  sw.Stop();
  std::cout << "--- End of event loop: "; sw.Print();
 
  ntree->Write("", TObject::kOverwrite);
  nfile2->Write("", TObject::kOverwrite);

  // --- Write histograms
  delete reader;
    
  std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}

void TMVARegressionApplicationProcessor(){

  TString s(gROOT->GetVersion()); 
  if (!(s.Contains("5.34"))){
    std::cout << "INCORRECT ROOT VERSION!" << std::endl;
    std::cout << "run " << "source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh" << std:: endl;
    std::cout << "source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.01/x86_64-slc5-gcc43-opt/root/bin/thisroot.csh" << std::endl;
    return ;
  } 

  
  const int size =1;

  string titles[size] = {
    
    "/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/ZH_125_summer12_33b.root"
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_ZZ_TuneZ2star_8TeV_pythia6_tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_WW_TuneZ2star_8TeV_pythia6_tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_WZ_TuneZ2star_8TeV_pythia6_tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_TTJets_Merged.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_T_s-channel_TuneZ2star_8TeV-powheg-tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root",
    //"/uscmst1b_scratch/lpc1/lpctrig/ingabu/Substructre/Regression/DiJetPt_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root"

  }

  for(int k=0;k<size;k++) {
    std::cout<< " file " << titles[k] << std::endl;
    TString *fTT  = new TString(titles[k]);
    Process(fTT, 1, "BDT");
    Process(fTT, 0 , "BDT");
  }


}
