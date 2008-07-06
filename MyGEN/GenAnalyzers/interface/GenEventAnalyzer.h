#ifndef Examples_GenEventAnalyzer_h
#define Examples_GenEventAnalyzer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

class TFile;

class GenEventAnalyzer : public edm::EDAnalyzer {
public:
  GenEventAnalyzer( const edm::ParameterSet & );

private:
  void beginJob( const edm::EventSetup & );
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();

  void getGENINFO(const reco::CandidateView&, const double);

  edm::InputTag mctruth_, genEventScale_;
  //std::string HistoFile_;

  edm::Service<TFileService> fs;

  TH1F *h_evtCounter, *h_ptHat;
  TH1F *h_mxElePt,*h_mxMuPt,*h_mxTauPt;
  TH1F *h_dimuonMass;
  //TFile* m_file;

  // store hlt information in a map
//    std::vector<bool> hlttrigs;
//    std::map <std::string,bool> hltTriggerMap;
//    std::map<std::string,bool>::iterator trig_iter;
//    
};

#endif
