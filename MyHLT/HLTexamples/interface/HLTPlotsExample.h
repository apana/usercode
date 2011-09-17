#ifndef Examples_HLTPlotsExample_h
#define Examples_HLTPlotsExample_h
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

class TFile;

class HLTPlotsExample : public edm::EDAnalyzer {
public:
  HLTPlotsExample( const edm::ParameterSet & );

private:
  void beginJob();
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();

  void getHLTResults(const edm::TriggerResults&, const edm::TriggerNames&);

  edm::Service<TFileService> fs;

  edm::InputTag HLTriggerResults1_,HLTriggerResults2_,CaloJetAlgorithm,muonCollection;
  std::string MyTrigger1,MyTrigger2;

  TH1F *h_ptCal, *h_ptCalLeading, *h_ptCalTrig, *h_etaCalLeading, *h_phiCalLeading;
  TH1F *h_ptMuon, *h_ptMuonLeading, *h_ptMuonTrig;

  TH1F *h_TriggerResults, *h_evtCounter, *h_lumi;

  // store hlt information in a map
  std::vector<bool> hlttrigs;
  std::map <std::string,bool> hltTriggerMap;
  std::map<std::string,bool>::iterator trig_iter;

  edm::TriggerNames triggerNames_;  // TriggerNames class

  int errCnt,oRun;
  bool HLTinit_;

  const int errMax(){return 20;}

};

#endif
