#include <iostream>

#include "HLTanalyzers/Jets/interface/HLTJetAnalysis.h"

#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

class HLTJetAnalyzer : public edm::EDAnalyzer {
public:
  explicit HLTJetAnalyzer(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
  virtual void endJob();
private:
  // variables persistent across events should be declared here.
  //
  ///Default analyses
  HLTJetAnalysis jet_analysis_;
  std::string recjets_,genjets_,recmet_,genmet_,calotowers_,hltobj_,hltresults_;
  edm::InputTag l1CollectionsTag_;
  int errCnt1,errCnt2;
  const int errMax(){return 100;}
};
