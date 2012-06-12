#ifndef Examples_CompJets_h
#define Examples_CompJets_h
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

class TFile;

class CompJets : public edm::EDAnalyzer {
public:
  CompJets( const edm::ParameterSet & );

private:
  void beginJob();
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();

  void getHLTResults(const edm::TriggerResults&, const edm::TriggerNames&);
  void plotJetPtandEta(const reco::PFJetCollection &, const std::string&);

  void fillHist(const TString& histName, const Double_t& value, const Double_t& wt=1.0);
  void fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt=1.0);
  void fill3DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& z,const Double_t& wt=1.0);

  void bookHistograms();

  TH1F* Book1dHist(const char*, const char*, Int_t, Double_t, Double_t, bool);
  TH2F* Book2dHist(const char*, const char*, Int_t, Double_t, Double_t, Int_t , Double_t , Double_t ,bool );
  TH3F* Book3dHist(const char*, const char*, Int_t, Double_t, Double_t, Int_t , Double_t ,
		   Double_t ,Int_t , Double_t , Double_t ,bool );

  edm::Service<TFileService> fs;

  edm::InputTag HLTriggerResults,PFJetCollection1_,PFJetCollection2_, GenJetCollection_;
  std::string MyTrigger;
  bool Debug_;

  TH1F *h_TriggerResults;

  // store hlt information in a map
  std::vector<bool> hlttrigs;
  std::map <std::string,bool> hltTriggerMap;
  std::map<std::string,bool>::iterator trig_iter;

  edm::TriggerNames triggerNames_;  // TriggerNames class

  std::map<TString, TH1*> m_HistNames;
  std::map<TString, TH1*>::iterator hid;

  std::map<TString, TH2*> m_HistNames2D;
  std::map<TString, TH2*>::iterator hid2D;

  std::map<TString, TH3*> m_HistNames3D;
  std::map<TString, TH3*>::iterator hid3D;

  int errCnt;
  bool HLTinit_;

  const int errMax(){return 20;}
  const float drMatch() {return 0.25;}

};

#endif
