#ifndef HLTJETS_H
#define HLTJETS_H

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNamed.h"
#include <vector>
#include <map>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/ProductID.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetfwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetfwd.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef std::vector<std::string> MyStrings;

/** \class HLTJetAnalysis
  *  
  * $Date: $
  * $Revision: $
  * \author L. Apanasevich - UIC
  */
class HLTJetAnalysis {
public:
  HLTJetAnalysis(); 

  void setup(const edm::ParameterSet& pSet);
  void fillHist(const TString& histName, const Double_t& value, const Double_t& wt=1.0);
  void fillHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt=1.0);
  void extractPtHat(edm::EventSetup const& isetup);

  /** Analyze the Data */
  void analyze(const CaloJetCollection& rjets,
	       const GenJetCollection& gjets,
	       const CaloMETCollection& rmets,
	       const GenMETCollection& gmets,
	       const CaloTowerCollection& caloTowers,
	       const HepMC::GenEvent mctruth,
	       const HLTFilterObjectWithRefs& hltobj,
	       const edm::TriggerResults& hltresults,
	       const CaloGeometry& geom);

  void dummyAnalyze(
	       const CaloGeometry& geom);
  /** Finalization (close files, etc) */
  void done();


  void bookHistograms();
  void bookGeneralHistograms();

  void bookCaloTowerHists();
  void fillCaloTowerHists(const CaloTowerCollection& caloTowers);

  void bookMetHists(const TString& prefix);
  template <typename T> void fillMetHists(const T& mets, const TString& prefx);

  void bookJetHistograms(const TString& prefix);

  template <typename T> void fillJetHists(const T& jets, const TString& prefx);

  void bookMCParticles();
  void fillMCParticles(const HepMC::GenEvent mctruth);

  void bookHLTHistograms();
  void getHLTResults(const edm::TriggerResults& hltResults);
  void getHLTParticleInfo(const HLTFilterObjectWithRefs& hltobj);

private:

  // input variables
  string _HistName; // Name of histogram file
  string _HLTPath; // Name of trigger path to analyze
  bool _Monte,_Debug;
  double _EtaMin,_EtaMax;
  double _CKIN3, _CKIN4;


  int evtCounter;
  bool doGenJets, doCaloJets;

  const float etaBarrel() {return 1.4;}

  TFile* m_file; // pointer to Histogram file

  // histogram declarations
  TH1* m_Cntr; // Simple single histogram

  // use the map function to access the rest of the histograms
  std::map<TString, TH1*> m_HistNames;
  std::map<TString, TH1*>::iterator hid;

  std::map<TString, TH2*> m_HistNames2D;
  std::map<TString, TH2*>::iterator hid2D;

  //create maps linking histogram pointers to HCAL Channel hits and digis

  TString gjetpfx, rjetpfx,gmetpfx, rmetpfx,calopfx;

  vector<bool> hlttrigs;
  std::map <string,bool> hltTriggerMap;
  std::map<string,bool>::iterator trig_iter;
  bool hlttrig;
  bool hltInfoExists;

};

class PtGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.pt() > j.pt());
  }
};

#endif
