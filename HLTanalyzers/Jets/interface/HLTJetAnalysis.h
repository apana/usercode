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

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

typedef std::vector<std::string> MyStrings;

/** \class HLTJetAnalysis
  *  
  * $Date: 2007/04/10 18:38:16 $
  * $Revision: 1.3 $
  * \author L. Apanasevich - UIC
  */
class HLTJetAnalysis {
public:
  HLTJetAnalysis(); 

  void setup(const edm::ParameterSet& pSet);
  void fillHist(const TString& histName, const Double_t& value, const Double_t& wt=1.0);
  void fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt=1.0);
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
	       const l1extra::L1JetParticleCollection& l1jets,
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

  void L1Analysis( const CaloJetCollection& caloJets,
		   const GenJetCollection& genJets,
		   const l1extra::L1JetParticleCollection& l1Jets);

  template <typename T> void mtchL1(const Double_t&, const Double_t&, const Double_t&, 
				    const T& jets, const TString&);

  void bookMCParticles();
  void fillMCParticles(const HepMC::GenEvent mctruth);

  void bookHLTHistograms();
  void bookL1Histograms();

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
  bool doGenJets, doCaloJets, doL1Jets;

  const float etaBarrel() {return 1.4;}
  const float etaEndcap() {return 2.5;}

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
  bool evtTriggered;

};

#endif
