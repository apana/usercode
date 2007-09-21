#ifndef MyL1_L1Analyzer_L1JetPlots_h
#define MyL1_L1Analyzer_L1JetPlots_h
/* \class L1JetPlots
 *
 * \author Leonard Apanasevich
 *
 * \version 1
 *
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MyL1/L1Analyzer/interface/L1JetPlots.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>
#include <map>

using namespace edm;
using namespace reco;
using namespace std;

class L1JetPlots : public edm::EDAnalyzer {
public:
  L1JetPlots( const edm::ParameterSet & );

private:
  void beginJob( const edm::EventSetup & );
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();

  void L1Analysis( const CaloJetCollection&,
		   const GenJetCollection&,
		   const l1extra::L1JetParticleCollection&);

  template <typename T> void mtchL1(const Double_t&, const Double_t&, const Double_t&, 
				    const T& jets, const TString&);

  std::string CaloJetAlgorithm, GenJetAlgorithm, recmet_,genmet_, histogram;
  edm::InputTag l1CollectionsTag_;
  int errCnt;
  const int errMax(){return 100;}

  void fillHist(const TString& histName, const Double_t& value, const Double_t& wt=1.0);
  void fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt=1.0);

  bool doCaloJets,doGenJets,doCaloMET,doGenMET;
  bool doL1Jets;

  TH1F ptCal, etaCal, phiCal;
  TH1F ptGen, etaGen, phiGen;

  TH1F MetPt, genMetPt;

  TFile* m_file;

  // use the map function to access the rest of the histograms
  std::map<TString, TH1*> m_HistNames;
  std::map<TString, TH1*>::iterator hid;

  std::map<TString, TH2*> m_HistNames2D;
  std::map<TString, TH2*>::iterator hid2D;

};

#endif
