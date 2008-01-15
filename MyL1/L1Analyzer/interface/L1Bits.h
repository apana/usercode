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

#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>

using namespace edm;
using namespace reco;
using namespace std;

class L1Bits : public edm::EDAnalyzer {
public:
  L1Bits( const edm::ParameterSet & );

private:
  void beginJob( const edm::EventSetup & );
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();

  void L1Analysis(const l1extra::L1ParticleMapCollection& l1mapcoll);


  std::string histogram, particleMapSource_, text_output;

  int errCnt;
  const int errMax(){return 100;}

  void fillHist(const TString& histName, const Double_t& value, const Double_t& wt=1.0);
  void fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt=1.0);

  bool doCaloJets,doGenJets,doCaloMET,doGenMET;
  bool doL1Jets;

  TFile* m_file;

  TH1F* evtCounter;


  // use the map function to access the rest of the histograms
  std::map<TString, TH1*> m_HistNames;
  std::map<TString, TH1*>::iterator hid;

  std::map<TString, TH2*> m_HistNames2D;
  std::map<TString, TH2*>::iterator hid2D;

  std::map<string, int>m_bits;
  std::map<string, int>::iterator m_iter;

};

#endif
