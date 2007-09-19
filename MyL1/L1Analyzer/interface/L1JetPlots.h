#ifndef MyL1_L1Analyzer_L1JetPlots_h
#define MyL1_L1Analyzer_L1JetPlots_h
#include <TH1.h>
/* \class L1JetPlots
 *
 * \author Leonard Apanasevich
 *
 * \version 1
 *
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"

class L1JetPlots : public edm::EDAnalyzer {
public:
  L1JetPlots( const edm::ParameterSet & );

private:
  void beginJob( const edm::EventSetup & );
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();
  std::string CaloJetAlgorithm, GenJetAlgorithm, recmet_,genmet_, calotowers_, histogram;

  double globalThreshold;

  TH1F h_ptCal, h_etaCal, h_phiCal;
  TH1F h_ptGen, h_etaGen, h_phiGen;

  TH1F h_MetPt,h_genMetPt,h_MetPt_cT,h_MetPt_cT_be,h_MetPt_cT_f;


  TFile* m_file;
};

#endif
