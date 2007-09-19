// JetPlotsExample.cc
// Description:  Example of simple EDAnalyzer for jets.
// Author: Robert M. Harris
// Date:  28 - August - 2006
//
#include "MyL1/L1JetPlots/JetPlotsExample.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>
using namespace edm;
using namespace reco;
using namespace std;

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
JetPlotsExample::JetPlotsExample( const ParameterSet & cfg ) :
  CaloJetAlgorithm( cfg.getParameter<string>( "CaloJetAlgorithm" ) ),
  GenJetAlgorithm( cfg.getParameter<string>( "GenJetAlgorithm" ) ),
  recmet_(cfg.getParameter< string > ("recmet")),
  genmet_(cfg.getParameter< string > ("genmet")),
  calotowers_(cfg.getParameter< string > ("calotowers")),
  histogram( cfg.getParameter<string>( "Histogram" ) )
  {
    globalThreshold=0.;
}

void JetPlotsExample::beginJob( const EventSetup & ) {

  // Open the histogram file and book some associated histograms
  m_file=new TFile(histogram.c_str(),"RECREATE");
  h_ptCal =  TH1F( "ptCal",  "p_{T} of leading CaloJets", 50, 0, 1000 );
  h_etaCal = TH1F( "etaCal", "#eta of leading CaloJets", 50, -3, 3 );
  h_phiCal = TH1F( "phiCal", "#phi of leading CaloJets", 50, -M_PI, M_PI );
  h_ptGen =  TH1F( "ptGen",  "p_{T} of leading GenJets", 50, 0, 1000 );
  h_etaGen = TH1F( "etaGen", "#eta of leading GenJets", 50, -3, 3 );
  h_phiGen = TH1F( "phiGen", "#phi of leading GenJets", 50, -M_PI, M_PI );

  int nbins=50;
  Double_t min=0.,max=150.;

  h_MetPt = TH1F( "h_MetPt", "Reconstructed MET ", nbins, min, max );
  h_genMetPt = TH1F( "h_genMetPt", "Generated MET ", nbins, min, max );

  h_MetPt_cT = TH1F( "h_MetPt_cT", "Reconstructed MET -- from CaloTowers", nbins, min, max );
  h_MetPt_cT_be = TH1F( "h_MetPt_cT_be", "Reconstructed MET -- from CaloTowers -- Barrel and Endcaps", nbins, min, max );
  h_MetPt_cT_f = TH1F( "h_MetPt_cT_f", "Reconstructed MET -- from CaloTowers -- Forward", nbins, min, max );

}

void JetPlotsExample::analyze( const Event& evt, const EventSetup& es ) {

  //Get the CaloJet collection
  Handle<CaloJetCollection> caloJets;
  evt.getByLabel( CaloJetAlgorithm, caloJets );

  //Loop over the two leading CaloJets and fill some histograms
  int jetInd = 0;
  for( CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end() && jetInd<2; ++ cal ) {
    // std::cout << "CALO JET #" << jetInd << std::endl << cal->print() << std::endl;
    h_ptCal.Fill( cal->pt() );
    h_etaCal.Fill( cal->eta() );
    h_phiCal.Fill( cal->phi() );
    jetInd++;
  }

  //Get the GenJet collection
  Handle<GenJetCollection> genJets;
  evt.getByLabel( GenJetAlgorithm, genJets );

  //Loop over the two leading GenJets and fill some histograms
  jetInd = 0;
  for( GenJetCollection::const_iterator gen = genJets->begin(); gen != genJets->end() && jetInd<2; ++ gen ) {
    // std::cout << "GEN JET #" << jetInd << std::endl << gen->print() << std::endl;
    h_ptGen.Fill( gen->pt() );
    h_etaGen.Fill( gen->eta() );
    h_phiGen.Fill( gen->phi() );
    jetInd++;
  }

  Handle<CaloMETCollection> recmet;
  evt.getByLabel (recmet_,recmet);
  Handle<GenMETCollection>  genmet;
  evt.getByLabel (genmet_,genmet);

  for( CaloMETCollection::const_iterator met = recmet->begin(); met != recmet->end() ; ++ met ) {
    h_MetPt.Fill(met->pt());
  }

  for( GenMETCollection::const_iterator met = genmet->begin(); met != genmet->end() ; ++ met ) {
    h_genMetPt.Fill(met->pt());
  }

  string errMsg("");
  try {
    Handle<CaloTowerCollection> caloTowers;
    evt.getByLabel (calotowers_,caloTowers);

    double sum_et = 0.0;
    double sum_ex = 0.0;
    double sum_ey = 0.0;

    double sum_et_be = 0.0;
    double sum_ex_be = 0.0;
    double sum_ey_be = 0.0;

    double sum_et_f = 0.0;
    double sum_ex_f = 0.0;
    double sum_ey_f = 0.0;

    double sum_ez = 0.0;
    for ( CaloTowerCollection::const_iterator tower=caloTowers->begin();
	  tower!=caloTowers->end(); tower++) {

      if( tower->et() > globalThreshold )
	{
	  double phi   = tower->phi();
	  //double theta = tower->theta();
	  double e     = tower->energy();
	  double et    = tower->et();
	  double eta   = tower->eta();

	  //sum_ez += e*cos(theta);
	  sum_et += et;
	  sum_ex += et*cos(phi);
	  sum_ey += et*sin(phi);

	  if (fabs(eta) < 2.5){
	    sum_et_be += et;
	    sum_ex_be += et*cos(phi);
	    sum_ey_be += et*sin(phi);
	  }else{
	    sum_et_f += et;
	    sum_ex_f += et*cos(phi);
	    sum_ey_f += et*sin(phi);
	  }
	}
    }
    double mex   = -sum_ex;
    double mey   = -sum_ey;
    //double mez   = -sum_ez;
    double met   = sqrt( sum_ex*sum_ex + sum_ey*sum_ey );
    double met_be   = sqrt( sum_ex_be*sum_ex_be + sum_ey_be*sum_ey_be );
    double met_f   = sqrt( sum_ex_f*sum_ex_f + sum_ey_be*sum_ey_f );

    double sumet = sum_et;
    double phi   = atan2( -sum_ey, -sum_ex );

    h_MetPt_cT.Fill(met);
    h_MetPt_cT_be.Fill(met_be);
    h_MetPt_cT_f.Fill(met_f);

  } catch (const cms::Exception& e){
    errMsg=errMsg + "  -- No CaloTowers\n"+e.what();
    std::cout << errMsg << std::endl;
  }

}

void JetPlotsExample::endJob() {

  //Write out the histogram file.
  m_file->Write();

}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetPlotsExample);
