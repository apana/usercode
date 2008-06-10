// HLTPlotsExample.cc
// Description:  Example of simple EDAnalyzer for jets.
// Author: Robert M. Harris
// Date:  28 - August - 2006
//
#include "MyHLT/HLTexamples/interface/HLTPlotsExample.h"

using namespace edm;
using namespace reco;
using namespace std;

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
HLTPlotsExample::HLTPlotsExample( const ParameterSet & cfg ) :
  HLTriggerResults( cfg.getParameter<InputTag>( "HLTriggerResults" ) ),
  CaloJetAlgorithm( cfg.getParameter<InputTag>( "CaloJetAlgorithm" ) ),
  GenJetAlgorithm( cfg.getParameter<InputTag>( "GenJetAlgorithm" ) ),
  HistoFile( cfg.getParameter<string>( "Histogram" ) ),
  MyTrigger( cfg.getParameter<string>( "MyTrigger" ) ),
  HLTinit_(false)
  {
}

void HLTPlotsExample::beginJob( const EventSetup & ) {

  // Open the histogram file and book some associated histograms
  m_file=new TFile(HistoFile.c_str(),"RECREATE");

  h_evtCounter    =  TH1F( "evtCounter",  "Event Counter", 10, -0.5, 9.5 );
  h_ptCal         =  TH1F( "ptCalAll",  "p_{T} of CaloJets", 100, 0, 500 );
  h_ptCalLeading  =  TH1F( "ptCal",  "p_{T} of leading CaloJets", 100, 0, 500 );
  h_ptCalTrig     =  TH1F( "ptCalTrig",  "p_{T} of leading CaloJets -- Triggered", 100, 0, 500 );

  h_etaCalLeading = TH1F( "etaCal", "#eta of leading CaloJets", 50, -3, 3 );
  h_phiCalLeading = TH1F( "phiCal", "#phi of leading CaloJets", 50, -M_PI, M_PI );
}

void HLTPlotsExample::analyze( const Event& evt, const EventSetup& es ) {

  bool gotHLT=true;
  bool myTrig=false;

  Handle<TriggerResults> hltresults,hltresultsDummy;
  evt.getByLabel(HLTriggerResults,hltresults);
  if (! hltresults.isValid() ) { cout << "  -- No HLTRESULTS"; gotHLT=false;}

  if (gotHLT) {
    getHLTResults(*hltresults);
    trig_iter=hltTriggerMap.find(MyTrigger);
    if (trig_iter==hltTriggerMap.end()){
      cout << "Could not find trigger path with name: " << MyTrigger << endl;
    }else{
      myTrig=trig_iter->second;
    }
  }

  h_evtCounter.Fill(0.); // count number of events processed
  if (myTrig) h_evtCounter.Fill(1.); // count number of events that fired my Trigger

  //Get the CaloJet collection
  Handle<CaloJetCollection> caloJets,caloJetsDummy;
  evt.getByLabel( CaloJetAlgorithm, caloJets );
  if (caloJets.isValid()) { 
    //Loop over the CaloJets and fill some histograms
    int jetInd = 0;
    for( CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++ cal ) {
      // std::cout << "CALO JET #" << jetInd << std::endl << cal->print() << std::endl;
      h_ptCal.Fill( cal->pt() );
      if (jetInd == 0){
	h_ptCalLeading.Fill( cal->pt() );
	h_etaCalLeading.Fill( cal->eta() );
	h_phiCalLeading.Fill( cal->phi() );
      
	if (myTrig) h_ptCalTrig.Fill( cal->pt() );
	jetInd++;
      }
    }
  }else{
    cout << "  -- No CaloJets" << endl;
  }

}

void HLTPlotsExample::getHLTResults( const edm::TriggerResults& hltresults) {


  int ntrigs=hltresults.size();

  if (! HLTinit_){
    HLTinit_=true;
    triggerNames_.init(hltresults);
    
    cout << "Number of HLT Paths: " << ntrigs << endl;

    // book histogram and label axis with trigger names
    h_TriggerResults = TH1F( "TriggerResults", "HLT Results", ntrigs, 0, ntrigs );

    for (int itrig = 0; itrig != ntrigs; ++itrig){
      string trigName = triggerNames_.triggerName(itrig);
      h_TriggerResults.GetXaxis()->SetBinLabel(itrig+1,trigName.c_str());
    }
  }

  
  for (int itrig = 0; itrig != ntrigs; ++itrig){
    string trigName = triggerNames_.triggerName(itrig);
     bool accept=hltresults.accept(itrig);

     if (accept) h_TriggerResults.Fill(float(itrig));

     // fill the trigger map
     typedef std::map<string,bool>::value_type valType;
     trig_iter=hltTriggerMap.find(trigName);
     if (trig_iter==hltTriggerMap.end())
       hltTriggerMap.insert(valType(trigName,accept));
     else
       trig_iter->second=accept;
  }
}

void HLTPlotsExample::endJob() {
  if (m_file !=0)
    {
    //Write the histograms to the file.
    m_file->cd();
    h_ptCal.Write();
    h_ptCalLeading.Write();
    h_ptCalTrig.Write();

    h_etaCalLeading.Write();
    h_phiCalLeading.Write();
    h_TriggerResults.Write();
    h_evtCounter.Write();
    delete m_file;
    m_file=0;
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTPlotsExample);
