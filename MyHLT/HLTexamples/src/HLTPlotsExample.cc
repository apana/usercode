// HLTPlotsExample.cc
// Description:  Example of simple EDAnalyzer for HLT.
// Author: L. Apanasevich
// Date:  28 - August - 2006
//
#include "MyHLT/HLTexamples/interface/HLTPlotsExample.h"

using namespace edm;
using namespace reco;
using namespace std;

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
HLTPlotsExample::HLTPlotsExample( const ParameterSet & cfg ) :
  HLTriggerResults1_( cfg.getParameter<InputTag>( "HLTriggerResults1" ) ),
  HLTriggerResults2_( cfg.getParameter<InputTag>( "HLTriggerResults2" ) ),
  CaloJetAlgorithm( cfg.getParameter<InputTag>( "CaloJetAlgorithm" ) ),
  muonCollection( cfg.getParameter<InputTag>( "MuonCollection" ) ),
  MyTrigger1( cfg.getParameter<string>( "MyTrigger1" ) ),
  MyTrigger2( cfg.getParameter<string>( "MyTrigger2" ) ),
  HLTinit_(false)
  {

  errCnt=0;
  oRun=0;
}

void HLTPlotsExample::beginJob() {

  h_lumi          =  fs->make<TH1F>( "lumi",  "Lumi Sections", 1000, 0.,1000. );
  h_lumi->Sumw2();
  h_evtCounter    =  fs->make<TH1F>( "evtCounter",  "Event Counter", 10, -0.5, 9.5 );
  h_ptCal         =  fs->make<TH1F>( "ptCal",  "CaloJet p_{T}", 50, 0, 50 );
  h_ptCalLeading  =  fs->make<TH1F>( "ptCalL",  "p_{T} of leading CaloJets", 50, 0, 50 );
  h_ptCalTrig     =  fs->make<TH1F>( "ptCalTrig",  "p_{T} of leading CaloJets -- Triggered", 50, 0, 50 );

  h_etaCalLeading = fs->make<TH1F>( "etaCal", "#eta of leading CaloJets", 50, -3, 3 );
  h_phiCalLeading = fs->make<TH1F>( "phiCal", "#phi of leading CaloJets", 50, -M_PI, M_PI );

  h_ptMuon        =  fs->make<TH1F>( "ptMuon",  "Muon p_{T}", 100, 0, 20 );
  h_ptMuonLeading =  fs->make<TH1F>( "ptMuonL", "p_{T} of Leading Muon", 100, 0, 20 );
  h_ptMuonTrig    =  fs->make<TH1F>( "ptMuonTrig", "p_{T} of Leading Muon -- Triggered", 100, 0, 20 );
}

void HLTPlotsExample::analyze( const Event& evt, const EventSetup& es ) {

  bool gotHLT=true;
  bool myTrig1=false, myTrig2=false;
  string errMsg("");


  Handle<TriggerResults> hltresults1,hltresults2;
  evt.getByLabel(HLTriggerResults1_,hltresults1);
  if (! hltresults1.isValid() ) { 
    gotHLT=false; errCnt+=1; errMsg=errMsg + "  -- No HLTriggerResults with process name: " + HLTriggerResults1_.process();
  }

  evt.getByLabel(HLTriggerResults2_,hltresults2);
  if (! hltresults2.isValid() ) { 
    gotHLT=false; errCnt+=1; errMsg=errMsg + "  -- No HLTriggerResults with process name: " + HLTriggerResults2_.process();
  }


  int iLumi = evt.luminosityBlock();
  int iRun = evt.id().run();
  int iEvent = evt.id().event();
  h_lumi->Fill(float(iLumi));

  if (iRun != oRun){
    cout << "Beginning to Process new run:" << iRun << endl;
    oRun=iRun;
    gotHLT=false;
  }
  if (gotHLT) {
    const TriggerNames & triggerNames_ = evt.triggerNames(*hltresults1);
    getHLTResults(*hltresults1, triggerNames_);
    trig_iter=hltTriggerMap.find(MyTrigger1);
    if (trig_iter==hltTriggerMap.end()){
      cout << "Could not find trigger path with name: " << MyTrigger1 << endl;
    }else{
      myTrig1=trig_iter->second;
    }
  }

  h_evtCounter->Fill(0.); // count number of events processed
  if (! myTrig1) return; 
  h_evtCounter->Fill(1.); // count number of events that fired my Trigger

  if (gotHLT) {
    const TriggerNames & triggerNames_ = evt.triggerNames(*hltresults2);
    getHLTResults(*hltresults2, triggerNames_);
    trig_iter=hltTriggerMap.find(MyTrigger2);
    if (trig_iter==hltTriggerMap.end()){
      cout << "Could not find trigger path with name: " << MyTrigger2 << endl;
    }else{
      myTrig2=trig_iter->second;
    }
  }

  //Get the CaloJet collection
  Handle<CaloJetCollection> caloJets,caloJetsDummy;
  evt.getByLabel( CaloJetAlgorithm, caloJets );
  if (caloJets.isValid()) { 
    //Loop over the CaloJets and fill some histograms
    int jetInd = 0;
    for( CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); ++ cal ) {
      // std::cout << "CALO JET #" << jetInd << std::endl << cal->print() << std::endl;
      h_ptCal->Fill( cal->pt() );
      if (jetInd == 0){
	h_ptCalLeading->Fill( cal->pt() );
	h_etaCalLeading->Fill( cal->eta() );
	h_phiCalLeading->Fill( cal->phi() );
      
	if (myTrig2) h_ptCalTrig->Fill( cal->pt() );
	jetInd++;
      }
    }
  }else{
    errMsg=errMsg + "  -- No CaloJets";
  }

  // Get the Muon Collection
  Handle<MuonCollection> muon,muonDummy;
  evt.getByLabel(muonCollection,muon);
  if (muon.isValid()) { 
    double ptmax=-999.;
    typedef MuonCollection::const_iterator muiter,mumax;
    for (muiter i=muon->begin(); i!=muon->end(); i++) {
      double mupt=i->pt();
      if (mupt > ptmax) ptmax=mupt;
      h_ptMuon->Fill( mupt );
    }
    if (ptmax > 0.){
      h_ptMuonLeading->Fill( ptmax );
      if (myTrig2) h_ptMuonTrig->Fill( ptmax );
    }
  }else{
    errMsg=errMsg + "  -- No Reco Muons";
  }


  if ((errMsg != "") && (errCnt < errMax())){
    errCnt=errCnt+1;
    errMsg=errMsg + ".";
    std::cout << "%MyHLT-Warning" << errMsg << std::endl;
    if (errCnt == errMax()){
      errMsg="%MyHLT-Warning -- Maximum error count reached -- No more messages will be printed.\n";
      std::cout << errMsg << std::endl;    
    }
  }

}

void HLTPlotsExample::getHLTResults( const edm::TriggerResults& hltresults,
				     const edm::TriggerNames& triggerNames_) {


  int ntrigs=hltresults.size();

  if (! HLTinit_){
    HLTinit_=true;
    // triggerNames_.init(hltresults);
    
    cout << "\nNumber of HLT Paths: " << ntrigs << "\n\n";

    // book histogram and label axis with trigger names
    h_TriggerResults = fs->make<TH1F>( "TriggerResults", "HLT Results", ntrigs, 0, ntrigs );

    for (int itrig = 0; itrig != ntrigs; ++itrig){
      string trigName = triggerNames_.triggerName(itrig);
      h_TriggerResults->GetXaxis()->SetBinLabel(itrig+1,trigName.c_str());
    }
  }

  
  for (int itrig = 0; itrig != ntrigs; ++itrig){
    string trigName = triggerNames_.triggerName(itrig);
     bool accept=hltresults.accept(itrig);

     if (accept) h_TriggerResults->Fill(float(itrig));

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

  cout << "\n%%%%%%%%%%%%%%%%  Job Summary %%%%%%%%%%%%%%%%%\n\n" ;

  if (h_evtCounter){
    double ntot=h_evtCounter->GetBinContent(1);
    cout  << "\tNumber of events processed: " << int(ntot) << "\n\n";
  }

  if (h_TriggerResults){

    int nbins=h_TriggerResults->GetNbinsX(); 

    cout  << "\tHLT Algorithm \t\t # of Accepts" << "\n";
    cout  << "\t------------ \t\t ------------" << "\n";
    for (int ibin=0; ibin<nbins; ++ibin){
      float cont=h_TriggerResults->GetBinContent(ibin+1);
      string trigName = string (h_TriggerResults->GetXaxis()->GetBinLabel(ibin+1));
      
      if (!trigName.empty()){
	cout  << "\t" << trigName << ":\t" << cont << endl;
      }
    }
    cout << "\n" << endl;
  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HLTPlotsExample);
