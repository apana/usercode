// File: HLTJetAnalyzer.cc
// Description:  Example of Analysis driver originally from Jeremy Mans,
// Date:  13-October-2006

#include "HLTanalyzers/Jets/interface/HLTJetAnalyzer.h"

// Boiler-plate constructor definition of an analyzer module:
//
HLTJetAnalyzer::HLTJetAnalyzer(edm::ParameterSet const& conf) {

  // If your module takes parameters, here is where you would define
  // their names and types, and access them to initialize internal
  // variables. Example as follows:
  //
  std::cout << " Beginning HLTJetAnalyzer Analysis " << std::endl;

  recjets_    = conf.getParameter< std::string > ("recjets");
  genjets_    = conf.getParameter< std::string > ("genjets");
  recmet_     = conf.getParameter< std::string > ("recmet");
  genmet_     = conf.getParameter< std::string > ("genmet");
  l1CollectionsTag_     = conf.getParameter< edm::InputTag > ("l1collections");
  calotowers_ = conf.getParameter< std::string > ("calotowers");
  hltobj_    = conf.getParameter< std::string > ("hltobj");
  errCnt1=0;
  errCnt2=0;
  jet_analysis_.setup(conf);

}

// Boiler-plate "analyze" method declaration for an analyzer module.
//
void HLTJetAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


  // To get information from the event setup, you must request the "Record"
  // which contains it and then extract the object you need
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<IdealGeometryRecord>().get(geometry);

  // These declarations create handles to the types of records that you want
  // to retrieve from event "iEvent".
  //  
  edm::Handle<CaloJetCollection>  recjets;
  edm::Handle<GenJetCollection>  genjets;
  edm::Handle<CaloTowerCollection> caloTowers;
  edm::Handle<CaloMETCollection> recmet;
  edm::Handle<GenMETCollection> genmet;
  //edm::Handle<edm::HepMCProduct> mctruthHandle;
  edm::Handle<edm::HepMCProduct> genEventHandle;
  edm::Handle<CandidateCollection> genParticles;

  edm::Handle<HLTFilterObjectWithRefs> hltobj;
  //  edm::Handle<std::vector<edm::HLTPathStatus> > hltresults;
  edm::Handle<edm::TriggerResults> hltresults;


  edm::Handle<l1extra::L1JetParticleCollection> l1jets;


  
  // Extract Data objects (event fragments)
  // make sure to catch exceptions if they don't exist...
  string errMsg("");
  string trgMsg("");
  try {iEvent.getByLabel(recjets_,recjets);} catch (...) { errMsg=errMsg + "  -- No RecJets";}
  try {iEvent.getByLabel(recmet_,recmet);} catch (...) {errMsg=errMsg + "  -- No RecMET";}
  try {iEvent.getByLabel(calotowers_,caloTowers);} catch (...) {errMsg=errMsg + "  -- No CaloTowers";}

  try {
    iEvent.getByLabel(hltobj_,hltobj);
  } catch (...) {
    errMsg=errMsg + "  -- No HLTOBJ with name: " + hltobj_ ;
  }

  edm::InputTag L1JetTag(edm::InputTag(l1CollectionsTag_.label(),"Central"));
  try {
    iEvent.getByLabel(L1JetTag,l1jets);
  } catch (...) {
    errMsg=errMsg + "  -- No L1Jets with name: " + L1JetTag.label() ;
  }

  try {
    //iEvent.getByType(hltresults);
    //void Event::getManyByType( std::vector<edm::Handle<T> >& result) 
    std::vector<edm::Handle<edm::TriggerResults> > vTrigResults ;
    iEvent.getManyByType( vTrigResults ) ;
    //std::cout << "Number of TrigResults objects" << vTrigResults.size() << std::endl;
    if (vTrigResults.size()==1){
      hltresults=vTrigResults[0];
    }else{
      std::ostringstream nobj;
      nobj << vTrigResults.size() ;
      trgMsg="Multiple ( " + nobj.str() + " ) TriggerResults objects found\n" ;
      for (unsigned int i = 0; i < vTrigResults.size(); ++i) {
    	hltresults=vTrigResults[i];
    	//std::cout << hltresults.provenance()->moduleLabel() << std::endl;
    	//std::cout << hltresults.provenance()->moduleName() << std::endl;
    	//std::cout << hltresults.provenance()->processName() << std::endl;
    	//std::cout << "BranchName: " << hltresults.provenance()->branchName() << std::endl;
	trgMsg=trgMsg + "BranchName: " + hltresults.provenance()->branchName() + "\n";
    	//if ((HepMCEvt.provenance()->product).moduleLabel() == "VtxSmeared")
    	//  break;
      }     
      hltresults=vTrigResults[0];
      //std::cout << "Using: " << hltresults.provenance()->branchName() << std::endl;
      trgMsg=trgMsg + "Using: " + hltresults.provenance()->branchName() + "\n";
    }
  } catch (...) {
    errMsg=errMsg + "  -- No HLTRESULTS";
  }

  // MC objects
  HepMC::GenEvent genEvent;
  try {
    iEvent.getByLabel("source","",genEventHandle);
    genEvent = genEventHandle->getHepMCData();
    iEvent.getByLabel("genParticleCandidates",genParticles);
  } catch (...) {
    errMsg=errMsg + "  -- No MC truth";
  }

  try {
    iEvent.getByLabel (genjets_,genjets);
  } catch (...) {
    errMsg=errMsg + "  -- No GenJets";
  }

  try {
    iEvent.getByLabel (genmet_,genmet);
  } catch (...) {
    errMsg=errMsg + "  -- No GenMet";
  }

  if ((errMsg != "") && (errCnt1 < errMax())){
    errCnt1=errCnt1+1;
    errMsg=errMsg + ".";
    std::cout << "%HLTJetAnalyzer-Warning" << errMsg << std::endl;
    if (errCnt1 == errMax()){
      errMsg="%HLTJetAnalyzer-Warning -- Maximum error1 count reached -- No more messages will be printed.";
      std::cout << errMsg << std::endl;    
    }
  }

  if ((trgMsg != "") && (errCnt2 < errMax())){
    errCnt2=errCnt2+1;
    std::cout << trgMsg << std::endl;
    if (errCnt2 == errMax()){
      errMsg="%HLTJetAnalyzer-Warning -- Maximum error2 count reached -- No more messages will be printed.";
      std::cout << errMsg << std::endl;    
    }
  }

  // run the jet analysis, passing required event fragments
  //
  jet_analysis_.analyze(*recjets,*genjets,
			*recmet,*genmet,
			*caloTowers,genParticles,genEvent,
			*hltobj,*hltresults,
			*l1jets,
			*geometry);
}

// "endJob" is an inherited method that you may implement to do post-EOF processing
// and produce final output.
//
void HLTJetAnalyzer::endJob() {
  jet_analysis_.done();
}

