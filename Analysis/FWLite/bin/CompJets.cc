#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CommonTools/Utils/interface/PtComparator.h"

typedef std::map<std::string,bool> trigmap_t;
typedef std::map<std::string,bool>::iterator trigmapiter_t;

trigmap_t getHLTResults(const edm::TriggerResults&, const edm::TriggerNames&);


using std::cout;
using std::endl;
using std::string;
using reco::PFJet;

int main(int argc, char* argv[]) 
{

  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // parse arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  if( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ){
    std::cout << " ERROR: ParametersSet 'process' is missing in your configuration file" << std::endl; exit(0);
  }
  // get the python configuration
  const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
  fwlite::InputSource inputHandler_(process); fwlite::OutputFiles outputHandler_(process);


  // now get each parameter
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("compJets");
  edm::InputTag TriggerResults_( ana.getParameter<edm::InputTag>("TriggerResults") );
  edm::InputTag InputJets1_( ana.getParameter<edm::InputTag>("InputJets1") );
  edm::InputTag InputJets2_( ana.getParameter<edm::InputTag>("InputJets2") );
  string Trigger_( ana.getParameter<string>( "Trigger" ) );
  int TriggerPS_( ana.getParameter<int>( "TriggerPS" ) );
  uint BeginLumi_( ana.getParameter<uint>("BeginLumi") );
  uint EndLumi_  ( ana.getParameter<uint>("EndLumi") );

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputHandler_.file().c_str());
  TFileDirectory dir = fs.mkdir("compJets");
  TH1F* hLumiBlocks  = dir.make<TH1F>("LumiBlocks"  , "Events per Luminosity Block"  ,   2000,   0.,  2000.);
  TH1F* jetPt_  = dir.make<TH1F>("jetPt"  , "pt"  ,   100,   0., 500.); 
  jetPt_->Sumw2();

  TH1F* jetEta_pt30_40_ = dir.make<TH1F>("jetEta"  , "jet #eta 30 < p_{T} <40"  ,   20,   -5., 5.); 
  jetEta_pt30_40_->Sumw2();

  // loop the events
  int ievt=0;  
  int maxEvents_( inputHandler_.maxEvents() );
  for(unsigned int iFile=0; iFile<inputHandler_.files().size(); ++iFile){
    // open input file (can be located on castor)
    cout << inputHandler_.files()[iFile] << endl;
    TFile* inFile = TFile::Open(inputHandler_.files()[iFile].c_str());
    if( inFile ){
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	if(inputHandler_.reportAfter()!=0 ? (ievt>0 && ievt%inputHandler_.reportAfter()==0) : false) {
	  cout << "Processing " << ievt<< "th event: "
	       << "run " << ev.id().run() 
	       << ", lumi " << ev.luminosityBlock() 
	       << ", evt " << ev.id().event() << endl;
	}
	
	if (ev.luminosityBlock() < BeginLumi_ || ev.luminosityBlock() > EndLumi_) continue;

	edm::Handle<edm::TriggerResults> TriggerResults;
	event.getByLabel(TriggerResults_,TriggerResults);
	if (! TriggerResults.isValid()){
	  std::cout << " ERROR: Could not read TriggerResults object -- Exiting" << std::endl; 
	  exit(0);
	}

	const edm::TriggerNames triggerNames = event.triggerNames(*TriggerResults);
	trigmap_t hltTriggers=getHLTResults(*TriggerResults,triggerNames);


	// filter on Trigger Bits from Original Trigger Collections: RefTriggerResults
	trigmapiter_t trig_iter=hltTriggers.find(Trigger_);

	bool myTrig=false;
	if (trig_iter==hltTriggers.end()){
	  cout << "Could not find trigger path with name: " << Trigger_ << endl;
	  exit(0);
	}else{
	  myTrig=trig_iter->second;
	}
	if (myTrig){

	  hLumiBlocks->Fill(ev.luminosityBlock(),1); // keep track of the lumiBlocks analysed

	  edm::Handle<std::vector<PFJet> > PFjets;
	  event.getByLabel(InputJets1_, PFjets);
	
	  // Sort the jet collections
	  reco::PFJetCollection sortedPFjets;
	  sortedPFjets = *PFjets;
	  GreaterByPt<PFJet>   pTComparator_PF;
	  std::sort(sortedPFjets.begin(),sortedPFjets.end(),pTComparator_PF);

	  // loop jet collection and fill histograms

	  for(std::vector<PFJet>::const_iterator jet=sortedPFjets.begin(); jet!=sortedPFjets.end(); ++jet){

	    if ( jet->pt () > 5. && fabs( jet->eta() )<5.1) {

	      // check the jetID
	      double chf   = jet->chargedHadronEnergyFraction();
	      double nhf   = (jet->neutralHadronEnergy() + jet->HFHadronEnergy())/jet->energy();
	      double phf   = jet->photonEnergyFraction();
	      double elf   = jet->electronEnergyFraction();
	      double chm   = jet->chargedHadronMultiplicity();
	      //int nhm   = jet->neutralHadronMultiplicity();
	      //int phm   = jet->photonMultiplicity();
	      //int elm   = jet->electronMultiplicity();
	      int npr   = jet->chargedMultiplicity() + jet->neutralMultiplicity();
	      bool looseID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(jet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(jet->eta())>2.4)) ;
	      if (looseID) {
		double jpt = jet->pt ();
		double jeta = jet->eta ();
		jetPt_ ->Fill( jpt,1. );
		if (jpt >30 && jpt<=40){
		  jetEta_pt30_40_ ->Fill( jeta,1. );
		}
	      }
	    }
	  }



	} // end of myTrig block
      }  
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }

  return 0;
}

//std::map <std::string,bool> getHLTResults(const edm::TriggerResults& hltresults,
//					  const edm::TriggerNames& triggerNames_){

trigmap_t getHLTResults(const edm::TriggerResults& hltresults,
					  const edm::TriggerNames& triggerNames_){

  trigmap_t hltTriggerMap;
  std::map<std::string,bool>::iterator trig_iter;

  int ntrigs=hltresults.size();
  // std::cout << "\nNumber of HLT Paths: " << ntrigs << "\n\n";

  for (int itrig = 0; itrig != ntrigs; ++itrig){
    string trigName = triggerNames_.triggerName(itrig);
    bool accept=hltresults.accept(itrig);
    
    //if (accept) h_TriggerResults->Fill(float(itrig));

    // fill the trigger map
    typedef std::map<string,bool>::value_type valType;
    trig_iter=hltTriggerMap.find(trigName);
    if (trig_iter==hltTriggerMap.end())
      hltTriggerMap.insert(valType(trigName,accept));
    else
      trig_iter->second=accept;
  }

  return hltTriggerMap;
}
