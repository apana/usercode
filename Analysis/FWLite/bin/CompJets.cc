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

#include "DataFormats/Math/interface/deltaR.h"
#include "Analysis/FWLite/interface/myRootIO.h"

typedef std::map<std::string,bool> trigmap_t;
typedef std::map<std::string,bool>::iterator trigmapiter_t;

trigmap_t getHLTResults(const edm::TriggerResults&, const edm::TriggerNames&);
bool  checkJetID(const std::vector<reco::PFJet>::const_iterator);
void  BookHistograms(TFileDirectory);

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
  bool Debug_( ana.getParameter<bool>( "Debug" ) );
  int TriggerPS_( ana.getParameter<int>( "TriggerPS" ) );
  uint BeginLumi_( ana.getParameter<uint>("BeginLumi") );
  uint EndLumi_  ( ana.getParameter<uint>("EndLumi") );

  // Prepare output file and book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputHandler_.file().c_str());
  TFileDirectory dir = fs.mkdir("compJets");
  BookHistograms(dir);


  // loop the events
  int ievt=0;  
  int maxEvents_( inputHandler_.maxEvents() );
  cout << "\nNumber of input files: " << inputHandler_.files().size() << "\n";
  cout << "Number of events to process: " << maxEvents_ << "\n";
  cout << "Output histograms written to: " << outputHandler_.file() << "\n" << endl;

  for(unsigned int iFile=0; iFile<inputHandler_.files().size(); ++iFile){
    // open input file (can be located on castor)
    cout << iFile << ": " << inputHandler_.files()[iFile] << endl;
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
	if (Trigger_ == "Any")
	  {
	    myTrig =true;
	  }
	else
	  {
	    if (trig_iter==hltTriggers.end()){
	      cout << "Could not find trigger path with name: " << Trigger_ << endl;
	      exit(0);
	    }else{
	      myTrig=trig_iter->second;
	    }
	  }
	if (myTrig){

	  fillHist("LumiBlocks",ev.luminosityBlock(),1.); // keep track of the lumiBlocks analysed
	  edm::Handle<std::vector<PFJet> > PFjets1,PFjets2;
	  event.getByLabel(InputJets1_, PFjets1);
	  event.getByLabel(InputJets2_, PFjets2);

	  if ( (! PFjets1.isValid()) or (! PFjets2.isValid()) ){
	    // std::cout << " " << std::endl; 
	    continue;
	  }

	
	  // loop over jet collection and fill histograms

	  for(std::vector<PFJet>::const_iterator jet1=PFjets1->begin(); jet1!=PFjets1->end(); ++jet1){

	    if ( jet1->pt () > 5. && fabs( jet1->eta() )<5.1) {

	      // check the jetID
	      if (checkJetID(jet1)) {
		double jpt = jet1->pt ();
		double jeta = jet1->eta ();
		fillHist("jetPt",jpt,1.);
		if (jeta<-4.){
		  if (jeta<-4.5)
		    {
		      fillHist("jetPt_e1",jpt,1.);
		    }
		  else
		    {
		      fillHist("jetPt_e2",jpt,1.);
		    }
		}
		double drmin=99.;
		std::vector<PFJet>::const_iterator miter;
		for(std::vector<PFJet>::const_iterator jet2=PFjets2->begin(); jet2!=PFjets2->end(); ++jet2){
		  double dr=deltaR(jet1->eta(),jet1->phi(),jet2->eta(),jet2->phi());
		  if (dr<drmin) {
		    drmin=dr;
		    miter=jet2;
		  }
		}
		fillHist("drMin",drmin,1.);
		double rat=(miter->pt())/jpt;

		if (Debug_ and (jeta<-4.5 and jpt>500.)){
		  cout << "jpt1:jpt2:eta: " << jpt << "\t" << miter->pt() << "\t" << jeta <<endl;
		}
		if (jpt >30 && jpt<=40){
		  fillHist("jetEta_pt30_40",jeta);
		  if (drmin<0.1){
		    if (fabs(jeta)<1.1){
		      fillHist("jetResp_pt30_40_eta1",rat);
		    }else if (fabs(jeta)<2.5){
		      fillHist("jetResp_pt30_40_eta2",rat);
		    }else if (fabs(jeta)<4.5){
		      fillHist("jetResp_pt30_40_eta3",rat);
		    }
		  }
		}else if (jpt >40 && jpt<=50){
		  if (drmin<0.1){
		    if (fabs(jeta)<1.1){
		      fillHist("jetResp_pt40_50_eta1",rat);
		    }else if (fabs(jeta)<2.5){
		      fillHist("jetResp_pt40_50_eta2",rat);
		    }else if (fabs(jeta)<4.5){
		      fillHist("jetResp_pt40_50_eta3",rat);
		    }
		  }
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
  cout << "\nTotal number of events processed: " << endl;
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

bool checkJetID(const std::vector<reco::PFJet>::const_iterator jet){

  // loose id for PFJets

  bool jid=false;

  double chf   = jet->chargedHadronEnergyFraction();
  double nhf   = (jet->neutralHadronEnergy() + jet->HFHadronEnergy())/jet->energy();
  double phf   = jet->photonEnergyFraction();
  double elf   = jet->electronEnergyFraction();
  double chm   = jet->chargedHadronMultiplicity();
  //int nhm   = jet->neutralHadronMultiplicity();
  //int phm   = jet->photonMultiplicity();
  //int elm   = jet->electronMultiplicity();
  int npr   = jet->chargedMultiplicity() + jet->neutralMultiplicity();

  jid  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(jet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(jet->eta())>2.4)) ;

  return jid;

}

void  BookHistograms(TFileDirectory dir){

  TString hname;
  TString htitle;

  hname="LumiBlocks"; htitle="Events per Luminosity Block";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle,   2000,   0.,  2000.,false);


  //hists vs pT
  int npt=1000.;
  double ptmin=0.,ptmax=5000.;

  hname="jetPt"; htitle="jet p_{T}";
  m_HistNames[hname]=Book1dHist(dir,hname,htitle, npt,   ptmin, ptmax); 

  hname="jetPt_e1"; htitle="jet p_{T} -4.5< #eta < -4.";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, npt,   ptmin, ptmax); 

  hname="jetPt_e2"; htitle="jet p_{T} -5.0< #eta < -4.5";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, npt,   ptmin, ptmax); 

  // DeltaR
  hname="drMin"; htitle="Min dr between jets";
  m_HistNames[hname]=Book1dHist(dir, hname, htitle, 100, 0., 2.); 

  //hists vs eta
  int neta=20;
  double etamin=-5.,etamax=5.;

  hname="jetEta_pt30_40"; htitle="jet #eta 30 < p_{T} <  40";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);


  //jet Resp vs eta

  int njr=20;
  double jrmin=0.,jrmax=2.;

  hname="jetResp_pt30_40_eta1"; htitle="jetResp |#eta|<1.1  30 < p_{T} <40";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, njr, jrmin, jrmax);

  hname="jetResp_pt30_40_eta2"; htitle="jetResp 1.1<|#eta|<2.5  30 < p_{T} <40";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, njr, jrmin, jrmax);

  hname="jetResp_pt30_40_eta3"; htitle="jetResp |#eta|>2.5  30 < p_{T} <40";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, njr, jrmin, jrmax);


  hname="jetResp_pt40_50_eta1"; htitle="jetResp |#eta|<1.1  40 < p_{T} <50";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, njr, jrmin, jrmax);

  hname="jetResp_pt40_50_eta2"; htitle="jetResp 1.1<|#eta|<2.5  40 < p_{T} <50";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, njr, jrmin, jrmax);

  hname="jetResp_pt40_50_eta3"; htitle="jetResp |#eta|>2.5  40 < p_{T} <50";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, njr, jrmin, jrmax);

}
