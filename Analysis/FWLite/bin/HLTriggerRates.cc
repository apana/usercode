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

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "Analysis/FWLite/interface/myRootIO.h"

typedef std::map<std::string,bool> trigmap_t;
typedef std::map<std::string,bool>::iterator trigmapiter_t;

trigmap_t getHLTResults(const edm::TriggerResults&, const edm::TriggerNames&);
void BookJetHistograms(TFileDirectory);

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
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("triggerAnalyzer");
  edm::InputTag RefTriggerResults_( ana.getParameter<edm::InputTag>("RefTriggerResults") );
  edm::InputTag TriggerResults_( ana.getParameter<edm::InputTag>("TriggerResults") );
  edm::InputTag InputJets_( ana.getParameter<edm::InputTag>("InputJets") );
  string RefTrigger_( ana.getParameter<string>( "RefTrigger" ) );
  string Trigger_( ana.getParameter<string>( "Trigger" ) );
  int TriggerPS_( ana.getParameter<int>( "TriggerPS" ) );
  uint BeginLumi_( ana.getParameter<uint>("BeginLumi") );
  uint EndLumi_  ( ana.getParameter<uint>("EndLumi") );
  double RefLumi_( ana.getParameter<double>("RefLumi"   ));
  double TargetLumi_( ana.getParameter<double>("TargetLumi"   ));

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputHandler_.file().c_str());
  TFileDirectory dir = fs.mkdir("triggerRates");
  TH1F* hLumiBlocks  = dir.make<TH1F>("LumiBlocks"  , "Events per Luminosity Block"  ,   2000,   0.,  2000.);
  TH1F* hLumiBlocksRef  = dir.make<TH1F>("LumiBlocksTrig"  , "Events per Luminosity Block -- RefTrigger passed"  ,   2000,   0.,  2000.);
  TH1F* hHLTRates=NULL;

  BookJetHistograms(dir);

  unsigned int hltFlg=0;
  
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

	hLumiBlocks->Fill(ev.luminosityBlock(),1); // keep track of the lumiBlocks analysed

	edm::Handle<edm::TriggerResults> RefTriggerResults, TriggerResults;
	event.getByLabel(RefTriggerResults_,RefTriggerResults);
	if (! RefTriggerResults.isValid()){
	  std::cout << " ERROR: Could not read RefTriggerResults object -- Exiting" << std::endl; 
	  exit(0);
	}
	event.getByLabel(TriggerResults_,TriggerResults);
	if (! TriggerResults.isValid()){
	  std::cout << " ERROR: Could not read TriggerResults object -- Exiting" << std::endl; 
	  exit(0);
	}

	const edm::TriggerNames reftriggerNames = event.triggerNames(*RefTriggerResults);
	trigmap_t hltRefTriggers=getHLTResults(*RefTriggerResults,reftriggerNames);

	const edm::TriggerNames triggerNames = event.triggerNames(*TriggerResults);
	trigmap_t hltTriggers=getHLTResults(*TriggerResults,triggerNames);


	// filter on Trigger Bits from Original Trigger Collections: RefTriggerResults
	trigmapiter_t trig_iter=hltRefTriggers.find(RefTrigger_);

	bool myTrig=false;
	if (trig_iter==hltRefTriggers.end()){
	  cout << "Could not find trigger path with name: " << RefTrigger_ << endl;
	  exit(0);
	}else{
	  myTrig=trig_iter->second;
	}
	if (myTrig){

	  hLumiBlocksRef->Fill(ev.luminosityBlock(),1); // keep track of the lumiBlocks analysed
	  int nTriggers=hltTriggers.size();
	  // bool accept[nTriggers];
	  if (hltFlg == 0){
	    hltFlg++;
	    
	    hHLTRates= dir.make<TH1F>("hltrates","HLT Rates",nTriggers,0.,nTriggers);
	    hHLTRates->Sumw2();

	    int itrig=0;
	    for ( trigmapiter_t iter = hltTriggers.begin(); iter != hltTriggers.end(); ++iter){
	      // cout << iter->first << '\t' << iter->second << '\n';
	      string trigName = iter->first;
	      //cout << trigName << endl;
	      hHLTRates->GetXaxis()->SetBinLabel(itrig+1,trigName.c_str());
	      itrig++;
	    }

	    // std::cout << "Number of triggers: " << nTriggers << std::endl;
	    //for (int itrig = 0; itrig != nTriggers; ++itrig) {
	    //  std::string trigName = triggerNames.triggerName(itrig);
	    //  hHLTRates->GetXaxis()->SetBinLabel(itrig+1,trigName.c_str());
	    //  std::cout << trigName << std::endl;
	    //}
	    
	  }

	  for ( trigmapiter_t iter = hltTriggers.begin(); iter != hltTriggers.end(); ++iter){
	    string trigName = iter->first;
	    bool accept = iter->second;
	    if(accept) hHLTRates->Fill(trigName.c_str(),1);
	  }

	  trigmapiter_t mytrig_iter=hltTriggers.find(Trigger_);
	  bool yesTrig=false;
	  if (mytrig_iter==hltTriggers.end()){
	    cout << "Could not find trigger path with name: " << Trigger_ << endl;
	    exit(0);
	  }else{
	    yesTrig=mytrig_iter->second;
	  }

	  edm::Handle<std::vector<PFJet> > PFjets;
	  event.getByLabel(InputJets_, PFjets);
	  if (! PFjets.isValid()){
	    // std::cout << " " << std::endl; 
	    continue;
	  }

	  int nj=PFjets->size();
	  fillHist("Njets",nj,1.);

	  //std::cout << " Number of jets: " << nj << std::endl; 
	  // Sort the jet collections
	  reco::PFJetCollection sortedPFjets;
	  sortedPFjets = *PFjets;
	  GreaterByPt<PFJet>   pTComparator_PF;
	  std::sort(sortedPFjets.begin(),sortedPFjets.end(),pTComparator_PF);

	  nj=sortedPFjets.size();
	  // std::cout << " Number of sorted jets: " << nj << std::endl; 
	  int ijet=0;
	  for(std::vector<PFJet>::const_iterator jet=sortedPFjets.begin(); jet!=sortedPFjets.end() && ijet < 1; ++jet){
	    ijet++;
	    if (ijet==1){
	      double jpt = jet->pt ();
	      double jeta = jet->eta ();
	      fillHist("jetPt",jpt,1.);
	      if (jpt >40 ){
		fillHist("jetEta40",jeta,1.);
		if (jpt >60 ) fillHist("jetEta60",jeta,1.);
		if (jpt >80 ) fillHist("jetEta80",jeta,1.);
		if (jpt >100) fillHist("jetEta100",jeta,1.);
		if (yesTrig) {
		  fillHist("jetEtaTrig40",jeta,1.);
		  if (jpt >60 ) fillHist("jetEtaTrig60",jeta,1.);
		  if (jpt >80 ) fillHist("jetEtaTrig80",jeta,1.);
		  if (jpt >100) fillHist("jetEtaTrig100",jeta,1.);
		  if (jeta<-4.){
		    if (jeta>-4.5){
		      fillHist("jetPt_e1",jpt,1.);
		    }else{
		      fillHist("jetPt_e2",jpt,1.);
		    }
		  }
		}
		//if (yesTrig and jeta<-4.5)
		//  cout << "XXX: " << jeta << "\t" << jpt << endl;
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


  int nLumiBins=hLumiBlocksRef->GetNbinsX(); 
  int nLumiBlocks=0;
  for (int ibin=0; ibin<nLumiBins; ++ibin){
    int cont=hLumiBlocksRef->GetBinContent(ibin+1);
    if (cont>0) nLumiBlocks++;
  }
  float fact=nLumiBlocks*23.3;
  float LumiFact=TargetLumi_/RefLumi_;
  hHLTRates->Scale((LumiFact*TriggerPS_)/fact);

  sclHist("jetPt",(LumiFact*TriggerPS_)/fact);
  sclHist("jetPt_e1",(LumiFact*TriggerPS_)/fact);
  sclHist("jetPt_e2",(LumiFact*TriggerPS_)/fact);

  sclHist("jetEta40",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEta60",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEta80",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEta100",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEtaTrig40",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEtaTrig60",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEtaTrig80",(LumiFact*TriggerPS_)/fact);
  sclHist("jetEtaTrig100",(LumiFact*TriggerPS_)/fact);

  double nEvtsTrg=hLumiBlocksRef->GetEntries();
  cout << "\tNumber of events passed by reference trigger: " << int(nEvtsTrg) << "\n";
  cout << "\tNumber of Lumi sections processed: " << nLumiBlocks << "\n";
  cout << "\tReference trigger rate (wrt original run): " << nEvtsTrg/(fact) << " Hz \n\n";



  
  if (hHLTRates){

    int nbins=hHLTRates->GetNbinsX(); 

    cout  << "\tHLT Algorithm \t\t # of Accepts" << "\n";
    cout  << "\t------------ \t\t ------------" << "\n";
    for (int ibin=0; ibin<nbins; ++ibin){
      float cont=hHLTRates->GetBinContent(ibin+1);
      float err =hHLTRates->GetBinError(ibin+1);
      int ncounts=(cont*fact)/TriggerPS_;
      string trigName = string (hHLTRates->GetXaxis()->GetBinLabel(ibin+1));
      
      if (!trigName.empty()){
	cout  << "\t" << trigName << ":\t" << ncounts << "\tRate:\t" << std::setprecision( 3 ) << cont << " +- " << err << endl;
      }
    }
    cout << "\n" << endl;
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

void BookJetHistograms(TFileDirectory dir){

  TString hname;
  TString htitle;

  hname="Njets"; htitle="Number of PF jets in event";
  m_HistNames[hname]=Book1dHist(dir,hname,htitle, 100,   0., 100.); 


  //hists vs pT
  hname="jetPt"; htitle="jet p_{T}";
  m_HistNames[hname]=Book1dHist(dir,hname,htitle, 1000,   0., 5000.); 

  hname="jetPt_e1"; htitle="jet p_{T} -4.5< #eta < -4.";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, 1000, 0., 5000.); 

  hname="jetPt_e2"; htitle="jet p_{T} -5.0< #eta < -4.5";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, 1000, 0., 5000.);

  //hists vs eta
  int neta=20;
  double etamin=-5.,etamax=5.;

  hname="jetEta40"; htitle="jet #eta p_{T} > 40";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);
  hname="jetEta60"; htitle="jet #eta p_{T} > 60";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);
  hname="jetEta80"; htitle="jet #eta p_{T} > 80";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);
  hname="jetEta100"; htitle="jet #eta p_{T} > 100";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);

  hname="jetEtaTrig40"; htitle="jet #eta p_{T} > 40  -- Triggered ";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);
  hname="jetEtaTrig60"; htitle="jet #eta p_{T} > 60  -- Triggered";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);
  hname="jetEtaTrig80"; htitle="jet #eta p_{T} > 80  -- Triggered";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);
  hname="jetEtaTrig100"; htitle="jet #eta p_{T} > 100  -- Triggered";
  m_HistNames[hname]=Book1dHist(dir,hname, htitle, neta, etamin, etamax);

}
