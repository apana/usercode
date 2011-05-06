#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "Analysis/FWLite/interface/FWLiteExample.h"
using namespace std;

int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  // using reco::Muon;
  using reco::BasicJet;
  using reco::PFJet;

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

  // initialize command line parser
  optutl::CommandLineParser parser ("Analyze FWLite MYHistograms");

  // set defaults
  parser.integerValue ("maxEvents"  ) = 1000;
  parser.integerValue ("outputEvery") = 1000;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";

  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFilesList_ = parser.stringVector("inputFiles");

  if (readFilelist(inputFilesList_[0]) < 0 ) {
    cout << "There was a problem reading the filelists" << endl;
    return 1;
  }
  cout << "Number of files to process: "<< inputFiles_.size() << endl;

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("analyzeBasicPat");
  TH1F* jetPt_  = dir.make<TH1F>("jetPt"  , "pt"  ,   100,   0., 500.); jetPt_->Sumw2();
  TH1F* subjetPt_  = dir.make<TH1F>("subjetPt"  , "SubJet pt"  ,   100,   0., 500.);  subjetPt_->Sumw2();
  TH1F* jetLPt_  = dir.make<TH1F>("jetLeadingPt"  , "Leading Jet p_{T}"  ,   100,   0., 500.); jetLPt_->Sumw2();

  TH1F* jetEta_ = dir.make<TH1F>("jetEta" , "eta" ,   100,  -3.,   3.);  jetEta_->Sumw2();
  TH1F* jetPhi_ = dir.make<TH1F>("jetPhi" , "phi" ,   100,  -5.,   5.);  jetPhi_->Sumw2();

  TH1F* jetMass_   = dir.make<TH1F>("jetMass" , "Pruned Jet mass" ,   100,  0.,   150.);   jetMass_->Sumw2();
  TH1F* jetMassN2_   = dir.make<TH1F>("jetMassN2" , "Pruned Jet mass -- 2 good constituents" ,   100,  0.,   150.);   jetMassN2_->Sumw2();
  TH1F* MassDrop_ = dir.make<TH1F>("MassDrop" , "mass_s1 / massJ" ,   50,  0.,   1.); MassDrop_->Sumw2();

  TH1F* Asymmetry_ = dir.make<TH1F>("Asymmetry" , "Subjet Asymmetry" ,   50,  0.,   1.); Asymmetry_->Sumw2();
  TH1F* DeltaR_    = dir.make<TH1F>("DeltaR" , "Subjet Delta R" ,   50,  0.,   1.); DeltaR_->Sumw2();

  TH1F* NConst_ = dir.make<TH1F>("nConst" , "Number of Constituents" ,   5,  -0.5,   4.5); NConst_->Sumw2();

  bool isMC(true);

  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    TString inputFile=inputFiles_.at(iFile);
    cout << "Processing file "<< iFile+1 << ": "<< inputFile<< endl;

    TFile* inFile = TFile::Open(inputFile);
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
      edm::Handle<LumiSummary> summary;

      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) {
	  std::cout << "  processing event: " << ievt << std::endl;
	  std::cout <<" Run ID " << ev.getRun().id()<< std::endl;
	}


	double wt(1.);
	if (isMC){
	  edm::Handle< GenEventInfoProduct > GenInfoHandle;
	  std::string Gen("generator");
	  event.getByLabel( Gen, GenInfoHandle );
	  double qScale = GenInfoHandle->qScale();
	  double pthat = ( GenInfoHandle->hasBinningValues() ? 
			   (GenInfoHandle->binningValues())[0] : 0.0);
	  
	  //cout << " qScale = " << qScale << " pthat = " << pthat << endl;
	  wt = GenInfoHandle->weight();
	}
	//std::cout << " integrated event weight = " << wt << std::endl;


	std::string inputColl("ak5PFJets");
	// determine this is a good event using the ak5 jet collections
	edm::Handle<std::vector<PFJet> > PFjets;
	event.getByLabel(inputColl, PFjets);
	
	// Sort the jet collections
	reco::PFJetCollection sortedPFjets;
	sortedPFjets = *PFjets;
	GreaterByPt<PFJet>   pTComparator_PF;
	std::sort(sortedPFjets.begin(),sortedPFjets.end(),pTComparator_PF);

	// loop jet collection and fill histograms
	int ngoodjet=0;
	for(std::vector<PFJet>::const_iterator jet=sortedPFjets.begin(); jet!=sortedPFjets.end(); ++jet){

	  if ( jet->pt () > 200. && fabs( jet->eta() )<2.5) {

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
	    if (looseID) ngoodjet++;
	  }
	}
	bool goodEVT=ngoodjet>1;

	// done
	inputColl="caPrunedPFlow";

	// Handle to the jet collection
	edm::Handle<std::vector<BasicJet> > jets;
	event.getByLabel(inputColl, jets);
	
	// Sort the jet collections
	reco::BasicJetCollection myjets;
	myjets = *jets;
	GreaterByPt<BasicJet>   pTComparator_;
	std::sort(myjets.begin(),myjets.end(),pTComparator_);

	// loop jet collection and fill histograms

	int ijet=0;
	//for(std::vector<PFJet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
	for(std::vector<BasicJet>::const_iterator jet=myjets.begin(); jet!=myjets.end(); ++jet){
	  jetPt_ ->Fill( jet->pt (),wt );
	  jetEta_->Fill( jet->eta(),wt );
	  jetPhi_->Fill( jet->phi(),wt );	  
	  // require some jet cleaning
	  if (jet->pt () > 250. && fabs(jet->eta())<2.5 && goodEVT) {

	      ijet++;

	      // std::cout << "jet: "<<ijet << "  pT: " << jet->pt () << std::endl;

	      const std::vector<edm::Ptr<reco::Candidate> > iparticles = jet->getJetConstituents();
	      int nConstituents = iparticles.size();
	      NConst_->Fill(nConstituents,1.);

	      if (ijet==1 ){
		double mfat = jet->mass();
		jetLPt_ ->Fill( jet->pt (),wt );
		jetMass_ ->Fill( mfat,wt );
		
		double maxSubjetPt=0.;
		double maxSubjetMass=0;
		int nGoodConstituents(0);

		for (int i = 0; i <nConstituents ; i++){
		  double subpt  = iparticles[i]->pt();
		  double submass= iparticles[i]->mass();

		  if (subpt>5.) nGoodConstituents++;
		  subjetPt_ ->Fill( subpt,wt );
		  if (subpt > maxSubjetPt) {
		    maxSubjetPt=subpt;
		    maxSubjetMass=submass;
		  }
		}
		bool GoodConstituents= (nGoodConstituents>1);
		if (nConstituents ==2 && GoodConstituents){

		  // if (mu>0.95){
		  //   
		  //   cout << "Mu: " << mu << " m1: " << iparticles[0]->mass() << " m2: " << iparticles[1]->mass() << endl;
		  //   cout << "Mfat: " << mfat << " pt1: " << iparticles[0]->pt() << " pt2: " << iparticles[1]->pt() << endl;
		  // }

		  double m1=iparticles[0]->mass();
		  double m2=iparticles[1]->mass();
		  double pt1=iparticles[0]->pt();
		  double pt2=iparticles[1]->pt();
		  double mu=m1/mfat;
		  if (pt2>pt1) mu=m2/mfat;

		  double dR = reco::deltaR<double>( iparticles[0] ->eta(),
						    iparticles[0] ->phi(),
						    iparticles[1] ->eta(),
						    iparticles[1] ->phi()  );
		  double asy = std::min( pt1*pt1, pt2*pt2) * dR*dR / (mfat*mfat);
		  jetMassN2_ ->Fill( mfat,wt );
		  DeltaR_->Fill( dR,wt );
		  Asymmetry_->Fill( asy,wt );
		  MassDrop_->Fill( mu,wt );
		  
		}
	      }// end ijet==1

	  }// end jetPT and eta if
	}
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
