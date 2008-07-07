// GenEventAnalyzer.cc
// Description:  Simple EDAnalyzer for Generated MC events.
// Author: Leonard Apanasevich
// Date:  5 - July - 2009
//
#include "MyGEN/GenAnalyzers/interface/GenEventAnalyzer.h"

using namespace edm;
using namespace reco;
using namespace std;

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
GenEventAnalyzer::GenEventAnalyzer( const ParameterSet & cfg ) :
  mctruth_(cfg.getParameter<InputTag> ("mctruth") ),
  genEventScale_( cfg.getParameter<InputTag> ("genEventScale") )
  //HistoFile_( cfg.getParameter<string>( "Histogram" ) )
  {
}

void GenEventAnalyzer::beginJob( const EventSetup & ) {

  // Open the histogram file and book some associated histograms
  //m_file=new TFile(HistoFile_.c_str(),"RECREATE");

  h_evtCounter    =  fs->make<TH1F>( "evtCounter",  "Event Counter", 10, -0.5, 9.5 );
  h_ptHat         =  fs->make<TH1F>( "ptHat",   "p_{T} Hat of Hard Scatter", 1000, 0.0, 1000. );

  int npbins=200;
  double ptmin=0., ptmax=100.;
  h_mxElePt       =  fs->make<TH1F>( "ptEleMX", "p_{T} of Leading Electron", npbins, ptmin, ptmax );
  h_mxMuPt        =  fs->make<TH1F>( "ptMuMX" , "p_{T} of Leading Muon"    , npbins, ptmin, ptmax );
  h_mxTauPt       =  fs->make<TH1F>( "ptTauMX", "p_{T} of Leading Tau"     , npbins, ptmin, ptmax );

  int nmbins=200;
  double mmin=0., mmax=200.;
  h_dimuonMass    =  fs->make<TH1F>( "dimuonMass", "Di-Muon Invariant Mass", nmbins, mmin, mmax );
  h_dieleMass    =  fs->make<TH1F>( "dieleMass", "Di-Electron Invariant Mass", nmbins, mmin, mmax );
  h_ditauMass    =  fs->make<TH1F>( "ditauMass", "Di-Tau Invariant Mass", nmbins, mmin, mmax );
}

void GenEventAnalyzer::analyze( const Event& evt, const EventSetup& es ) {

  // cout << "Beginning GenEventAnalyzer::analyze" << endl;
  bool gotGEN=true;

  Handle<CandidateView> mctruth,mctruthDummy;
  Handle< double > genEventScale;

  evt.getByLabel(genEventScale_, genEventScale );
  evt.getByLabel(mctruth_,mctruth);

  double pthat=*genEventScale;
  if (! mctruth.isValid() ) { cout << "  -- No Gen Particles"; gotGEN=false;}

  h_evtCounter->Fill(0.); // count number of events processed
  if (gotGEN) {
    h_evtCounter->Fill(1.); 
    getGENINFO(*mctruth,pthat);
  }



}

void GenEventAnalyzer::getGENINFO(const CandidateView& mctruth ,const double pthat) {

  //cout << "Beginning GenEventAnalyzer::getGENINFO" << endl;
  h_ptHat->Fill(pthat); 


  if (&mctruth){

    double ptEleMX=-1.,ptMuMX=-1.,ptTauMX=-1.;
    vector<math::XYZTLorentzVector> muonP4,eleP4,tauP4;

    for (size_t i = 0; i < mctruth.size(); ++ i) {
      const Candidate & p = (mctruth)[i];
      int mcpid = p.pdgId();
      double mcpt = p.pt();
      int Status =  p.status();
      if (Status == 1){
	if (abs(mcpid) == 11 ) {
	  if (mcpt > ptEleMX) ptEleMX=mcpt;
	  eleP4.push_back(p.p4());
	}else if (abs(mcpid) == 13 ) {
	  if (mcpt > ptMuMX) ptMuMX=mcpt;
	  muonP4.push_back(p.p4());
	  //cout << "Pt: " << mcpt << "\tStatus: " << Status << "\n";
	}else if (abs(mcpid) == 15 ) {
	  if (mcpt > ptTauMX) ptTauMX=mcpt;
	  tauP4.push_back(p.p4());
	} 
      }
    } // end for loop over generated particles
    if (ptEleMX > 0.) h_mxElePt->Fill(ptEleMX); 
    if (ptMuMX > 0.) h_mxMuPt->Fill(ptMuMX); 
    if (ptTauMX > 0.) h_mxTauPt->Fill(ptTauMX); 

    plotInvMass(eleP4,h_dieleMass);
    plotInvMass(muonP4,h_dimuonMass);
    plotInvMass(tauP4,h_ditauMass);

  }
}

void GenEventAnalyzer::plotInvMass(vector<math::XYZTLorentzVector>& P4vec, TH1F* h) {

  vector<math::XYZTLorentzVector>::iterator iter1, iter2;

  if (P4vec.size() >1 ){
    for(iter1=P4vec.begin(); iter1<P4vec.end()-1; iter1++){
      math::XYZTLorentzVector m1 = *iter1;
      for(iter2=iter1+1; iter2<P4vec.end(); iter2++){
	math::XYZTLorentzVector m2 = *iter2;
	double mass = (m1+m2).mass();
	h->Fill(mass); 
      }// end loop2
    }// end loop1
  }

}

void GenEventAnalyzer::endJob() {
  /*
  if (m_file !=0)
    {
    //Write the histograms to the file.
    m_file->cd();
    h_evtCounter.Write();
    h_ptHat.Write();

    h_mxElePt.Write();
    h_mxTauPt.Write();
    h_mxMuPt.Write();

    h_dimuonMass.Write();

    delete m_file;
    m_file=0;
    }
  */
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventAnalyzer);
