// JetPlotsExample.cc
// Description:  Example of simple EDAnalyzer for jets.
// Author: Robert M. Harris
// Date:  28 - August - 2006
//

#include "MyL1/L1Analyzer/interface/L1JetPlots.h"

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
L1JetPlots::L1JetPlots( const ParameterSet & cfg ) {
  cout << " Beginning L1Jet Analysis " << endl;
  CaloJetAlgorithm = cfg.getParameter<string>( "CaloJetAlgorithm" );
  GenJetAlgorithm  = cfg.getParameter<string>( "GenJetAlgorithm" );
  recmet_          = cfg.getParameter< string > ("recmet");
  genmet_          = cfg.getParameter< string > ("genmet");
  l1CollectionsTag_= cfg.getParameter< InputTag > ("l1collections");
  histogram        = cfg.getParameter<string>( "Histogram" );

  errCnt=0;
}

void L1JetPlots::beginJob( const EventSetup & ) {

  // Open the histogram file and book some associated histograms
  m_file=new TFile(histogram.c_str(),"RECREATE");
  ptCal =  TH1F( "ptCal",  "p_{T} of leading CaloJets", 50, 0, 1000 );
  etaCal = TH1F( "etaCal", "#eta of leading CaloJets", 50, -3, 3 );
  phiCal = TH1F( "phiCal", "#phi of leading CaloJets", 50, -M_PI, M_PI );
  ptGen =  TH1F( "ptGen",  "p_{T} of leading GenJets", 50, 0, 1000 );
  etaGen = TH1F( "etaGen", "#eta of leading GenJets", 50, -3, 3 );
  phiGen = TH1F( "phiGen", "#phi of leading GenJets", 50, -M_PI, M_PI );

  int nbins=50;
  Double_t min=0.,max=150.;

  MetPt = TH1F( "MetPt", "Reconstructed MET ", nbins, min, max );
  genMetPt = TH1F( "genMetPt", "Generated MET ", nbins, min, max );

  TString hname,htitle;

  hname="L1JetCollSize"; htitle="L1 Jet Collection Size";
  m_HistNames[hname] = new TH1F( hname , htitle  , 10, -0.5, 9.5 );

  Int_t nptbins=40;
  Double_t ptmin=0,ptmax=200.;
  hname="L1JetPt"; htitle="L1 Jet p_{T} -- Central";
  m_HistNames[hname] = new TH1F( hname , htitle  , nptbins, ptmin, ptmax );
  m_HistNames[hname]->Sumw2();

  Int_t nphibins=20;
  Double_t phimin=-M_PI,phimax=M_PI;
  hname="L1JetPhi"; htitle="L1 Jet #phi -- Central";
  m_HistNames[hname] = new TH1F( hname , htitle  , nphibins, phimin, phimax );
  m_HistNames[hname]->Sumw2();

  Int_t netabins=80;
  Double_t etamin=-5,etamax=5.;
  hname="L1JetEta"; htitle="L1 Jet #eta -- Central";
  m_HistNames[hname] = new TH1F( hname , htitle  , netabins, etamin, etamax );  m_HistNames[hname]->Sumw2();

  hname="L1DeltaR_Calo"; htitle="Delta R -- L1 and CaloJets";
  m_HistNames[hname] = new TH1F( hname , htitle  , 100, 0., 10. );
  m_HistNames[hname]->Sumw2();

  hname="L1DeltaR_Gen"; htitle="Delta R -- L1 and GenJets";
  m_HistNames[hname] = new TH1F( hname , htitle  , 100, 0., 10. );
  m_HistNames[hname]->Sumw2();

  hname="L1PtOverCaloPt_Calo";htitle=" L1 Jet p_{T} over CaloJet p_{T}";
  m_HistNames[hname] = new TH1F( hname , htitle  , 40, 0., 2.0 );
  m_HistNames[hname]->Sumw2();

  hname="L1PtOverCaloPt_Gen";htitle=" L1 Jet p_{T} over GenJet p_{T}";
  m_HistNames[hname] = new TH1F( hname , htitle  , 40, 0., 2.0 );
  m_HistNames[hname]->Sumw2();

  Int_t n2dbins=50;
  hname="L1etaVSJeteta_Calo"; htitle="L1 Jet #eta vs CaloJet #eta -- Central";
  m_HistNames2D[hname] = new TH2F( hname , htitle  , n2dbins, etamin, etamax, n2dbins, etamin, etamax );
  hname="L1etaVSJeteta_Gen"; htitle="L1 Jet #eta vs GenJet #eta -- Central";
  m_HistNames2D[hname] = new TH2F( hname , htitle  , n2dbins, etamin, etamax, n2dbins, etamin, etamax );


  hname="L1phiVSJetphi_Calo"; htitle="L1 Jet #phi vs CaloJet #phi -- Central";
  m_HistNames2D[hname] = new TH2F( hname , htitle  , n2dbins, phimin, phimax, n2dbins, phimin, phimax );
  hname="L1phiVSJetphi_Gen"; htitle="L1 Jet #phi vs GenJet #phi -- Central";
  m_HistNames2D[hname] = new TH2F( hname , htitle  , n2dbins, phimin, phimax, n2dbins, phimin, phimax );

  hname="L1ptVSJetpt_Calo"; htitle="L1 Jet p_{T} vs CaloJet p_{T} -- Central";
  m_HistNames2D[hname] = new TH2F( hname , htitle  , n2dbins, ptmin, ptmax, n2dbins, ptmin, ptmax );
  hname="L1ptVSJetpt_Gen"; htitle="L1 Jet p_{T} vs GenJet p_{T} -- Central";
  m_HistNames2D[hname] = new TH2F( hname , htitle  , n2dbins, ptmin, ptmax, n2dbins, ptmin, ptmax );

}

void L1JetPlots::analyze( const Event& evt, const EventSetup& es ) {

  string errMsg("");
  doCaloJets=true; doGenJets=true; doCaloMET=true; doGenMET=true, doL1Jets=true;

  //Get the collections
  Handle<CaloJetCollection> caloJets;
  Handle<GenJetCollection> genJets;
  Handle<CaloMETCollection> recmet;
  Handle<GenMETCollection>  genmet;
  Handle<l1extra::L1JetParticleCollection> l1jets;

  try { evt.getByLabel( CaloJetAlgorithm, caloJets );} catch (...) { errMsg=errMsg + "  -- No RecJets"; doCaloJets=false;}
  try { evt.getByLabel( GenJetAlgorithm, genJets );} catch (...) { errMsg=errMsg + "  -- No GenJets"; doGenJets=false;}
  try { evt.getByLabel( recmet_,recmet );} catch (...) { errMsg=errMsg + "  -- No RecMET"; doCaloMET=false;}
  try { evt.getByLabel( genmet_,genmet );} catch (...) { errMsg=errMsg + "  -- No GenMET"; doGenMET=false;}

  InputTag L1JetTag(edm::InputTag(l1CollectionsTag_.label(),"Central"));
  try { evt.getByLabel(L1JetTag,l1jets); } catch (...) 
    { errMsg=errMsg + "  -- No L1Jets with name: " + L1JetTag.label() ; doL1Jets=false;}


  if (doCaloJets){
    //Loop over the two leading CaloJets and fill some histograms
    int jetInd = 0;
    for( CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end() && jetInd<2; ++ cal ) {
      // std::cout << "CALO JET #" << jetInd << std::endl << cal->print() << std::endl;
      ptCal.Fill( cal->pt() );
      etaCal.Fill( cal->eta() );
      phiCal.Fill( cal->phi() );
      jetInd++;
    }
  }

  if (doGenJets){
    //Loop over the two leading GenJets and fill some histograms
    int jetInd = 0;
    for( GenJetCollection::const_iterator gen = genJets->begin(); gen != genJets->end() && jetInd<2; ++ gen ) {
      // std::cout << "GEN JET #" << jetInd << std::endl << gen->print() << std::endl;
      ptGen.Fill( gen->pt() );
      etaGen.Fill( gen->eta() );
      phiGen.Fill( gen->phi() );
      jetInd++;
    }
  }

  if (doCaloMET){
    for( CaloMETCollection::const_iterator met = recmet->begin(); met != recmet->end() ; ++ met ) {
      MetPt.Fill(met->pt());
    }
  }

  if (doGenMET){
    for( GenMETCollection::const_iterator met = genmet->begin(); met != genmet->end() ; ++ met ) {
      genMetPt.Fill(met->pt());
    }
  }

  if (doL1Jets && doCaloJets && doGenJets)
    L1Analysis(*caloJets,*genJets,*l1jets);

  if ((errMsg != "") && (errCnt < errMax())){
    errCnt=errCnt+1;
    errMsg=errMsg + ".";
    cout << "%L1JetPlots-Warning" << errMsg << endl;
    if (errCnt == errMax()){
      errMsg="%L1JetPLots-Warning -- Maximum error count reached -- No more messages will be printed.";
      cout << errMsg << endl;
    }
  }

}

void L1JetPlots::L1Analysis(const reco::CaloJetCollection& caloJets,
				const reco::GenJetCollection& genJets,
				const l1extra::L1JetParticleCollection& l1Jets) {

  //CalJetIter cIter;
  //GenJetIter gIter;

  //cout << "%doL1Analysis -- Number of l1jets:   " << l1Jets.size() << endl;
  //cout << "%doL1Analysis -- Number of calojets: " << caloJets.size() << endl;


  fillHist("L1JetCollSize",l1Jets.size());

  for(l1extra::L1JetParticleCollection::const_iterator l1 = l1Jets.begin(); l1 != l1Jets.end(); ++l1) {
    
    Double_t pt_l1=l1->pt();
    Double_t eta_l1=l1->eta();
    Double_t phi_l1=l1->phi();
    
    fillHist("L1JetPt",pt_l1);
    fillHist("L1JetPhi",phi_l1);
    fillHist("L1JetEta",eta_l1);
    
    // match L1 Jets with Jets
    mtchL1(eta_l1,phi_l1,pt_l1,caloJets,"_Calo");
    mtchL1(eta_l1,phi_l1,pt_l1,genJets,"_Gen");
  }
}

template <typename T> void L1JetPlots::mtchL1(const Double_t& eta_l1, 
						  const Double_t& phi_l1, 
						  const Double_t& pt_l1, 
						  const T& jets,
						  const TString& WhichJets) {

  // kick out if Jet collection does not exist
  if (! &jets) return;

  double etaCut=2.5;
  double drCut=0.5;


  typedef typename T::const_iterator iter;
  iter mIter;

  float rmin(99.);
  for ( iter jiter=jets.begin(); jiter!=jets.end(); jiter++) {
    
    if (fabs(jiter->eta()) < etaCut){

      Double_t pt_jet=jiter->pt();
      Double_t eta_jet=jiter->eta();
      Double_t phi_jet=jiter->phi();
      
      float dr=deltaR(eta_jet,phi_jet,eta_l1,phi_l1);
      if (pt_l1>20 && pt_jet>10) fillHist("L1DeltaR" + WhichJets,dr);
      
      if(dr<rmin){rmin=dr;mIter=jiter;}
      
    }  
  }

  if (rmin < drCut){
    fill2DHist("L1etaVSJeteta" + WhichJets,mIter->eta(),eta_l1,1.);
    fill2DHist("L1phiVSJetphi" + WhichJets,mIter->phi(),phi_l1,1.);
    fill2DHist("L1ptVSJetpt" + WhichJets,mIter->pt(),pt_l1,1.);
    
    if (mIter->pt()>75.) {
      fillHist("L1PtOverCaloPt" + WhichJets,pt_l1/(mIter->pt()));
    }
  }
}

void L1JetPlots::endJob() {

  //Write out the histogram file.
  m_file->Write();

}

void L1JetPlots::fillHist(const TString& histName, const Double_t& value, const Double_t& wt) {

  hid=m_HistNames.find(histName);
  if (hid==m_HistNames.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value,wt); 

}

void L1JetPlots::fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt) {

  hid2D=m_HistNames2D.find(histName);
  if (hid2D==m_HistNames2D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid2D->second->Fill(x,y,wt); 

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1JetPlots);
