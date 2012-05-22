// CompJets.cc
// Description:  Example of simple EDAnalyzer for HLT.
// Author: L. Apanasevich
// Date:  28 - August - 2006
//
#include "MyHLT/HLTexamples/interface/CompJets.h"

using namespace edm;
using namespace reco;
using namespace std;

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings PFJetAlgorithm and GenJetAlgorithm.
CompJets::CompJets( const ParameterSet & cfg ) :
  HLTriggerResults( cfg.getParameter<InputTag>( "HLTriggerResults" ) ),
  PFJetCollection1_( cfg.getParameter<InputTag>( "PFJetCollection1" ) ),
  PFJetCollection2_( cfg.getParameter<InputTag>( "PFJetCollection2" ) ),
  GenJetCollection_( cfg.getParameter<InputTag>( "GenJetCollection" ) ),
  MyTrigger( cfg.getParameter<string>( "MyTrigger" ) ),
  HLTinit_(false)
  {

  errCnt=0;
}

void CompJets::beginJob() {

  bookHistograms();

}

void CompJets::analyze( const Event& evt, const EventSetup& es ) {

  bool gotHLT=true;
  bool myTrig=false;
  string errMsg("");

  Handle<TriggerResults> hltresults,hltresultsDummy;
  evt.getByLabel(HLTriggerResults,hltresults);
  if (! hltresults.isValid() ) { 
    gotHLT=false; errCnt+=1; errMsg=errMsg + "  -- No HLTriggerResults";
  }

  int iLumi = evt.luminosityBlock();
  fillHist("lumi",float(iLumi));


  if (gotHLT) {
    const TriggerNames & triggerNames_ = evt.triggerNames(*hltresults);
    getHLTResults(*hltresults, triggerNames_);
    trig_iter=hltTriggerMap.find(MyTrigger);

    if (MyTrigger == "NoTrigger"){
      myTrig=true;
    }else{
      if (trig_iter==hltTriggerMap.end()){
	cout << "Could not find trigger path with name: " << MyTrigger << endl;
      }else{
	myTrig=trig_iter->second;
      }
    }
  }

  fillHist("evtCounter",0.); // count number of events processed
  if (myTrig) fillHist("evtCounter",1.); // count number of events that fired my Trigger

  if (! myTrig) return;

  //Get the PFJet collections
  Handle<PFJetCollection> Jets1,Jets2;
  evt.getByLabel( PFJetCollection1_, Jets1 );
  evt.getByLabel( PFJetCollection2_, Jets2 );

  // Compare with GenJets if available
  Handle<GenJetCollection>  genJets;
  evt.getByLabel( GenJetCollection_, genJets );

  if (Jets1.isValid()) plotJetPtandEta(*Jets1,"coll1");
  if (Jets2.isValid()) plotJetPtandEta(*Jets2,"coll2");


  if (Jets1.isValid() && Jets2.isValid() ) { 
    //Loop over the PFJets and fill some histograms
    int jetInd = 0;
    // std::cout << "Got Jets: " << Jets1->size() << " " << Jets2->size() << std::endl;
    for( PFJetCollection::const_iterator jet1 = Jets1->begin(); jet1 != Jets1->end(); ++ jet1 ) {
      // std::cout << "CALO JET #" << jetInd << std::endl << cal->print() << std::endl;

      // double scale = corrector->correction(jet1->p4());
      double scale=1.;
      double jpt=scale*jet1->pt();

      fillHist("PFJetPt",jpt );
      if (jpt>30 && jpt <=40) fillHist("jetEta_pt30_40",jet1->eta());

      if (jetInd == 0){  // plot leading jet kinematic distributions
	fillHist("PFJetPt_Leading", jpt);
	fillHist("PFJetEta_Leading", jet1->eta());
	fillHist("PFJetPhi_Leading", jet1->phi());
      
	if (myTrig) fillHist("PFJetPt_Trig", jpt );
	jetInd++;
      }

      //find the best matched jet
      double drmin=99.;
      PFJetCollection::const_iterator miter;
      for( PFJetCollection::const_iterator jet2 = Jets2->begin(); jet2 != Jets2->end(); ++ jet2 ) {
	double dr=deltaR(jet1->eta(),jet1->phi(),jet2->eta(),jet2->phi());
	if (dr<drmin) {
	  drmin=dr;
	  miter=jet2;
	}	  
      }
	
      fillHist("drMin",drmin,1.);
      if (jpt>10.) fillHist("drMin_pt10",drmin,1.);
      if (jpt>20.) fillHist("drMin_pt20",drmin,1.);
      if (jpt>40.) fillHist("drMin_pt40",drmin,1.);
      fill2DHist("drMinVSpT",jpt,drmin,1.);
      
      if (drmin<drMatch()){
	double rat=miter->pt()/jpt;
	fill3DHist("Reco2overReco1_vsReco1",jpt,jet1->eta(),rat,1.);
	fill3DHist("Reco2overReco1_vsReco2",miter->pt(),miter->eta(),rat,1.);

	fill3DHist("Reco1overReco2_vsReco1",jpt,jet1->eta(),1./rat,1.);
	fill3DHist("Reco1overReco2_vsReco2",miter->pt(),miter->eta(),1./rat,1.);
	
	if (jpt>30 && jpt<=40){
	  if (fabs(jet1->eta())<1.1){
	    fillHist("jetResp_pt30_40_eta1",rat);
	  }else if (fabs(jet1->eta())<2.5){
	    fillHist("jetResp_pt30_40_eta2",rat);
	  }else {
	    fillHist("jetResp_pt30_40_eta3",rat);
	  }
	}
      }
    }

    if (genJets.isValid()){

      for( GenJetCollection::const_iterator gjet = genJets->begin(); gjet != genJets->end(); ++ gjet ) {
      //find the best matched jet

	double genpt =gjet->pt();
	double geneta=gjet->eta();
	double genphi=gjet->phi();

	double drmin1=99.,drmin2=99.;
	PFJetCollection::const_iterator miter1,miter2;
	for( PFJetCollection::const_iterator jet1 = Jets1->begin(); jet1 != Jets1->end(); ++ jet1 ) {
	  double dr=deltaR(geneta,genphi,jet1->eta(),jet1->phi());
	  if (dr<drmin1) {
	    drmin1=dr;
	    miter1=jet1;
	  }	  
	}

	for( PFJetCollection::const_iterator jet2 = Jets2->begin(); jet2 != Jets2->end(); ++ jet2 ) {
	  double dr=deltaR(geneta,genphi,jet2->eta(),jet2->phi());
	  if (dr<drmin2) {
	    drmin2=dr;
	    miter2=jet2;
	  }	  
	}
	//fill histograms

	fill2DHist("GenReco1_drMinVSpT",genpt,drmin1,1.);
	fill2DHist("GenReco2_drMinVSpT",genpt,drmin2,1.);
	
	if (drmin1<drMatch()){
	  double rat=miter1->pt()/genpt;
	  fill3DHist("GenReco1_jetResp3D",genpt,geneta,rat,1.);
	}

	if (drmin2<drMatch()){
	  double rat=miter2->pt()/genpt;
	  fill3DHist("GenReco2_jetResp3D",genpt,geneta,rat,1.);
	}

      }
    }else{
      errMsg=errMsg + "  -- No " + GenJetCollection_.label();
    }

  }else{
    if (! Jets1.isValid())
      errMsg=errMsg + "  -- No " + PFJetCollection1_.label();
    if (! Jets2.isValid())
      errMsg=errMsg + "  -- No " + PFJetCollection2_.label();
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

void CompJets::getHLTResults( const edm::TriggerResults& hltresults,
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

void CompJets::plotJetPtandEta(const PFJetCollection& Jets, const string& suffix){

  for( PFJetCollection::const_iterator jet = Jets.begin(); jet != Jets.end(); ++ jet ) {

      double pt=jet->pt();
      double eta=jet->eta();

      fill2DHist("JetPtAndEta_" + suffix ,pt,eta,1.);
  }
}

void CompJets::endJob() {

  cout << "\n%%%%%%%%%%%%%%%%  Job Summary %%%%%%%%%%%%%%%%%\n\n" ;

  //if (h_evtCounter){
  //  double ntot=h_evtCounter->GetBinContent(1);
  //  cout  << "\tNumber of events processed: " << int(ntot) << "\n\n";
  //}

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

TH1F* CompJets::Book1dHist(const char* name, const char* title, Int_t nbins, Double_t xmin, Double_t xmax,bool DoSumw2=true){


  TH1F *h= fs->make<TH1F>(name, title, nbins, xmin, xmax); 
  if (DoSumw2) h->Sumw2();

  return h;
}

TH2F* CompJets::Book2dHist(const char* name, const char* title, Int_t nbinsx, Double_t xmin, Double_t xmax, 
		 Int_t nbinsy, Double_t ymin, Double_t ymax ,bool DoSumw2=true){


  TH2F *h= fs->make<TH2F>(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax); 
  if (DoSumw2) h->Sumw2();

  return h;
}

TH3F* CompJets::Book3dHist(const char* name, const char* title, Int_t nbinsx, Double_t xmin, Double_t xmax , 
			   Int_t nbinsy, Double_t ymin, Double_t ymax , Int_t nbinsz, Double_t zmin, 
			   Double_t zmax , bool DoSumw2=true){

  TH3F *h= fs->make<TH3F>(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax); 
  if (DoSumw2) h->Sumw2();

  return h;
}

void CompJets::bookHistograms() {

  TString hname,htitle;

  hname="evtCounter"; htitle="Event Counter";
  m_HistNames[hname] =  Book1dHist(hname, htitle, 10, -0.5, 9.5, false );

  hname="lumi"; htitle="Luminosity Sections";
  m_HistNames[hname] =  Book1dHist( hname, htitle, 1000, 0.,1000., false );

  int njbins=200;
  float minj=0., maxj=1000.;
  hname="PFJetPt"; htitle="PFJet p_{T}";
  m_HistNames[hname] = Book1dHist(hname, htitle, njbins,minj,maxj);

  hname="PFJetPt_Leading"; htitle="p_{T} of leading PFJets";
  m_HistNames[hname] =   Book1dHist( hname, htitle, njbins,minj,maxj);

  hname="PFJetPt_Trig"; htitle="p_{T} of leading PFJets -- Triggered";
  m_HistNames[hname] =  Book1dHist( hname, htitle, njbins,minj,maxj);

  hname="PFJetEta_Leading"; htitle="#eta of leading PFJets";
  m_HistNames[hname] = Book1dHist(hname, htitle , 102, -5.1, 5.1 );

  hname="PFJetPhi_Leading"; htitle="#phi of leading PFJets";
  // m_HistNames[hname] = fs->make<TH1F>(hname, htitle , 50, -M_PI, M_PI );
  m_HistNames[hname] = Book1dHist(hname, htitle , 50, -M_PI, M_PI );


  // DeltaR
  hname="drMin"; htitle="Min dr between jets";
  m_HistNames[hname]=Book1dHist(hname, htitle, 100, 0., 1.); 

  hname="drMin_pt10"; htitle="Min dr between jets -- p_{T} > 10 GeV";
  m_HistNames[hname]=Book1dHist(hname, htitle, 100, 0., 1.); 

  hname="drMin_pt20"; htitle="Min dr between jets -- p_{T} > 20 GeV";
  m_HistNames[hname]=Book1dHist(hname, htitle, 100, 0., 1.); 

  hname="drMin_pt40"; htitle="Min dr between jets -- p_{T} > 40 GeV";
  m_HistNames[hname]=Book1dHist(hname, htitle, 100, 0., 1.); 

  hname="drMinVSpT"; htitle="Min dr between rjet1 and rjet2";
  m_HistNames2D[hname]=Book2dHist(hname, htitle, 100, 0. ,1000., 100, 0., 1.); 

  hname="GenReco1_drMinVSpT"; htitle="Min dr between rjet1 and genjets";
  m_HistNames2D[hname]=Book2dHist(hname, htitle, 100, 0. ,1000., 100, 0., 1.); 

  hname="GenReco2_drMinVSpT"; htitle="Min dr between rjet2 an genjets";
  m_HistNames2D[hname]=Book2dHist(hname, htitle, 100, 0. ,1000., 100, 0., 1.); 

  //hists vs eta
  int neta=20;
  double etamin=-5.,etamax=5.;

  hname="jetEta_pt30_40"; htitle="jet #eta 30 < p_{T} <  40";
  m_HistNames[hname]=Book1dHist(hname, htitle, neta, etamin, etamax);


  //jet Resp vs eta

  int njr=20;
  double jrmin=0.,jrmax=2.;

  hname="JetPtAndEta_coll1"; htitle="Jet p_{T} and #eta; collection 1 ";  // pt, eta, resp
  m_HistNames2D[hname]=Book2dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax);
  hname="JetPtAndEta_coll2"; htitle="Jet p_{T} and #eta; collection 2 ";  // pt, eta, resp
  m_HistNames2D[hname]=Book2dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax);

  hname="Reco2overReco1_vsReco1"; htitle="p_{T} (coll2) / p_{T} (coll1) vs jet collection 1 pt,eta";  // pt, eta, resp
  m_HistNames3D[hname]=Book3dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);
  hname="Reco2overReco1_vsReco2"; htitle="p_{T} (coll2) / p_{T} (coll1) vs jet collection 2 pt,eta";
  m_HistNames3D[hname]=Book3dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);

  hname="Reco1overReco2_vsReco1"; htitle="p_{T} (coll1) / p_{T} (coll2) vs jet collection 1 pt,eta";  // pt, eta, resp
  m_HistNames3D[hname]=Book3dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);
  hname="Reco1overReco2_vsReco2"; htitle="p_{T} (coll1) / p_{T} (coll2) vs jet collection 2 pt,eta";
  m_HistNames3D[hname]=Book3dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);

  hname="GenReco1_jetResp3D"; htitle="Jet response Coll1/Gen pt,eta";
  m_HistNames3D[hname]=Book3dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);
  hname="GenReco2_jetResp3D"; htitle="Jet response Coll2/Gen pt,eta";
  m_HistNames3D[hname]=Book3dHist(hname, htitle, 100, 0. ,1000., neta, etamin, etamax, njr, jrmin, jrmax);


  hname="jetResp_pt30_40_eta1"; htitle="jetResp |#eta|<1.1  30 < p_{T} <40";
  m_HistNames[hname]=Book1dHist(hname, htitle, njr, jrmin, jrmax);

  hname="jetResp_pt30_40_eta2"; htitle="jetResp 1.1<|#eta|<2.5  30 < p_{T} <40";
  m_HistNames[hname]=Book1dHist(hname, htitle, njr, jrmin, jrmax);

  hname="jetResp_pt30_40_eta3"; htitle="jetResp |#eta|>2.5  30 < p_{T} <40";
  m_HistNames[hname]=Book1dHist(hname, htitle, njr, jrmin, jrmax);



}

void CompJets::fillHist(const TString& histName, const Double_t& value, const Double_t& wt) {

  hid=m_HistNames.find(histName);
  if (hid==m_HistNames.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value,wt); 

}

void CompJets::fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt) {

  hid2D=m_HistNames2D.find(histName);
  if (hid2D==m_HistNames2D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid2D->second->Fill(x,y,wt); 

}

void CompJets::fill3DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& z,const Double_t& wt) {

  hid3D=m_HistNames3D.find(histName);
  if (hid3D==m_HistNames3D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid3D->second->Fill(x,y,z,wt); 

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CompJets);
