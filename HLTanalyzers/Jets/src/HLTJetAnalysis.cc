#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <stdlib.h>
#include <string.h>

#include "HLTanalyzers/Jets/interface/HLTJetAnalysis.h"
#include "HLTanalyzers/Jets/interface/JetUtil.h"

typedef CaloJetCollection::const_iterator CalJetIter;
typedef GenJetCollection::const_iterator GenJetIter;

HLTJetAnalysis::HLTJetAnalysis() {
  m_file=0; // set to null
  evtCounter=0;

  //set parameter defaults 
  _Monte=false;
  _Debug=false;
  _EtaMin=-5.2;
  _EtaMax=5.2;
  _HistName="test.root";
  _HLTPath="xxx";

  // initialize some variables
  _CKIN3=-999.;
  _CKIN4=-999.;
  hlttrig=false;
  hltInfoExists=false;
    
}

/*  Setup the analysis to put the histograms into HistoFile. */
void HLTJetAnalysis::setup(const edm::ParameterSet& pSet) {


  edm::ParameterSet myJetParams = pSet.getParameter<edm::ParameterSet>("RunParameters") ;
  vector<std::string> parameterNames = myJetParams.getParameterNames() ;
  
  for ( vector<std::string>::iterator iParam = parameterNames.begin();
	iParam != parameterNames.end(); iParam++ ){
    if  ( (*iParam) == "Monte" ) _Monte =  myJetParams.getParameter<bool>( *iParam );
    else if ( (*iParam) == "Debug" ) _Debug =  myJetParams.getParameter<bool>( *iParam );
    else if ( (*iParam) == "EtaMin" ) _EtaMin =  myJetParams.getParameter<double>( *iParam );
    else if ( (*iParam) == "EtaMax" ) _EtaMax =  myJetParams.getParameter<double>( *iParam );
    else if ( (*iParam) == "HLTPath" ) _HLTPath =  myJetParams.getParameter<string>( *iParam );
    else if ( (*iParam) == "HistogramFile" ) _HistName =  myJetParams.getParameter<string>( *iParam );
  }

  cout << "---------- Input Parameters ---------------------------" << endl;
  cout << "  Monte:  " << _Monte << endl;    
  cout << "  Debug:  " << _Debug << endl;    
  cout << "  EtaMin: " << _EtaMin << endl;    
  cout << "  EtaMax: " << _EtaMax << endl;    
  cout << "  HLT Path: " << _HLTPath << endl;    
  cout << "  Output histograms written to: " << _HistName << std::endl;
  cout << "-------------------------------------------------------" << endl;  
  // open the histogram file

  m_file=new TFile(_HistName.c_str(),"RECREATE");
  m_file->cd();
  bookHistograms();

}

void HLTJetAnalysis::fillHist(const TString& histName, const Double_t& value, const Double_t& wt) {

  hid=m_HistNames.find(histName);
  if (hid==m_HistNames.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value,wt); 

}

void HLTJetAnalysis::fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt) {

  hid2D=m_HistNames2D.find(histName);
  if (hid2D==m_HistNames2D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid2D->second->Fill(x,y,wt); 

}

void HLTJetAnalysis::bookHistograms() {

  bookGeneralHistograms();

  bookHLTHistograms();
  bookL1Histograms();

  bookCaloTowerHists();

  if (_Monte) bookJetHistograms("Gen");
  bookJetHistograms("Calo");

  if (_Monte) bookMetHists("Gen");
  bookMetHists("Calo");

  if (_Monte) bookMCParticles();
  
}
/* **Analyze the event** */
void HLTJetAnalysis::analyze( const CaloJetCollection& calojets,
			      const GenJetCollection& genjets,
			      const CaloMETCollection& recmets,
			      const GenMETCollection& genmets,
			      const CaloTowerCollection& caloTowers,
			      const HepMC::GenEvent mctruth,
			      const HLTFilterObjectWithRefs& hltobj,
			      const edm::TriggerResults& hltresults,
			      const l1extra::L1JetParticleCollection& l1jets,
			      const CaloGeometry& geom) {

  // std::cout << " Beginning HLTJetAnalysis -- Pthat= " << _CKIN3 << ":" << _CKIN4 <<std::endl;

  m_file->cd();

  (&calojets) ? doCaloJets=true : doCaloJets=false;
  (&genjets) ? doGenJets=true : doGenJets=false;
  (&l1jets) ? doL1Jets=true : doL1Jets=false;


  float ptHat=-1.;
  if (_Monte){
    ptHat=mctruth.event_scale();
    //cout << "Pt of hard scatter: " << ptHat << endl;
  }

  fillHist("Nevents",0.0);

  evtTriggered=false;
  getHLTResults(hltresults);
  if (evtTriggered)fillHist("Nevents",1.0); // fill Nevents histogram if at least one trigger fired
  trig_iter=hltTriggerMap.find(_HLTPath);
  if (trig_iter==hltTriggerMap.end()){
    std::cout << "%HLTJetAnalysis -- Could not find trigger with pathname: " << _HLTPath << std::endl;
    //return;
  }else{
    hlttrig=trig_iter->second;
  }

  getHLTParticleInfo(hltobj);


  if (hlttrig)fillHist("Nevents",2.0); // fill Nevents histogram if desired trigger fired

  // Make a copy, so that you can sort

  CaloJetCollection mycalojets;
  if (doCaloJets) {
    mycalojets=calojets;
    std::sort(mycalojets.begin(),mycalojets.end(),PtGreater());
  }

  GenJetCollection mygenjets;
  if (doGenJets) {
    mygenjets=genjets;
    std::sort(mygenjets.begin(),mygenjets.end(),PtGreater());
  }

  // fill calojet and genjet hists 
  fillJetHists(mycalojets,"Calo");
  if (_Monte) fillJetHists(mygenjets,"Gen");

  // fill recmet and genjet hists
  fillMetHists(recmets,"Calo");
  if (_Monte) fillMetHists(genmets,"Gen");

  // fill CaloTower hists
  fillCaloTowerHists(caloTowers);

  if (doL1Jets && doCaloJets && doGenJets)
    L1Analysis(mycalojets,mygenjets,l1jets);
  
  if (_Monte) fillMCParticles(mctruth);

}

void HLTJetAnalysis::getHLTResults(const edm::TriggerResults& hltResults) {

  hltInfoExists=false;
  if (! &hltResults) return;
  hltInfoExists=true;
  
  // TriggerResults is derived from from HLTGlobalStatus
  int ntrigs=hltResults.size();  
  if (_Debug) std::cout << "%getHLTResults --  Number of HLT Triggers: " << ntrigs << std::endl;

  TString hname,htitle;
  hname="HLTrigger"; htitle="HLT Trigger Status";
  const TString hnames[]={ "Calopt_","Calopt_Leading_","Genpt_","Genpt_Leading_" };

  bool justBooked=false;
  // check if this histogram exists, if not book it
  hid=m_HistNames.find(hname);
  if (hid==m_HistNames.end()){
    std::cout << "Trigger Summary Histogram does not exist -- creating" << std::endl;
    m_HistNames[hname]= new TH1F(hname,htitle,ntrigs,0.,ntrigs);
    justBooked=true;
  }
  //hid->second->Fill(value,wt); 

  for (int itrig = 0; itrig != ntrigs; ++itrig){
    string trigName=hltResults.name(itrig);
    if ( justBooked ) {
      m_HistNames[hname]->GetXaxis()->SetBinLabel(itrig+1,trigName.c_str());
      // book jet histograms
      for (int ilabel=0; ilabel != 4; ++ilabel){
	TString h_jet=hnames[ilabel]+"Untriggered";
	hid=m_HistNames.find(h_jet);
	if (hid==m_HistNames.end())
	  std::cout << "%fillHist -- Could not find histogram with name: " << h_jet << std::endl;
	else
	  {
	    TString hnewname= hnames[ilabel] + trigName.c_str();
	    //	    TriggeredJetHists.push_back(hnewname);
	    TH1F* h =  (TH1F*)hid->second->Clone();
	    h->SetName(hnewname);
	    h->Reset();
	    m_HistNames[hnewname]=h;
	    //hid->second->Fill(value,wt); 
	  }
      }
    }

    bool accept=hltResults.accept(itrig);
    if (accept) {
      fillHist("HLTrigger",float(itrig));
      evtTriggered=true;
    }

    if (_Debug){
      std::cout << "%getHLTResults --  HLTTrigger(" << itrig << "): " 
		<< trigName << " = " << accept << std::endl;
    }

    typedef std::map<string,bool>::value_type valType;
    trig_iter=hltTriggerMap.find(trigName);
    if (trig_iter==hltTriggerMap.end())
      hltTriggerMap.insert(valType(trigName,accept));
    else
      trig_iter->second=accept;
  }
}

void HLTJetAnalysis::getHLTParticleInfo(const HLTFilterObjectWithRefs& hltObj) {

  if (! &hltObj) return;
  
  int mod=-1,path=-1,npart=-1;

  mod = hltObj.module();
  path = hltObj.path();
  //  npart = hltObj.numberParticles();
  npart = hltObj.size();

  fillHist("Npart",npart);
  if (hlttrig) fillHist("Npart_triggered",npart);

  for (int ipart = 0; ipart != npart; ++ipart){
    const edm::RefToBase<Candidate> ref_ = hltObj.getParticleRef(ipart);
    Double_t pt=ref_->pt();
    Double_t eta=ref_->eta();
    fillHist("HLT_pt",pt); 
    fillHist("HLT_eta",eta); 
  }

  if (_Debug){
    std::cout << "%getHLTParticleInfo --  HLTobj module: " 
	      << mod << "   path: " << path << "   Npart:" << npart << std::endl;
  }

}

void HLTJetAnalysis::bookGeneralHistograms() {

  TString hname,htitle;

  hname="Nevents"; htitle="Number of events";
  m_HistNames[hname]= new TH1F(hname,htitle,7,0.0,7.0);

}

void HLTJetAnalysis::bookJetHistograms(const TString& prefix) {

  TString hname;
  TString htitle;

  std::ostringstream ch_eta; ch_eta << etaBarrel();
  std::ostringstream ch_etamin; ch_etamin << _EtaMin;
  std::ostringstream ch_etamax; ch_etamax << _EtaMax;

  TString h_EtaRng= ch_etamin.str() + " < #eta < " + ch_etamax.str();
  
  // book rec and gen jet histograms
  Int_t netbins=40, nengbins=100;
  Double_t etmin=0.,etmax=400.,engmin=0.,engmax=500.;

  hname=prefix + "et"; htitle=prefix+" Jet E_{T} -- " + h_EtaRng;
  m_HistNames[hname]= new TH1F(hname,htitle,netbins,etmin,etmax);
  hname=prefix + "pt"; htitle=prefix+" Jet p_{T} -- " + h_EtaRng;
  m_HistNames[hname]= new TH1F(hname,htitle,netbins,etmin,etmax);

  // Book Untriggered Jet Histograms for HLT studies -- triggered jet hists
  // booked in routine getHLTResults
  hname=prefix + "pt_Untriggered"; htitle=prefix+" Jet p_{T} " + h_EtaRng;
  m_HistNames[hname]= new TH1F(hname,htitle,netbins,etmin,etmax);  
  hname=prefix + "pt_Leading_Untriggered"; htitle=prefix+" Jet p_{T} -- Leading " + h_EtaRng;
  m_HistNames[hname]= new TH1F(hname,htitle,netbins,etmin,etmax);  
  //end

  hname=prefix + "energy"; htitle=prefix+" Jet Energy -- " + h_EtaRng;
  m_HistNames[hname]= new TH1F(hname,htitle,nengbins,engmin,engmax);

  float deltaEta=0.1,deltaPhi=0.1;
  Int_t netabins=int((_EtaMax-_EtaMin)/deltaEta);
  Int_t nphibins=int(2.*M_PI/deltaPhi);

  hname=prefix + "phi" ;  htitle=prefix+" Jet #phi -- " + h_EtaRng;
  m_HistNames[hname] = new TH1F(hname,htitle,nphibins,-M_PI,M_PI);

  hname=prefix + "eta"; htitle=prefix+" Jet #eta -- " + h_EtaRng;
  m_HistNames[hname] = new TH1F(hname,htitle,netabins,_EtaMin,_EtaMax);

  TString PtCut1=" -- p_{T} > 10 GeV/c";
  TString PtCut2=" -- p_{T} > 25 GeV/c";
  TString PtCut3=" -- p_{T} > 50 GeV/c";

  hname=prefix + "eta1"; htitle=prefix+" Jet #eta -- " + h_EtaRng + PtCut1;
  m_HistNames[hname] = new TH1F(hname,htitle,netabins,_EtaMin,_EtaMax);
  hname=prefix + "eta2"; htitle=prefix+" Jet #eta -- " + h_EtaRng + PtCut2;
  m_HistNames[hname] = new TH1F(hname,htitle,netabins,_EtaMin,_EtaMax);
  hname=prefix + "eta3"; htitle=prefix+" Jet #eta -- " + h_EtaRng + PtCut3;
  m_HistNames[hname] = new TH1F(hname,htitle,netabins,_EtaMin,_EtaMax);


  hname=prefix + "et_Barrel";  htitle=prefix+" Jet E_{T} -- |#eta| < " + ch_eta.str();
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);
  hname=prefix + "pt_Barrel";  htitle=prefix+" Jet p_{T} -- |#eta| < " + ch_eta.str();
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);

  hname=prefix + "pt_Endcap";  htitle=prefix+" Jet p_{T} -- 1.4 < |#eta| < 2.5 ";
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);
  hname=prefix + "pt_Forward";  htitle=prefix+" Jet p_{T} -- |#eta| > 2.5 ";
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);


  hname=prefix + "etmax"; htitle=prefix+" Max Jet E_{T} -- |#eta| < " + ch_eta.str();
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);
  hname=prefix + "ptmax"; htitle=prefix+" Max Jet p_{T} -- |#eta| < " + ch_eta.str();
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);


  hname=prefix + "etmax2"; htitle=prefix+" Max Jet E_{T} -- " + h_EtaRng;
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);
  hname=prefix + "ptmax2"; htitle=prefix+" Max Jet p_{T} -- " + h_EtaRng;
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);

  hname=prefix + "etmax2_hlt"; htitle=prefix+" Max Jet E_{T} -- HLT -- " + h_EtaRng;
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);
  hname=prefix + "ptmax2_hlt"; htitle=prefix+" Max Jet p_{T} -- HLT --" + h_EtaRng;
  m_HistNames[hname] = new TH1F(hname,htitle,netbins,etmin,etmax);


}

template <typename T> void HLTJetAnalysis::fillJetHists(const T& jets, const TString& prefix) {

  // kick out if Jet collection does not exist
  if (! &jets) return;

  typedef typename T::const_iterator iter;
  
  Double_t maxEt=0.,maxPt=0.;
  Double_t maxEt2=0.,maxPt2=0.;

  int ijet=0;
  for ( iter i=jets.begin(); i!=jets.end(); i++) {

    Double_t jetEng = i->energy();
    Double_t jetEt = i->et();
    Double_t jetPt = i->pt();
    Double_t jetEta = i->eta();
    Double_t jetPhi = i->phi();

    if (jetEta > _EtaMin && jetEta < _EtaMax){

      ijet++; 
      bool leading=(ijet==1); // get leading jet in eta range for rate histograms
      
      fillHist(prefix + "energy",jetEng);
      fillHist(prefix + "et",jetEt);
      fillHist(prefix + "pt",jetPt);

      // HLT histograms
      if (hltInfoExists){
	fillHist(prefix + "pt_Untriggered",jetPt);
	if (leading) fillHist(prefix + "pt_Leading_Untriggered",jetPt);
	//loop over triggers
	std::map<std::string,bool>::const_iterator titer=hltTriggerMap.begin();
	while (titer != hltTriggerMap.end()){
	  TString hname=prefix + "pt_" + titer->first;
	  TString hname_leading=prefix + "pt_Leading_" + titer->first;
	  bool triggered=titer->second;
	  if (triggered) {
	    fillHist(hname,jetPt);
	    if (leading) fillHist(hname_leading,jetPt);
	  }
	  ++titer;
	}

      }

      fillHist(prefix + "eta",jetEta);
      if (jetPt >10.)fillHist(prefix + "eta1",jetEta);
      if (jetPt >25.)fillHist(prefix + "eta2",jetEta);
      if (jetPt >50.)fillHist(prefix + "eta3",jetEta);

      fillHist(prefix + "phi",jetPhi);

      if (jetEt > maxEt2) maxEt2 = jetEt;
      if (jetPt > maxPt2) maxPt2 = jetPt;

      if (fabs(jetEta) < etaBarrel()){
	fillHist(prefix + "et_Barrel",jetEt);
	fillHist(prefix + "pt_Barrel",jetPt);
	if (jetEt > maxEt) maxEt = jetEt;
	if (jetPt > maxPt) maxPt = jetPt;
      }else if (fabs(jetEta) < etaEndcap()){
	fillHist(prefix + "pt_Endcap",jetPt);
      }else{
	fillHist(prefix + "pt_Forward",jetPt);
      }
    }
  }
  fillHist(prefix + "etmax",maxEt);
  fillHist(prefix + "ptmax",maxPt);

  fillHist(prefix + "etmax2",maxEt2);
  fillHist(prefix + "ptmax2",maxPt2);

  if (hlttrig) {
    fillHist(prefix + "etmax2_hlt",maxEt2);
    fillHist(prefix + "ptmax2_hlt",maxPt2);
  }

}

void HLTJetAnalysis::bookHLTHistograms() {

  TString hname,htitle;

  hname="Npart"; htitle="Number of HLT particles";
  m_HistNames[hname]= new TH1F(hname,htitle,7,0.0,7.0);

  hname="Npart_triggered"; htitle="Number of HLT particles -- HLT Satisfied";
  m_HistNames[hname]= new TH1F(hname,htitle,7,0.0,7.0);

  hname="HLT_pt"; htitle="p_{T} Spectrum of HLT particles";
  m_HistNames[hname]= new TH1F(hname,htitle,100,0.0,300.0);

  hname="HLT_eta"; htitle="#eta Spectrum of HLT particles";
  m_HistNames[hname]= new TH1F(hname,htitle,26,-5.2,5.2);

}

void HLTJetAnalysis::bookL1Histograms() {

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
  hname="L1JetEta2030"; htitle="L1 Jet #eta -- Central -- 20 < L1 p_{T} < 30";
  m_HistNames[hname] = new TH1F( hname , htitle  , netabins, etamin, etamax );  m_HistNames[hname]->Sumw2();
  hname="L1JetEta3040"; htitle="L1 Jet #eta -- Central -- 30 < L1 p_{T} < 40";
  m_HistNames[hname] = new TH1F( hname , htitle  , netabins, etamin, etamax );  m_HistNames[hname]->Sumw2();
  hname="L1JetEta4050"; htitle="L1 Jet #eta -- Central -- 40 < L1 p_{T} < 50";
  m_HistNames[hname] = new TH1F( hname , htitle  , netabins, etamin, etamax );  m_HistNames[hname]->Sumw2();
  hname="L1JetEta5060"; htitle="L1 Jet #eta -- Central -- 50 < L1 p_{T} < 60";
  m_HistNames[hname] = new TH1F( hname , htitle  , netabins, etamin, etamax );  m_HistNames[hname]->Sumw2();
  hname="L1JetEta6080"; htitle="L1 Jet #eta -- Central -- 60 < L1 p_{T} < 80";
  m_HistNames[hname] = new TH1F( hname , htitle  , netabins, etamin, etamax );  m_HistNames[hname]->Sumw2();
  hname="L1JetEta80100"; htitle="L1 Jet #eta -- Central -- 80 < L1 p_{T} < 100";
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

void HLTJetAnalysis::bookCaloTowerHists() {

  TString hname;
  int netbins=50;
  float etmin=0.; float etmax=50.;

  hname="CaloTowerEt";
  m_HistNames[hname] = new TH1F(hname,"CaloTower E_{T}",netbins,etmin,etmax);

  hname="CaloTowerEta";
  m_HistNames[hname] = new TH1F(hname,"CaloTower #eta",110,-5.5,5.5);


  hname="CaloTowerEnergy";
  m_HistNames[hname] = new TH1F(hname,"CaloTower Energy",100,0.0,100.);


  hname="CaloTowerEmEnergy";
  m_HistNames[hname] = new TH1F(hname,"CaloTower EmEnergy",100,0.0,100.);

  hname="CaloTowerHadEnergy";
  m_HistNames[hname] = new TH1F(hname,"CaloTower HadEnergy",100,0.0,100.);

  hname="CaloTowerOuterEnergy";
  m_HistNames[hname] = new TH1F(hname,"CaloTower Outer Energy",100,0.0,100.);

  hname="CaloTowerEnergyEtaPhi";
  //  TH2F* f2d = new TH2F(hname,"CaloTower Outer Energy",100,0.0,100.,100,-3,+3);
  m_HistNames2D[hname] = new TH2F(hname,"CaloTower Energy",110,-5.5,5.5,72,-M_PI,+M_PI);


}
void HLTJetAnalysis::fillCaloTowerHists(const CaloTowerCollection& caloTowers) {

  // kick out if there is no CaloTower collection
  if (! &caloTowers) return;

  for ( CaloTowerCollection::const_iterator tower=caloTowers.begin(); 
	tower!=caloTowers.end(); tower++) {

    Double_t et=tower->et();
    Double_t eta=tower->eta();
    Double_t phi=tower->phi();

    if (et<1.) continue; 

    Double_t  totEnergy= tower->energy();
    Double_t  emEnergy= tower->emEnergy();
    Double_t  hadEnergy= tower->hadEnergy();
    Double_t  outerEnergy= tower->outerEnergy();


    fillHist("CaloTowerEt",et);
    fillHist("CaloTowerEta",eta);
    fillHist("CaloTowerEnergy",totEnergy);
    fillHist("CaloTowerEmEnergy",emEnergy);
    fillHist("CaloTowerHadEnergy",hadEnergy);
    fillHist("CaloTowerOuterEnergy",outerEnergy);

    fill2DHist("CaloTowerEnergyEtaPhi",eta,phi,totEnergy);

  }
}
void HLTJetAnalysis::bookMetHists(const TString& prefix) {

  TString hname;
  TString htitle;

  hname=prefix + "num";
  htitle = prefix+" Number of MET objects";
  m_HistNames[hname] = new TH1I(hname,htitle,10,0.,10.);

  hname=prefix+"etMiss";
  htitle=prefix+" Missing Et";
  m_HistNames[hname] = new TH1F(hname,htitle,500,0.0,500.);

  hname=prefix+"etMissX";
  htitle=prefix+" Missing Et-X";
  m_HistNames[hname] = new TH1F(hname,htitle,2000,-1000.0,1000.);

  hname=prefix+"etMissPhi";
  htitle=prefix+" Phi of Missing Et";
  m_HistNames[hname] = new TH1F(hname,htitle,100,-M_PI,M_PI);


  hname=prefix+"sumEt";
  htitle=prefix+" Sum Et";
  m_HistNames[hname] = new TH1F(hname,htitle,1400,0.0,14000.);

}
template <typename T> void HLTJetAnalysis::fillMetHists(const T& mets, const TString& prefix) {

  // kick out if Met collection does not exist
  if (! &mets) return;

  Int_t metnum=mets.size();
  fillHist(prefix + "num",metnum);

  typedef typename T::const_iterator iter;
  for ( iter met=mets.begin(); met!=mets.end(); met++) {

    Double_t mEt=met->pt();
    Double_t sumEt=met->sumEt();
    Double_t mEtPhi=met->phi();
    Double_t mEtX=met->px();

    fillHist(prefix + "etMiss",mEt);
    fillHist(prefix + "sumEt",sumEt);
    fillHist(prefix + "etMissX",mEtX);
    fillHist(prefix + "etMissPhi",mEtPhi);
  }
}

void HLTJetAnalysis::bookMCParticles(){

  TString hname;
  TString htitle;
 
  const int imax=1;
  const int jmax=1;


  for(int i=0;i<imax;i++){
    std::ostringstream oi; oi << i;
 
    for(int j=0;j<jmax;++j){
      std::ostringstream oj; oj << j;
 
      istringstream ints(oj.str());

      int k;
      ints>>k;
      cout << k << endl;

      hname="VertexZ"+oi.str()+oi.str();
      m_HistNames[hname] = new TH1F(hname,hname,100,-50.,50.);

      hname="VertexX"+oi.str()+oj.str();
      m_HistNames[hname] = new TH1F(hname,hname,200,-5.,5.);

      hname="VertexY"+oi.str()+oj.str();
      m_HistNames[hname] = new TH1F(hname,hname,200,-5.,5.);

      hname="Pt"+oi.str()+oj.str();
      m_HistNames[hname] = new TH1F(hname,hname,500,0.0,500.);

      hname="Pid"+oi.str()+oj.str();
      m_HistNames[hname] = new TH1F(hname,hname,10000,0.0,10000.);

    }
  }
}

void HLTJetAnalysis::fillMCParticles(const HepMC::GenEvent mctruth){

  // return for null mctruth collection
  if (! &mctruth) return;

  for (HepMC::GenEvent::particle_const_iterator partIter = mctruth.particles_begin(); partIter != mctruth.particles_end();
         ++partIter) {

    // Find the end vertex
    //   for (HepMC::GenEvent::vertex_const_iterator vertIter = mctruth.vertices_begin();
    //   vertIter != mctruth.vertices_end();
    //    ++vertIter) {
       CLHEP::HepLorentzVector creation = (*partIter)->CreationVertex();
       CLHEP::HepLorentzVector momentum = (*partIter)->Momentum();
       HepPDT::ParticleID id = (*partIter)->particleID();  // electrons and positrons are 11 and -11
       //   cout << "MC particle id " << id.pid() << ", creationVertex " << creation << " cm, initialMomentum " << momentum << " GeV/c" << endl;   
       fillHist("Pid00",id.pid()); 
       fillHist("VertexX00",creation.x());  
       fillHist("VertexY00",creation.y());  
       fillHist("VertexZ00",creation.z());  

       fillHist("Pt00",momentum.perp());
  }
}

void HLTJetAnalysis::extractPtHat(edm::EventSetup const& isetup) {
  edm::Service<edm::ConstProductRegistry> reg;

  // initialize pthat variables
  _CKIN3=-1.; _CKIN4=-1.;

  bool foundPars(false);

  const bool tracked = true;

  //std::vector<std::string> PythiaVString;
  edm::ParameterSet PythiaPars;
  // Loop over provenance of products in registry.
  for (edm::ProductRegistry::ProductList::const_iterator it = reg->productList().begin();
       it != reg->productList().end(); ++it) {
    edm::BranchDescription desc = it->second;

    /*
    std::cout << "  ModuleLabel:       " << desc.moduleLabel() << "\n"
	      << "  processName:       " << desc.processName() << "\n"
	      << "  className:         " << desc.className() << "\n"
	      << "  friendlyClassName: " << desc.friendlyClassName() << "\n"
	      << std::endl;

    */

    if (desc.friendlyClassName()=="edmHepMCProduct" && desc.moduleLabel()=="source" ) {
      //std::cout << "\tFound: " << desc.friendlyClassName() << std::endl;
      edm::ParameterSet result = getParameterSet(desc.psetID());
      // std::cout << result << std::endl; // dumps the full set

      vector<string> names;
      result.getParameterSetNames(names, tracked);
      vector<string>::const_iterator it = names.begin();
      vector<string>::const_iterator end = names.end();
      //cout << "Number of Pythia Parameter Sets: " << names.size() << endl; 
      for( ; it != end; ++it ){
	if ( *it =="PythiaParameters" ){
	  PythiaPars=result.getParameter<edm::ParameterSet>(*it);
	  foundPars=true;
	  goto FOUND;
	}
      }
    }
  }
 FOUND:
  if (foundPars){
    vector<string> names = PythiaPars.getParameterNames();
    vector<string>::const_iterator it = names.begin();
    vector<string>::const_iterator end = names.end();
    //cout << "Number of Pythia Parameter Sets2: " << names.size() << endl; 
    for( ; it != end; ++it ){
      //cout << " Name: " << *it << endl;

      std::vector<std::string> PythiaVString = 
	PythiaPars.getParameter<std::vector<std::string> >(*it);

      //cout << " CCLA -- Size of pythiaVstring: " << PythiaVString.size() << endl;
      for (unsigned int ipar = 0 ; ipar != PythiaVString.size() ; ++ipar){
	string pypar=PythiaVString[ipar];
	//cout << "Par[" << ipar << "]: " <<  pypar << endl;
	if (! pypar.compare(0,4,"CKIN")){
	  //cout << "!!!!!!!!!!!!!!!!!!!!!!!!  " << pypar <<endl;
	  if (! pypar.compare(0,7,"CKIN(3)")){
	    int i1=pypar.find_first_of("=")+1;
	    int i2=pypar.find_first_of(" ");
	    _CKIN3 = atof (pypar.substr(i1,i2-i1+1).c_str());
	    //cout << " Value: " << pypar.substr(i1,i2-i1+1) << endl; 
	  }
	  if (! pypar.compare(0,7,"CKIN(4)")){
	    int i1=pypar.find_first_of("=")+1;
	    int i2=pypar.find_first_of(" ");
	    _CKIN4 = atof (pypar.substr(i1,i2-i1+1).c_str());
	    //cout << " Value: " << pypar.substr(i1,i2-i1+1) << endl; 
	  }
	}
      }
    }
  }

}

void HLTJetAnalysis::L1Analysis(const reco::CaloJetCollection& caloJets,
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
    
    if (pt_l1 >= 20. && pt_l1 < 30. ) fillHist("L1JetEta2030",eta_l1);
    if (pt_l1 >= 30. && pt_l1 < 40. ) fillHist("L1JetEta3040",eta_l1);
    if (pt_l1 >= 40. && pt_l1 < 50. ) fillHist("L1JetEta4050",eta_l1);
    if (pt_l1 >= 50. && pt_l1 < 60. ) fillHist("L1JetEta5060",eta_l1);
    if (pt_l1 >= 60. && pt_l1 < 80. ) fillHist("L1JetEta6080",eta_l1);
    if (pt_l1 >= 80. && pt_l1 < 100. ) fillHist("L1JetEta80100",eta_l1);

    // match L1 Jets with Jets
    mtchL1(eta_l1,phi_l1,pt_l1,caloJets,"_Calo");
    mtchL1(eta_l1,phi_l1,pt_l1,genJets,"_Gen");
  }
}

template <typename T> void HLTJetAnalysis::mtchL1(const Double_t& eta_l1, 
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
      
      float dr=radius(eta_jet,phi_jet,eta_l1,phi_l1);
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

void HLTJetAnalysis::dummyAnalyze(
			   const CaloGeometry& geom) {
  std::cout << "Inside dummyAnalyse routine" << std::endl;
}

/* Finalization (close files, etc) */
void HLTJetAnalysis::done() {
  std::cout << "Closing up.\n";

  if (m_file!=0) { // if there was a histogram file...
    m_file->Write(); // write out the histrograms
    delete m_file; // close and delete the file
    m_file=0; // set to zero to clean up
  }
}
