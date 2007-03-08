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

void HLTJetAnalysis::fillHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt) {

  hid2D=m_HistNames2D.find(histName);
  if (hid2D==m_HistNames2D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid2D->second->Fill(x,y,wt); 

}

void HLTJetAnalysis::bookHistograms() {

  bookGeneralHistograms();

  bookHLTHistograms();

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
			   const CaloGeometry& geom) {

  // std::cout << " Beginning HLTJetAnalysis -- Pthat= " << _CKIN3 << ":" << _CKIN4 <<std::endl;

  m_file->cd();

  (&calojets) ? doCaloJets=true : doCaloJets=false;
  (&genjets) ? doGenJets=true : doGenJets=false;


  getHLTResults(hltresults);
  trig_iter=hltTriggerMap.find(_HLTPath);
  if (trig_iter==hltTriggerMap.end()){
    std::cout << "%HLTJetAnalysis -- Could not find trigger with pathname: " << _HLTPath << std::endl;
    return;
  }else{
    hlttrig=trig_iter->second;
  }

  getHLTParticleInfo(hltobj);


  fillHist("Nevents",0.0);
  if (hlttrig)fillHist("Nevents",1.0);

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
    if (accept) fillHist("HLTrigger",float(itrig));

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

    fillHist("CaloTowerEnergyEtaPhi",eta,phi,totEnergy);

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
