// JetNtupler.cc
// Description:  Example of simple EDAnalyzer for HLT.
// Author: L. Apanasevich
// Date:  28 - August - 2006
//
#include "MyHLT/HLTexamples/interface/JetNtupler.h"

using namespace edm;
using namespace reco;
using namespace std;

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings PFJetAlgorithm and GenJetAlgorithm.
JetNtupler::JetNtupler( const ParameterSet & cfg ) :
  HLTriggerResults( cfg.getParameter<InputTag>( "HLTriggerResults" ) ),
  HLTJets_( cfg.getParameter<InputTag>( "HLTJets" ) ),
  RecoJets_( cfg.getParameter<InputTag>( "RecoJets" ) ),
  GenJets_( cfg.getParameter<InputTag>( "GenJets" ) ),
  Rho_( cfg.getParameter<InputTag>( "Rho" ) ),
  Vertices_( cfg.getParameter<InputTag>( "Vertices" ) ),
  Debug_( cfg.getParameter<bool>( "Debug" ) ),
  Monte_( cfg.getParameter<bool>( "Monte" ) ),
  HLTinit_(false)
  {

  errCnt=0;
}

void JetNtupler::beginJob() {

  bookTree();

}

void JetNtupler::analyze( const Event& evt, const EventSetup& es ) {

  bool gotHLT=true;
  string errMsg("");

  Handle<TriggerResults> hltresults,hltresultsDummy;
  evt.getByLabel(HLTriggerResults,hltresults);
  if (! hltresults.isValid() ) { 
    gotHLT=false; errCnt+=1; errMsg=errMsg + "  -- No HLTriggerResults";
  }

  Handle<double> rhoHandle;
  evt.getByLabel(Rho_,rhoHandle);   
  rho= *rhoHandle;

  EVENT.run = evt.id().run();
  EVENT.lumi = evt.id().luminosityBlock();
  EVENT.event = evt.id().event();

  nhJets=0;
  nrJets=0;
  ntrig=0;
  NVrtx=0;

  hJets.reset();
  rJets.reset(); 
  for(int i=0; i < MAXTRIG; ++i){
    triggerBits[i]=0;
  }

  if (gotHLT) {
    const TriggerNames & triggerNames_ = evt.triggerNames(*hltresults);
    getHLTResults(*hltresults, triggerNames_);
    ntrig=hltresults->size();
    for(int i=0; i < ntrig; ++i){
      if (hltresults->accept(i)) triggerBits[i]=1;
    }
  }

  //Get the PFJet collections
  Handle<PFJetCollection> HLTJets,RecoJets;
  evt.getByLabel( HLTJets_, HLTJets );
  evt.getByLabel( RecoJets_, RecoJets );


  if (HLTJets.isValid()) { 
    //Loop over the PFJets and fill some histograms
    // std::cout << "Got Jets: " << HLTJets->size() << " " << RecoJets->size() << std::endl;


    for( PFJetCollection::const_iterator jet = HLTJets->begin(); jet != HLTJets->end(); ++ jet ) {

      double scale=1.;
      double jpt=scale*jet->pt();
      if (jpt>5. && nhJets<MAXJ){
	hJets.set(jet,nhJets);
	nhJets++;
      }

    }
  }else{
    errMsg=errMsg + "  -- No " + HLTJets_.label();
  }

  if (RecoJets.isValid()) { 
    for( PFJetCollection::const_iterator jet = RecoJets->begin(); jet != RecoJets->end(); ++ jet ) {
      double scale=1.;
      double jpt=scale*jet->pt();
      if (jpt>5. && nrJets<MAXJ){
	rJets.set(jet,nrJets);
	nrJets++;
      }
    }
  }else{
    errMsg=errMsg + "  -- No " + RecoJets_.label();
  }

  if (Monte_){
    // Compare with GenJets if requested
    Handle<GenJetCollection>  genJets;
    evt.getByLabel( GenJets_, genJets );

    if (genJets.isValid()){

      for( GenJetCollection::const_iterator gjet = genJets->begin(); gjet != genJets->end(); ++ gjet ) {
	//find the best matched jet
	
	double genpt =gjet->pt();
	double geneta=gjet->eta();
	double genphi=gjet->phi();

      }
    }else{
      errMsg=errMsg + "  -- No " + GenJets_.label();
    }
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

  Handle<VertexCollection> recoVertexs;
  evt.getByLabel( Vertices_, recoVertexs );

  if ( recoVertexs.isValid() ) {
    const reco::VertexCollection* vertexs = recoVertexs.product();
    reco::VertexCollection::const_iterator vertex_i;

    size_t size = std::min(MAXVRT, size_t(vertexs->size()) ); 
    NVrtx= size;

    int nVertexCand=0;
    if (Debug_)  std::cout << "Found " << vertexs->size() << " vertices" << std::endl;  
    for (vertex_i = vertexs->begin(); vertex_i != vertexs->end(); vertex_i++){
      if (nVertexCand>=NVrtx) break;
      VertexCand_x[nVertexCand] = vertex_i->x();
      VertexCand_y[nVertexCand] = vertex_i->y();
      VertexCand_z[nVertexCand] = vertex_i->z();
      VertexCand_tracks[nVertexCand] = vertex_i->tracksSize();
      VertexCand_chi2[nVertexCand] = vertex_i->chi2();
      VertexCand_ndof[nVertexCand] = vertex_i->ndof();
      if (Debug_) { 
        std::cout << "RECOVertex -- VX, VY VZ   = " 
                  << VertexCand_x[nVertexCand] << " "
                  << VertexCand_y[nVertexCand] << " "
                  << VertexCand_z[nVertexCand]
                  << std::endl;
        std::cout << "RECOVertex -- Ntracks, Chi2/Dof   = " 
                  << VertexCand_tracks[nVertexCand] << " "
                  << VertexCand_chi2[nVertexCand] << " / " << VertexCand_ndof[nVertexCand]
                  << std::endl;
      }
      nVertexCand++;
      
    }
  }

  HltJets->Fill();
}

void JetNtupler::getHLTResults( const edm::TriggerResults& hltresults,
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

void JetNtupler::endJob() {

  cout << "\n%%%%%%%%%%%%%%%%  Job Summary %%%%%%%%%%%%%%%%%\n\n" ;


}


void JetNtupler::bookTree() {

  // Initialize the tree
  // HltJets = new TTree("HltJets", "");
  HltJets = fs->make<TTree>("HltJets", "");

  HltJets->Branch("EVENT",  &EVENT,   "run/I:lumi/I:event/I");

  std::stringstream s;
  triggerBits = new int[MAXTRIG];
  HltJets->Branch("ntrig"		,  &ntrig	            ,  "ntrig/I");
  HltJets->Branch("triggerBits", triggerBits, "triggerBits[ntrig]/I"); 

  HltJets->Branch("nhJets"		,  &nhJets	            ,  "nhJets/I");
  HltJets->Branch("hJet_pt",hJets.pt ,"pt[nhJets]/F");
  HltJets->Branch("hJet_eta",hJets.eta ,"eta[nhJets]/F");
  HltJets->Branch("hJet_phi",hJets.phi ,"phi[nhJets]/F");
  HltJets->Branch("hJet_e",hJets.e ,"e[nhJets]/F");
  HltJets->Branch("hJet_chf",hJets.chf ,"chf[nhJets]/F");
  HltJets->Branch("hJet_nhf",hJets.nhf ,"nhf[nhJets]/F");
  HltJets->Branch("hJet_cef",hJets.cef ,"cef[nhJets]/F");
  HltJets->Branch("hJet_nef",hJets.nef ,"nef[nhJets]/F");
  HltJets->Branch("hJet_nch",hJets.nch ,"nch[nhJets]/F");
  HltJets->Branch("hJet_nconstituents",hJets.nconstituents ,"nconstituents[nhJets]");
  HltJets->Branch("hJet_id",hJets.id ,"id[nhJets]/b");

  HltJets->Branch("nrJets"		,  &nrJets	            ,  "nrJets/I");
  HltJets->Branch("rJet_pt",rJets.pt ,"pt[nrJets]/F");
  HltJets->Branch("rJet_eta",rJets.eta ,"eta[nrJets]/F");
  HltJets->Branch("rJet_phi",rJets.phi ,"phi[nrJets]/F");
  HltJets->Branch("rJet_e",rJets.e ,"e[nrJets]/F");
  HltJets->Branch("rJet_chf",rJets.chf ,"chf[nrJets]/F");
  HltJets->Branch("rJet_nhf",rJets.nhf ,"nhf[nrJets]/F");
  HltJets->Branch("rJet_cef",rJets.cef ,"cef[nrJets]/F");
  HltJets->Branch("rJet_nef",rJets.nef ,"nef[nrJets]/F");
  HltJets->Branch("rJet_nch",rJets.nch ,"nch[nrJets]/F");
  HltJets->Branch("rJet_nconstituents",rJets.nconstituents ,"nconstituents[nrJets]");
  HltJets->Branch("rJet_id",rJets.id ,"id[nrJets]/b");

  HltJets->Branch("rho",&rho,"rho/F");

  VertexCand_x          = new float[MAXVRT];
  VertexCand_y          = new float[MAXVRT];
  VertexCand_z          = new float[MAXVRT];
  VertexCand_tracks     = new int[MAXVRT];
  VertexCand_chi2       = new float[MAXVRT];
  VertexCand_ndof       = new float[MAXVRT];

  HltJets->Branch("NVrtx",& NVrtx,"NVrtx/I");

  TString br_recoVrtY = "recoVrtY";
  HltJets->Branch(br_recoVrtY,     VertexCand_y,          "recoVrtY[NVrtx]/F");

  TString br_recoVrtZ = "recoVrtZ";
  HltJets->Branch(br_recoVrtZ,     VertexCand_z,          "recoVrtZ[NVrtx]/F");

  TString br_recoVrtNtrk = "recoVrtNtrk";
  HltJets->Branch(br_recoVrtNtrk,  VertexCand_tracks,     "recoVrtNtrk[NVrtx]/I");

  TString br_recoVrtChi2 = "recoVrtChi2";
  HltJets->Branch(br_recoVrtChi2,  VertexCand_chi2,       "recoVrtChi2[NVrtx]/F");

  TString br_recoVrtNdof = "recoVrtNdof";
  HltJets->Branch(br_recoVrtNdof,  VertexCand_ndof,       "recoVrtNdof[NVrtx]/F");

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetNtupler);
