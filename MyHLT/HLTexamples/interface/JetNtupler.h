#ifndef Examples_CompJets_h
#define Examples_CompJets_h
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <cmath>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/VertexReco/interface/VertexFwd.h" 

#include "DataFormats/Math/interface/deltaR.h"

class TFile;

#define MAXJ 130
#define MAXTRIG 500
// #define MAXVRT 50

static const size_t MAXVRT = 50;

typedef struct 
{
  int run;
  int lumi;
  int event;
} EventInfo;

typedef struct 
{
  void set(reco::PFJetCollection::const_iterator j, int i) 
  {
    pt[i]=j->pt();
    eta[i]=j->eta();
    phi[i]=j->phi();
    e[i]=j->energy();

    chf[i]  =j->chargedHadronEnergyFraction();
    nhf[i]  =(j->neutralHadronEnergy() + j->HFHadronEnergy())/j->energy();
    //nhf[i]  =j-> neutralHadronEnergyFraction ();
    cef[i]  =j-> chargedEmEnergyFraction ();
    nef[i]  =j-> neutralEmEnergyFraction ();
    nconstituents[i]  =j->getPFConstituents().size();
    //nch[i]=j->associatedTracks().size();

    id[i]=jetId(i);
 
  }
  bool jetId(int i)
  {
    if(nhf[i] > 0.99) return false;
    if(nef[i] > 0.99) return false;
    if(nconstituents[i]  <= 1) return false;
    if(fabs(eta[i])<2.5) {
      if(cef[i] > 0.99) return false;
      if(chf[i] == 0) return false;
      //if(nch[i]== 0) return false;
    }
    return true;
  }
  void reset()
  {
    for(int i=0;i<MAXJ;i++) {
      pt[i]=-99; eta[i]=-99; phi[i]=-99;e[i]=-99;
      chf[i]=-99; nhf[i]=-99; cef[i]=-99; nef[i]=-99; nch[i]=-99; nconstituents[i]=-99; 
      id[i]=false;
    }
  }
  float pt[MAXJ];
  float eta[MAXJ];
  float phi[MAXJ];
  float e[MAXJ];
  float chf[MAXJ];
  float nhf[MAXJ];
  float cef[MAXJ];
  float nef[MAXJ];
  float nch[MAXJ];
  float nconstituents[MAXJ];
  bool id[MAXJ];

} JetInfo;



class JetNtupler : public edm::EDAnalyzer {
public:
  JetNtupler( const edm::ParameterSet & );

private:
  void beginJob();
  void analyze( const edm::Event& , const edm::EventSetup& );
  void endJob();

  void getHLTResults(const edm::TriggerResults&, const edm::TriggerNames&);

  void bookTree();

  edm::Service<TFileService> fs;

  edm::InputTag HLTriggerResults,HLTJets_,RecoJets_, GenJets_, Rho_, Vertices_;
  bool Debug_,Monte_;

  TH1F *h_TriggerResults;

  // store hlt information in a map
  std::vector<bool> hlttrigs;
  std::map <std::string,bool> hltTriggerMap;
  std::map<std::string,bool>::iterator trig_iter;

  edm::TriggerNames triggerNames_;  // TriggerNames class

  int errCnt;
  bool HLTinit_;

  int nhJets,nrJets,ntrig;

  
  //Tree Variables
  TTree *HltJets;
  
  EventInfo EVENT;
  JetInfo hJets, rJets;

  int *triggerBits;
  float rho;

  int   NVrtx;
  float *VertexCand_x, *VertexCand_y, *VertexCand_z;
  int   *VertexCand_tracks;
  float *VertexCand_chi2;
  float *VertexCand_ndof;

  const int errMax(){return 20;}

};

#endif
