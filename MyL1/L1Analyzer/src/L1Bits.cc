// JetPlotsExample.cc
// Description:  Example of simple EDAnalyzer for jets.
// Author: Robert M. Harris
// Date:  28 - August - 2006
//

#include "MyL1/L1Analyzer/interface/L1Bits.h"

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
L1Bits::L1Bits( const ParameterSet & cfg ) {
  cout << " Beginning L1Jet Analysis " << endl;
  particleMapSource_= cfg.getParameter< edm::InputTag > ("particleMapSource");
  histogram        = cfg.getParameter<string>( "Histogram" );
  text_output      = cfg.getParameter<string>( "Outfile" );
  errCnt=0;
}

void L1Bits::beginJob( const EventSetup & ) {

  // Open the histogram file and book some associated histograms
  m_file=new TFile(histogram.c_str(),"RECREATE");
  evtCounter=new TH1F("EventCounter","Event Counter",5,0.,5.);
}

void L1Bits::analyze( const Event& evt, const EventSetup& es ) {

  string errMsg("");

  //Get the collections
  evtCounter->Fill(0.,1.);
  edm::Handle<l1extra::L1ParticleMapCollection> l1mapcoll;
  try {evt.getByLabel(particleMapSource_,l1mapcoll );} catch (...) { errMsg=errMsg + "  -- No L1 Map Collection";}

  
  L1Analysis(*l1mapcoll);

}

void L1Bits::L1Analysis(const l1extra::L1ParticleMapCollection& L1MapColl) {


  //cout << "%doL1Analysis -- Number of l1bits:   " << L1MapColl.size() << endl;
  //cout << "%doL1Analysis -- Number of l1bits:   " << l1extra::L1ParticleMap::kNumOfL1TriggerTypes << endl;

  int nacc=0;
  for (unsigned int itrig = 0; itrig != L1MapColl.size() ; ++itrig){
    const l1extra::L1ParticleMap& map = ( L1MapColl )[ itrig ] ;
    bool accept = map.triggerDecision();
    string trigName = map.triggerName();

    if (accept)nacc++;

    m_iter=m_bits.find(trigName);
    if (m_iter==m_bits.end()){
      typedef std::map<string,int>::value_type valType;
      m_bits.insert(valType(trigName,accept));
    }else{
      int nacc=m_iter->second+accept;
      m_iter->second=nacc;
    }
  }
  if (nacc>0) evtCounter->Fill(1.,1.);
}

void L1Bits::endJob() {

  double ntot=evtCounter->GetBinContent(1);
  double nacc=evtCounter->GetBinContent(2);

  ofstream ofile;
  ofile.open (text_output.c_str());

  cout << "Number of L1 Triggers: " << m_bits.size() << endl;
  ofile << "Number of L1 Triggers: " << m_bits.size() << endl;
  TH1F *h = new TH1F("TriggerBits","L1 Trigger Bits",m_bits.size(),0.,m_bits.size());

  m_iter = m_bits.begin();
  int ibin=0;
  while (m_iter != m_bits.end()){
    ibin++;
    cout << m_iter->first << ":\t" << m_iter->second << endl;
    ofile << m_iter->first << ":\t" << m_iter->second << endl;

    const char* trigName =  m_iter->first.c_str();
    h->GetXaxis()->SetBinLabel(ibin,trigName);
    h->SetBinContent(ibin,m_iter->second);
    h->SetBinError(ibin,sqrt(m_iter->second));
    ++m_iter;
  }

  ofile << "\n";
  ofile << "Number of Events Processed  : " << ntot << "\n";
  ofile << "Number of L1 Accepted Events: " << nacc << "\n";
  ofile << "                      Ratio : " << nacc/ntot << endl;

  ofile.close();

  //Write out the histogram file.
  m_file->Write();

}

void L1Bits::fillHist(const TString& histName, const Double_t& value, const Double_t& wt) {

  hid=m_HistNames.find(histName);
  if (hid==m_HistNames.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value,wt); 

}

void L1Bits::fill2DHist(const TString& histName, const Double_t& x,const Double_t& y,const Double_t& wt) {

  hid2D=m_HistNames2D.find(histName);
  if (hid2D==m_HistNames2D.end())
    std::cout << "%fillHist -- Could not find histogram with name: " << histName << std::endl;
  else
    hid2D->second->Fill(x,y,wt); 

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1Bits);
