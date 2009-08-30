#include "MyL1/L1Analyzer/interface/L1Bits.h"

// Get the algorithm of the jet collections we will read from the .cfg file
// which defines the value of the strings CaloJetAlgorithm and GenJetAlgorithm.
L1Bits::L1Bits( const ParameterSet & cfg ) {
  cout << " Beginning L1Jet Analysis " << endl;
  l1GtRecordInputTag = cfg.getParameter<InputTag>( "L1GtRecordInputTag" );
  l1GtObjectMap      = cfg.getParameter<InputTag> ("L1GtObjectMapRecord");
  l1AlgoName         = cfg.getParameter<string> ("L1AlgoName");
  text_output        = cfg.getParameter<string>( "Outfile" );
  errCnt=0;
  initL1=false;
  initAnalysis=false;
  h_L1Results=0;
}

void L1Bits::beginJob( const EventSetup & ) {

  evtCounter=fs->make<TH1F>("EventCounter","Event Counter",5,0.,5.);

}

void L1Bits::analyze( const Event& evt, const EventSetup& es ) {

  string errMsg("");
  bool gotL1=true;

  evtCounter->Fill(0.);
  //Get the collections
  Handle<L1GlobalTriggerReadoutRecord>   gtRecord;
  Handle<L1GlobalTriggerObjectMapRecord> gtOMRec;

  evt.getByLabel(l1GtRecordInputTag, gtRecord);
  evt.getByLabel(l1GtObjectMap     , gtOMRec);

  if (!gtRecord.isValid()) {

    cout << "\nL1GlobalTriggerReadoutRecord with \n  "
	 << l1GtRecordInputTag << "\nnot found"
      "\n  --> returning false by default!\n" << endl;

    gotL1=false;
    
  }

  if (!gtOMRec.isValid()) {

    cout << "\nL1GlobalTriggerObjectMapRecord with \n  "
	 << l1GtObjectMap << "\nnot found"
      "\n  --> returning false by default!\n" << endl;

    gotL1=false;
    
  }

  if (gotL1) {
    evtCounter->Fill(1.);
    getL1Results(*gtRecord.product(), *gtOMRec.product());
    doL1Analysis();
  }
}

void L1Bits::doL1Analysis(){

    bool l1Fired=false;


    // fill histogram with L1 results
    if (!initAnalysis){
      initAnalysis=true;
      h_L1Results = fs->make<TH1F>( "h_L1Results", "L1 Trigger Results", numberTriggerBits, 0, numberTriggerBits );
      trig_iter=l1TriggerMap.begin();
      for ( unsigned int ibin = 0; ibin < numberTriggerBits; ++ ibin){
	trig_iter=l1TriggerMap.find(algoBitToName[ibin]);
	const char* trigName =  trig_iter->first.c_str();
	h_L1Results->GetXaxis()->SetBinLabel(ibin+1,trigName);
      }
    }

    for ( unsigned int ibin = 0; ibin < numberTriggerBits; ++ ibin){
      trig_iter=l1TriggerMap.find(algoBitToName[ibin]);
      bool accept = trig_iter->second;
      if (accept) h_L1Results->Fill(float(ibin));
    }


    // get results for particular L1 algorithm

    trig_iter=l1TriggerMap.find(l1AlgoName);
    if (trig_iter==l1TriggerMap.end()){
      cout << "Could not find trigger path with name: " << l1AlgoName << endl;
    }else{
      l1Fired=trig_iter->second;
    }
    if (l1Fired) evtCounter->Fill(2.);
}

void L1Bits::getL1Results(const L1GlobalTriggerReadoutRecord& gtRecord,
			const L1GlobalTriggerObjectMapRecord& gtOMRec
			) {

  const DecisionWord dWord = gtRecord.decisionWord();  // this will get the decision word *before* masking disabled bits
  numberTriggerBits= dWord.size();

  if (numberTriggerBits > nL1BitsMax) {
    cout << "number of Trigger Bits exceed max allowed in trigger array -- truncating" << "\n";
    numberTriggerBits = nL1BitsMax;
  }

  if ( !initL1){
    initL1=true;
    cout << "\n  Number of Trigger bits " << numberTriggerBits << "\n\n";
    cout << "\tBit \t L1 Algorithm " << endl;
    // get ObjectMaps from ObjectMapRecord
    const std::vector<L1GlobalTriggerObjectMap>& objMapVec =  gtOMRec.gtObjectMap();
    for (std::vector<L1GlobalTriggerObjectMap>::const_iterator itMap = objMapVec.begin();
	 itMap != objMapVec.end(); ++itMap) {
      // Get trigger bits
      int itrig = (*itMap).algoBitNumber();
      // Get trigger names
      algoBitToName[itrig] = (*itMap).algoName();
      cout << "\t" << itrig << "\t" << algoBitToName[itrig] << endl;      
    } // end of for loop    
  } // end of if

  for (unsigned int iBit = 0; iBit < numberTriggerBits; ++iBit) {  

    bool accept = dWord[iBit];
    //if (accept) h_L1Results2->Fill(float(iBit));

      // fill the trigger map
    typedef std::map<string,bool>::value_type valType;
    trig_iter=l1TriggerMap.find(algoBitToName[iBit]);
    if (trig_iter==l1TriggerMap.end())
      l1TriggerMap.insert(valType(algoBitToName[iBit],accept));
    else
      trig_iter->second=accept;
    }
}

void L1Bits::endJob() {

  double ntot=evtCounter->GetBinContent(1);
  double nacc=evtCounter->GetBinContent(2);

  if (h_L1Results){

    cout << "XXX" << endl;
    ofstream ofile;
    ofile.open (text_output.c_str());

    int nbins=h_L1Results->GetNbinsX(); 

    cout  << "Number of L1 Triggers: " << nbins << "\n\n";
    ofile << "Number of L1 Triggers: " << nbins << "\n\n";

    cout  << "\tL1 Algorithm \t\t # of Accepts" << "\n";
    cout  << "\t------------ \t\t ------------" << "\n";
    ofile << "\tL1 Algorithm \t\t # of Accepts" << "\n";
    ofile << "\t------------ \t\t ------------" << "\n";
    for (int ibin=0; ibin<nbins; ++ibin){
      float cont=h_L1Results->GetBinContent(ibin+1);
      //const char* trigName =  h_L1Results->GetXaxis()->GetBinLabel(ibin+1);
      string trigName = string (h_L1Results->GetXaxis()->GetBinLabel(ibin+1));
      
      if (!trigName.empty()){
	cout  << "\t" << trigName << ":\t" << cont << endl;
	ofile << "\t" << trigName << ":\t" << cont << endl;
      }
    }
    
    ofile << "\n";
    ofile << "Number of Events Processed  : " << ntot << "\n";
    ofile << "Number of L1 Accepted Events: " << nacc << "\n";
    ofile << "                      Ratio : " << nacc/ntot << endl;
    
    ofile.close();
  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1Bits);
