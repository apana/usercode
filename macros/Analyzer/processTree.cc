/////////////////////////////////////////////////////////////////////////////////////////////////
//
//          Macro calculating overlaps between different triggers, L1 trigger individual- and
//          pure-rates, taking into account background and signal type samples
//
//          Send comments to bargassa@cern.ch
//
// Inputs : * Nfil : Number of files of different processes to be read
//          * ProcFil : Vector of files of different processes to be read (full path).
//                      These files are the outputs of HLTrigger/HLTanalyzers
//          * xsec : Cross-sections of different processes to be read (same order)
//          * Nqcd : Number of QCD pt-hat bins to be read
//          * ILumi : Instanteneous Luminosity
//          * Ntrig : Number of triggers
//          * trigname : Vector of triggers
//
// Outputs : * Individual rate of each trigger (as if alone)
//           * Pure rate (additional rate given already present triggers)
//             taking into account overlaps in each process, including each QCD pt-hat bin
//           * Total (pure) rate
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include "processTree.h"

using namespace std;


void processTree(string input_file, string _HistName, int nevts, bool doXSweighting=true, bool skimmedFiles=false){


  bool plotIt=false;

  cout << "\n Beginning program execution" << endl;

  if (readFilelist(input_file) < 0 ) {
    cout << "There was a problem reading the filelists" << endl;
    return ;
  }

  loadHLTNames();
  bookHistograms();

  unsigned int istart=0; unsigned int ifin=100;
  //ifin=2;
  if (ifin > ProcFil.size()) ifin=ProcFil.size();
  loadChain(istart,ifin);

  if (doXSweighting) cout << "   FileName                      Number of Entries    XSWeight       Skim Efficiency"  <<endl;
  // Loop over input files 
  for (unsigned int ip = istart; ip != ifin; ++ip){

    setBranchAdds(ip); // SetBranchAddress's
    
    Int_t   nevent = TabChain[ip]->GetEntries();
    TString filename=ProcFil[ip];
    int mcSample=mcSample_bin.at(ip);

    int n=nevent;
    if (nevts>0) n=nevts;
   
    for ( int index = 0; index < n; ++index ) { // loop over entries

      TabChain[ip]->GetEntry(index);

      double wt=1.;
      if (doXSweighting){ // get cross section weights and skim efficiencies

	double xswt=getXSWeight(mcpthat,mcSample);
	double skimeff=1.;
	if (skimmedFiles) skimeff=getSkimEff(mcpthat,mcSample);

	wt=(xswt*skimeff)/n;
	if (index == 0){
	  cout << filename << "\t" << n << "\t"  << xswt << "\t"  << skimeff << endl;
	}
      }
      // fill rate histograms
      for (int i=0; i<nhlt; i++ ){
	Double_t x=i+0.1;
	string hltname=hltnames[i];
	if (hltTriggers[hltname]>0) hlt_rates->Fill(x,wt);
      }


      if (i_etm20>0)met_rates->Fill(0.,wt);
      if (i_etm30>0)met_rates->Fill(1.,wt);
      if (i_etm40>0)met_rates->Fill(2.,wt);
      if (i_etm50>0)met_rates->Fill(3.,wt);
      if (hltTriggers[hlt1met]>0)met_rates->Fill(4.,wt);
      if (i_etm30>0 && metpt>50.) met_rates->Fill(5.,wt);
      if (i_etm30>0 && metpt>60.) met_rates->Fill(6.,wt);
      if (i_etm30>0 && metpt>70.) met_rates->Fill(7.,wt);
      if (i_etm30>0 && metpt>75.) met_rates->Fill(8.,wt);

      met_rat->Fill(metpt/metsum);
      if (metpt/metsum > .6){
      	//continue;
      }
      
      // ************** Level 1 Jet Information **********************
      nl1jetc->Fill(nl1cenjet);
      nl1jetf->Fill(nl1forjet);
      nl1jett->Fill(nl1taujet);

      double maxl1cenpt=0.,maxl1forpt=0.,maxl1taupt=0.,maxl1pt=0;
      for (int i=0;i<nl1cenjet; ++i){
	//printf("i=%d value=%f\n",i,l1cenjetpt[i]);
	if (l1cenjetpt[i] > maxl1cenpt) maxl1cenpt=l1cenjetpt[i];
      }
      if (maxl1cenpt>0.) mx_l1cenpt->Fill(maxl1cenpt);

      for (int i=0;i<nl1forjet; ++i){
	//printf("i=%d value=%f\n",i,l1forjetpt[i]);
	if (l1forjetpt[i] > maxl1forpt) maxl1forpt=l1forjetpt[i];
      }
      if (maxl1forpt>0.) mx_l1forpt->Fill(maxl1forpt);

      for (int i=0;i<nl1taujet; ++i){
	//printf("i=%d value=%f\n",i,l1taujetpt[i]);
	if (l1taujetpt[i] > maxl1taupt) maxl1taupt=l1taujetpt[i];
      }
      if (maxl1taupt>0.) mx_l1taupt->Fill(maxl1taupt);

      maxl1pt=maxl1cenpt;
      if (maxl1forpt>maxl1pt) maxl1pt=maxl1forpt;
      if (maxl1taupt>maxl1pt) maxl1pt=maxl1taupt;
      
      //cout << "l1 max pt " << maxl1pt << " " << i_l1bit150 << endl;

      h_genmet->Fill(genmetpt,wt);
      // L1 Missing ET and HT info
      //if (genmetpt > 20) cout << genmetpt << "\t"  << metpt << "\t" << l1met << endl;
      h_l1met->Fill(l1met,wt);
      if (genmetpt < 5.) h_l1met_3->Fill(l1met,wt);
      h_l1metvsMet->Fill(metpt,l1met,wt);

      float dphi=deltaPhi(metphi,l1metphi);
      if (abs(dphi) < 1.)h_l1met_1->Fill(l1met,wt);

      if ( metpt < 10 ){
	h_l1_metDphi_1->Fill(dphi,wt);
	if (abs(dphi) < 1.)h_l1_metRat_1->Fill(l1met/metpt,wt);
      }else if ( metpt < 20 ){
	h_l1_metDphi_2->Fill(dphi,wt);
	if (abs(dphi) < 1.)h_l1_metRat_2->Fill(l1met/metpt,wt);
      }else if ( metpt < 30 ){
	h_l1_metDphi_3->Fill(dphi,wt);
	if (abs(dphi) < 1.)h_l1_metRat_3->Fill(l1met/metpt,wt);
      }else if ( metpt < 50 ){
	h_l1_metDphi_4->Fill(dphi,wt);
	if (abs(dphi) < 1.)h_l1_metRat_4->Fill(l1met/metpt,wt);
      }else {
	h_l1_metDphi_5->Fill(dphi,wt);
	if (abs(dphi) < 1.)h_l1_metRat_5->Fill(l1met/metpt,wt);
      }

      h_l1sumet->Fill(l1mettot);
      h_l1sumht->Fill(l1mettothad);
      if (i_l1htt200 > 0) h_l1sumht_l1htt200->Fill(l1mettothad);

      // ************** Reconstructed Jet Info ***********************

      //cout << "Number of Jets in the event: " << njetcal << endl;
      double maxcaljetpt=0;
      double htsumjets=0;
      //double htsumpx=0., htsumpy=0.;
      //cout << "zzz: " << njetcal << endl;
      for (int i = 0; i < njetcal; ++i){
	//printf("j=%d value=%f\n",i,jetpt[i]);
	if (caljeteta[i]>mineta && caljeteta[i]< maxeta){
	  h_caljetpt->Fill(caljetpt[i],wt);
	  if (caljetpt[i] > maxcaljetpt) maxcaljetpt=caljetpt[i];
	  if (caljetpt[i] > 5.){ // calculate sumet from jets directly
	    htsumjets+=caljetpt[i];
	  }
	}
      }
      sumhtjets_Untriggered ->Fill(htsumjets,wt);

      if (maxcaljetpt>0.){
	mx_caljetpt_Untriggered->Fill(maxcaljetpt,wt);
	if (i_l1bit15>0) mx_caljetpt_l15->Fill(maxcaljetpt,wt);
	if (i_l1bit30>0) mx_caljetpt_l30->Fill(maxcaljetpt,wt);
	if (i_l1bit50>0) mx_caljetpt_l50->Fill(maxcaljetpt,wt);
	if (i_l1bit70>0) mx_caljetpt_l70->Fill(maxcaljetpt,wt);
	if (i_l1bit100>0) mx_caljetpt_l100->Fill(maxcaljetpt,wt);
	if (i_l1bit150>0) mx_caljetpt_l150->Fill(maxcaljetpt,wt);
	if (i_l1bit200>0) mx_caljetpt_l200->Fill(maxcaljetpt,wt);

	if (maxl1pt>14.) mx_caljetpt_m15->Fill(maxcaljetpt,wt);
	if (maxl1pt>29.) mx_caljetpt_m30->Fill(maxcaljetpt,wt);
	if (maxl1pt>49.) mx_caljetpt_m50->Fill(maxcaljetpt,wt);
	if (maxl1pt>69.) mx_caljetpt_m70->Fill(maxcaljetpt,wt);
	if (maxl1pt>99.) mx_caljetpt_m100->Fill(maxcaljetpt,wt);
	if (maxl1pt>149.) mx_caljetpt_m150->Fill(maxcaljetpt,wt);
	if (maxl1pt>199.) mx_caljetpt_m200->Fill(maxcaljetpt,wt);

	// loop over all hlt trigger bits and fill histograms
	for (int i=0; i<nhlt; i++ ){
	  string hltname=hltnames[i];
	  if (hltTriggers[hltname]>0) h_cal[i]->Fill(maxcaljetpt,wt);
	}
      }
      
     // ************** Generated Jet Info ***********************

      double maxgenjetpt=0,maxgenjeteta=-999.;
      for (int i = 0; i < njetgen; ++i){
	//printf("j=%d value=%f\n",i,jetpt[i]);
	if (genjeteta[i]>mineta && genjeteta[i]< maxeta){
	  h_genjetpt->Fill(genjetpt[i],wt);
	  if (genjetpt[i] > maxgenjetpt) {
	    maxgenjetpt=genjetpt[i];
	    maxgenjeteta=genjeteta[i];
	  }
	}
      }
      if (maxgenjetpt>0.){

	if (abs(maxgenjeteta)<1.4){
	  h_l1met_2->Fill(l1met,wt);
	}
	  
	if (maxgenjetpt>45. && maxgenjetpt<55. ) h_l1->Fill(maxl1pt/maxgenjetpt,1.);
	h_l1vsGen->Fill(maxgenjetpt,maxl1pt,1.);

	mx_genjetpt_Untriggered->Fill(maxgenjetpt,wt);
	if (i_l1bit15>0) mx_genjetpt_l15->Fill(maxgenjetpt,wt);
	if (i_l1bit30>0) mx_genjetpt_l30->Fill(maxgenjetpt,wt);
	if (i_l1bit50>0) mx_genjetpt_l50->Fill(maxgenjetpt,wt);
	if (i_l1bit70>0) mx_genjetpt_l70->Fill(maxgenjetpt,wt);
	if (i_l1bit100>0) mx_genjetpt_l100->Fill(maxgenjetpt,wt);
	if (i_l1bit150>0) mx_genjetpt_l150->Fill(maxgenjetpt,wt);
	if (i_l1bit200>0) mx_genjetpt_l200->Fill(maxgenjetpt,wt);

	if (maxl1pt>14.) mx_genjetpt_m15->Fill(maxgenjetpt,wt);
	if (maxl1pt>29.) mx_genjetpt_m30->Fill(maxgenjetpt,wt);
	if (maxl1pt>49.) mx_genjetpt_m50->Fill(maxgenjetpt,wt);
	if (maxl1pt>69.) mx_genjetpt_m70->Fill(maxgenjetpt,wt);
	if (maxl1pt>99.) mx_genjetpt_m100->Fill(maxgenjetpt,wt);
	if (maxl1pt>149.) mx_genjetpt_m150->Fill(maxgenjetpt,wt);
	if (maxl1pt>199.) mx_genjetpt_m200->Fill(maxgenjetpt,wt);

	// reconstruct hlt triggers
	if (i_l1bit15>0 && maxcaljetpt>30)   mx_genjetpt_hlt30->Fill(maxgenjetpt,wt);
	if (i_l1bit30>0 && maxcaljetpt>60)   mx_genjetpt_hlt60->Fill(maxgenjetpt,wt);
	if (i_l1bit30>0 && maxcaljetpt>65)   mx_genjetpt_hlt65->Fill(maxgenjetpt,wt);
	if (i_l1bit30>0 && maxcaljetpt>70)   mx_genjetpt_hlt70->Fill(maxgenjetpt,wt);
	if (i_l1bit70>0 && maxcaljetpt>110)  mx_genjetpt_hlt110->Fill(maxgenjetpt,wt);
	if (i_l1bit100>0 && maxcaljetpt>150) mx_genjetpt_hlt150->Fill(maxgenjetpt,wt);
	if (i_l1bit150>0 && maxcaljetpt>200) mx_genjetpt_hlt200->Fill(maxgenjetpt,wt);
	if (i_l1bit200>0 && maxcaljetpt>250) mx_genjetpt_hlt250->Fill(maxgenjetpt,wt);

	if (i_l1bit100>0 && maxl1pt<99.9)
	  cout << " l1bit100 set. Maximum l1pt= " << maxl1pt << endl;

	// loop over all hlt trigger bits and fill histograms
	for (int i=0; i<nhlt; i++ ){
	  string hltname=hltnames[i];
	  if (hltTriggers[hltname]>0) h_gen[i]->Fill(maxgenjetpt,wt);
	}

      }

      // ***************** MET Triggers *****************************8

      met_Untriggered->Fill(metpt,wt);
      if (genmetpt < 5.) met_Untriggered_gencut->Fill(metpt,wt);
      summet_Untriggered->Fill(metsum,wt);
      sumht_Untriggered ->Fill(htsum,wt);
      htvspt->Fill(maxcaljetpt,htsum,wt);
      htvsl1ht->Fill(l1mettothad,htsum,wt);

      if (hltTriggers[hlt1met]>0) met_trg->Fill(metpt,wt);
      if (hltTriggers[hlt1metPre1]>0) met_trgPre1->Fill(metpt,wt);
      if (hltTriggers[hlt1metPre2]>0) met_trgPre2->Fill(metpt,wt);
      if (hltTriggers[hlt1metPre3]>0) met_trgPre3->Fill(metpt,wt);

      if (i_etm20>0 )met_etm20->Fill(metpt,wt);
      if (i_etm30>0 )met_etm30->Fill(metpt,wt);
      if (i_etm40>0 )met_etm40->Fill(metpt,wt);
      if (i_etm50>0 )met_etm50->Fill(metpt,wt);
      if (i_etm60>0 )met_etm60->Fill(metpt,wt);

      double maxjetpt = maxcaljetpt;
      jetplusmetRates(metpt, maxjetpt, i_etm30, wt);

    }// end loop over entries
  }// end loop over files


  if (plotIt){
    h_genjetpt->Draw();
    //met_trg->SetLineColor(kRed);
    //met_trg->Draw("same");
  }

  // Create histo file and book histograms
  cout << "\n Writing histograms to: " << _HistName << endl;
  TFile histofile(_HistName.c_str(),"RECREATE");
  //TFile histofile("/dev/null","RECREATE");

  histofile.cd();
    // save histograms
  //  Hlist.Write();
  Hlist->Write();
  //histofile.ls();

  histofile.Close();

  cout << " " << endl;

}
