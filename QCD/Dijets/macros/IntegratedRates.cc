#include "rootfuncs.h"
#include "rateIntegrator.h"

// Integrated Luminosity in pb-1
//const float ILumi=10.; 
//const float ILumi=0.0386; // Jet 20 Lumi
//const float ILumi=0.1917; // Jet 30 Lumi
//const float ILumi=1.6145; // Jet 50 Lumi
//const float ILumi=8.9317; // Jet 80 Lumi
//const float ILumi=44.8552;  // Jet 110 Lumi
const float ILumi=243.7150; // Jet 150 Lumi

const float conv=1e-36; // conversion from pb to cm^2

void myrateIntegrator(TH1* hin, TH1* hout, float Lumi,bool reScale=true){

  Float_t dy=1.;  // take out average over rapidity

  Int_t nbins=hin->GetNbinsX();
  //cout << "Number of bins: " << nbins << endl;

//   for (Int_t ibin=0; ibin<nbins; ++ibin){
//     Float_t cont=h->GetBinContent(ibin+1);
//     Float_t err =h->GetBinError(ibin+1);
//     Float_t bw  =h->GetBinWidth(ibin+1);
//     bw=1.;

//     Float_t fact=bw*dy;

//     //cout << ibin << " " << cont << " +- " << err << " -- BinWidth: " << bw << endl;
//     cont=cont*fact;
//     err = err*fact;
    
//     h->SetBinContent(ibin+1,cont);
//     h->SetBinError(ibin+1,err);
//   }

  // now the integrated rate
  // error on rate is calculated assuming error on integral is dominated by current bin error
  for (Int_t ibin=0; ibin<nbins; ++ibin){
    Int_t ibin1=ibin+1;
    Int_t ibin2=nbins;

    Float_t cont =hin->GetBinContent(ibin+1);
    Float_t err  =hin->GetBinError(ibin+1);
    
    Float_t ferr=0.;
    if (cont>0.)ferr=err/cont;

    //Float_t fint=h->Integral(ibin1,ibin2,"width");
    Float_t fint=hin->Integral(ibin1,ibin2,"");

    //cout << ibin << "  -- Integral: " << fint << endl;
    cont=fint;
    err = ferr*fint;
    
    hout->SetBinContent(ibin+1,cont);
    hout->SetBinError(ibin+1,err);
  }
  //cout << "  -- Integral: " << copy_h1->Integral() << endl;


  Float_t fact=1.;
  if (reScale){
    const float conv=1e-36; // conversion from pb to cm^2
    //const  float conv=1e-27; // conversion from mb to cm^2
    
    cout << "Calculating rate for Luminosity: " << Lumi << endl;
    // now multiply by luminosity to get total rate
    Float_t fact=Lumi*conv;
  }

  cout << "%myRateIntegrator - Conversion factor: " << fact << endl;
  hout->Scale(fact);
}

void IntegratedRates(){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetOptLogy(1);

  gStyle->SetHistLineWidth(2);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.02);

  gROOT->ForceStyle();

  bool reScale=false;
  reScale=true;

  TString TrigStream="Jet150";
  TString WhichJets="Gen"; 
  WhichJets="CaloMCCor"; 

  TString vers="xxx";
  TString rootname = "../nhsts/DiJets_" + TrigStream + "_Scone7_" +WhichJets + "_eta3_wght0.root"; vers="207";

  double wt=1;
  if (reScale){
    if (TrigStream == "Jet20"){
      double Nevt=3926600;
      double xs= 101600000;    
    } else if (TrigStream == "Jet30"){
      double Nevt=4131600;
      double xs=21550000;
    } else if (TrigStream == "Jet50"){
      double Nevt=4010400;
      double xs=2484000;
    } else if (TrigStream == "Jet80"){
      double Nevt=2891200;
      double xs=323700;
    } else if (TrigStream == "Jet110"){
      double Nevt=3980000;
      double xs=88730;
    } else if (TrigStream == "Jet150"){
      double Nevt=4171600;
      double xs=17120;
    } else {
      cout << "Problems! " << endl;
      return;
    }
    wt=xs/Nevt;
  }
  TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;
  //rootfile->GetListOfKeys()->Print();


  Double_t ptmin=10,ptmax=499;


  TString xaxisLabel="Leading GenJet p_{T} (GeV/c)";
  if ( WhichJets == "CaloMCCor") xaxisLabel="Leading CaloJet p_{T} (GeV/c)";

  const int nhists=1;
  TH1 *hists[nhists];
  TH1F *clones[nhists];
  TString trignames[nhists]={"JetPtMax"};

  float ptcut[nhists]={250.};

  float prescl[nhists];
  prescl[0]=1.;

  Int_t hcols[nhists]={kBlue};
  Int_t ltyps[nhists]={   1};

  TString GenRate[nhists];
  TString HLTRate[nhists];

  TString chopt;

  int n=nhists;
  for (int ihist=0;ihist<n;ihist++){

    TString hname=trignames[ihist];
    hists[ihist]=GetHist(rootfile,hname); if (!hists[ihist]) return;

    Int_t rebin=1;
    cout << "Wt= " << wt << endl;
    hists[ihist]->Rebin(rebin); hists[ihist]->Scale(1./rebin);
    hists[ihist]->Scale(wt/prescl[ihist]);

    hname=hname+ "_clone";
    clones[ihist] = (TH1F*)hists[ihist]->Clone(); clones[ihist]->SetName(hname);

    clones[ihist]->Reset();
    clones[ihist]->SetLineColor(hcols[ihist]);
    clones[ihist]->SetLineStyle(ltyps[ihist]);
    float Lumi=ILumi/conv; // rateIntegrator expects lumi in cm-2
    //float Lumi=1.e32;
    cout << Lumi << endl;

    myrateIntegrator(hists[ihist],clones[ihist],Lumi,reScale);

    //chopt="hist,same";
    chopt="e,same";
    if (ihist == 0){
      //chopt="hist";
      clones[ihist]->GetYaxis()->SetTitle("Integrated Events");
      clones[ihist]->GetXaxis()->SetTitle(xaxisLabel);
      clones[ihist]->GetXaxis()->SetRangeUser(ptmin,ptmax);

      //clones[ihist]->Draw("hist,e");
      clones[ihist]->Draw("e");

      Int_t ibin= clones[ihist]->GetXaxis()->FindBin(ptcut[ihist]);
      Double_t pre=prescl[ihist];
      Double_t trgRate=clones[ihist]->GetBinContent(1);
      Double_t untRate=clones[0]->GetBinContent(ibin);

      pre=prescl[ihist];
      cout << "\n Rate for pt threshold Unt:" << ptcut[ihist] << " : " << untRate << "/" << prescl[ihist] << " =" << untRate/prescl[ihist] << endl;

    }else{
      Int_t ibin= clones[ihist]->GetXaxis()->FindBin(ptcut[ihist]);
      Double_t pre=prescl[ihist];
      Double_t trgRate=clones[ihist]->GetBinContent(1);
      Double_t untRate=clones[0]->GetBinContent(ibin);

      pre=prescl[ihist];
      cout << "\n Rate for pt threshold Unt:" << ptcut[ihist] << " : " << untRate << "/" << prescl[ihist] << " =" << untRate/prescl[ihist] << endl;
      cout << " Rate for pt threshold Trg: " << ptcut[ihist] << " : " << trgRate << "/" << prescl[ihist] << " =" << trgRate/prescl[ihist]<< endl;
      //printf("Total CPU time:\t%2.1f  Hz\n",untRate);
      ostringstream ch_z1,ch_z2;
      if (untRate < 10. ) ch_z1.flags ( ios_base::showpoint );
      ch_z1.precision(2);
      //ch_z1 << untRate/pre ;
      ch_z1 << untRate ;

      if (trgRate < 10. ) ch_z2.flags ( ios_base::showpoint );
      ch_z2.precision(2);
      //ch_z2 << trgRate/pre ;
      ch_z2 << trgRate ;
      
      TString Untoutput=ch_z1.str();
      TString Trgoutput=ch_z2.str();
      cout << " Untriggered: " << Untoutput << endl;
      cout << " Triggered: " << Trgoutput << endl;

      GenRate[ihist]=Untoutput;
      HLTRate[ihist]=Trgoutput;
      
    }

    clones[ihist]->Draw(chopt);
 
  }

  //clones[0]->SetFillColor(0);
  //clones[0]->Draw(chopt);
  cout << "\n  %Luminosity: " << ILumi << " pb-1" << endl;
  if (reScale)
    cout << " Weighted histograms by XS/NEVT " << endl;
  else
    cout << " No XS weighting done " << endl;

  return;

  Int_t txtsize=17;
  Int_t txtfnt=63;
  Int_t txtalign=13;


  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextAlign(txtalign);
  t1->SetTextFont(txtfnt);
  t1->SetTextSizePixels(txtsize);

  /*
  Double_t xtxt0=0.54, ytxt0=0.25;
  Double_t xtxt1=0.39, ytxt1=0.335;
  Double_t xtxt2=0.19, ytxt2=0.43;
  Double_t xtxt3=0.275, ytxt3=0.575;
  */

  Double_t xtxt0=0.58, ytxt0=0.28;
  Double_t xtxt1=0.4, ytxt1=0.39;
  Double_t xtxt2=0.19, ytxt2=0.49;
  Double_t xtxt3=0.31, ytxt3=0.65;
  Double_t xtxt4=0.20 , ytxt4=0.85;

  bool showRates=true;
  //showRates=false;
  if (showRates){
    bool showUntRate=false;
    //showUntRate=true;
    Float_t xtxt=xtxt0, ytxt=ytxt0;  
    if (showUntRate){
      t1->SetTextColor(kRed);
      t1->DrawLatex(xtxt,ytxt,GenRate[1] + "Hz");
      ytxt=ytxt-0.03;
    }
    t1->SetTextColor(kBlue);
    t1->DrawLatex(xtxt,ytxt,HLTRate[1] + "Hz");
    
    xtxt=xtxt1, ytxt=ytxt1;  
    if (showUntRate){
      t1->SetTextColor(kRed);
      t1->DrawLatex(xtxt,ytxt,GenRate[2] + "Hz");
      ytxt=ytxt-0.03;
    }
    t1->SetTextColor(kBlue);
    t1->DrawLatex(xtxt,ytxt,HLTRate[2] + "Hz");
    
    xtxt=xtxt2, ytxt=ytxt2;  
    if (showUntRate){
      t1->SetTextColor(kRed);
      t1->DrawLatex(xtxt,ytxt,GenRate[3] + "Hz");
      ytxt=ytxt-0.03;
    }
    t1->SetTextColor(kBlue);
    t1->DrawLatex(xtxt,ytxt,HLTRate[3] + "Hz");
    
    xtxt=xtxt3, ytxt=ytxt3;  
    if (showUntRate){
      t1->SetTextColor(kRed);
      t1->DrawLatex(xtxt,ytxt,GenRate[4] + "Hz");
      ytxt=ytxt-0.03;
    }
    t1->SetTextColor(kBlue);
    t1->DrawLatex(xtxt,ytxt,HLTRate[4] + "Hz");

    xtxt=xtxt4, ytxt=ytxt4;  
    if (showUntRate){
      t1->SetTextColor(kRed);
      t1->DrawLatex(xtxt,ytxt,GenRate[5] + "Hz");
      ytxt=ytxt-0.03;
    }
    t1->SetTextColor(kBlue);
    t1->DrawLatex(xtxt,ytxt,HLTRate[5] + "Hz");

    //  }
  }
  //printf("Total CPU time:\t%2.2fu  %2.2fs\t%2.2f (sec)\n",0., 0., 0.);

  bool saveIt=false;
  //saveIt=true;
  if (saveIt){
    TString psfile="IntegratedJetRates";
    if (useColor) psfile= psfile + "_color";
    psfile=psfile + vers + ".eps";

    cout << "Writing postscript file: " << psfile << endl;
    c1->Print(psfile);
  }
}
