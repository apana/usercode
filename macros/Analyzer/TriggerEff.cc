#include "rootfuncs.h"

void SetupCanvas(){

  TCanvas *c1 = new TCanvas("c1","Root Canvas",450,20,840,675); // My Default Canvas
  Int_t wtopx,wtopy; UInt_t ww, wh;
  c1->GetCanvasPar(wtopx,wtopy,ww,wh); // Gets Canvas Parameters
  //cout << wtopx << " " << wtopy << " " << ww << " " << wh << endl;
  //TCanvas *c2 = new TCanvas("c2","Root Canvas 2");
  //TCanvas *c3 = new TCanvas("c3","Root Canvas 3",2);

  pad1 = new TPad("pad1","This is pad1",0.0,0.5,0.5,1.00);
  pad2 = new TPad("pad2","This is pad1",0.5,0.5,1.00,1.00);
  pad3 = new TPad("pad3","This is pad3",0.0,0.0,0.5,0.5);
  pad4 = new TPad("pad4","This is pad4",0.5,0.0,1.00,0.5);

  pad1->UseCurrentStyle();
  pad2->UseCurrentStyle();
  pad3->UseCurrentStyle();
  pad4->UseCurrentStyle();

  bool useColor=false;
  useColor=true;
  if (!useColor){
    color=(TColor*)(gROOT->GetListOfColors()->At(kMagenta)); color->SetRGB(0,0,0);
  }
  // Set the right and left margins for all pads
  Double_t p1_lm=0.100, p3_lm=p1_lm;   // pad1 and 3 left margins
  Double_t p1_rm=0.05, p3_rm=p1_rm;   // pad1 and 3 right margins

  Double_t p2_lm=p1_rm, p4_lm=p1_rm;
  Double_t p2_rm=p1_lm, p4_rm=p1_lm;

  pad1->SetLeftMargin(p1_lm); pad1->SetRightMargin(p1_rm);
  pad2->SetLeftMargin(p2_lm); pad2->SetRightMargin(p2_rm);
  pad3->SetLeftMargin(p3_lm); pad3->SetRightMargin(p3_rm);
  pad4->SetLeftMargin(p4_lm); pad4->SetRightMargin(p4_rm);

  // Set the top and bottom margins for all pads
  Double_t p1_tm=0.125, p2_tm=p1_tm;   // pad1 and 2 top margins
  Double_t p1_bm=0.09, p2_bm=p1_bm;   // pad1 and 2 bottom margins

  Double_t p3_tm=p1_bm, p4_tm=p2_bm;
  Double_t p3_bm=p1_tm, p4_bm=p2_tm;

  pad1->SetTopMargin(p1_tm); pad1->SetBottomMargin(p1_bm);
  pad2->SetTopMargin(p2_tm); pad2->SetBottomMargin(p2_bm);
  pad3->SetTopMargin(p3_tm); pad3->SetBottomMargin(p3_bm);
  pad4->SetTopMargin(p4_tm); pad4->SetBottomMargin(p4_bm);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();

}

void TriggerEff(){

  gROOT->Reset();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetOptLogy(0);

  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  
  gStyle->SetLabelSize(0.04,"y");
  gStyle->SetLabelSize(0.045,"x");

  gStyle->SetHistFillColor(kYellow);
  gROOT->ForceStyle();

  bool plotL1=false;
  bool DrawL1Cut=true;
  if (!plotL1)
    DrawL1Cut=false;

  TString rootname = "./hlt_hists_RelVal_160_unwt.root";

  TFile *rootfile=OpenRootFile(rootname); if (!rootfile) return;
  rootfile->GetListOfKeys()->Print();


  //TString WhichJets="mx_caljetpt_"; // Calo or Gen
  TString WhichJets="mx_genjetpt_"; // Calo or Gen
  
  TString xaxis="Leading HLT Corrected Jet p_{T} (GeV)";
  if (WhichJets=="mx_genjetpt_"){
    xaxis="Leading HLT Generated Jet p_{T} (GeV)";
  }
  //  TString hname="Genpt_Leading_Untriggered_combined";
  TString hname=WhichJets + "Untriggered";
  TH1F *hden=GetHist(rootfile,hname); if (!hden) return;

  TH1 *hnum1,*hnum2,*hnum3,*hnum4;

  if (plotL1){
    hname=WhichJets + "l150";     hnum1=GetHist(rootfile,hname); if (!hnum1) return;
    hname=WhichJets + "l100";  hnum2=GetHist(rootfile,hname); if (!hnum2) return;
    hname=WhichJets + "l70";  hnum3=GetHist(rootfile,hname); if (!hnum3) return;
    hname=WhichJets + "l30";  hnum4=GetHist(rootfile,hname); if (!hnum4) return;
  }else{
    hname=WhichJets + "HLT1jet";     hnum1=GetHist(rootfile,hname); if (!hnum1) return;
    //hname=WhichJets + "HLT1jetPE1";  hnum2=GetHist(rootfile,hname); if (!hnum2) return;
    cout << "Using open-HLT1jetPE1 trigger since I blew the prescale factor" << endl;
    hname=WhichJets + "hlt150";  hnum2=GetHist(rootfile,hname); if (!hnum2) return;
    hname=WhichJets + "HLT1jetPE3";  hnum3=GetHist(rootfile,hname); if (!hnum3) return;
    hname=WhichJets + "HLT1jetPE5";  hnum4=GetHist(rootfile,hname); if (!hnum4) return;
    //hname=WhichJets + "hlt60";  hnum4=GetHist(rootfile,hname); if (!hnum4) return;
  }

  Double_t rebin=2;
  hnum1->Rebin(rebin);
  hnum2->Rebin(rebin);
  hnum3->Rebin(rebin);
  hnum4->Rebin(rebin);
  hden->Rebin(rebin);

  //hden->Sumw2();
  //hnum1->Sumw2();
  //hnum2->Sumw2();
  //hnum3->Sumw2();
  //hnum4->Sumw2();


  bool useColor=true;
  if (!useColor){
    color=(TColor*)(gROOT->GetListOfColors()->At(kMagenta)); color->SetRGB(0,0,0);
    color=(TColor*)(gROOT->GetListOfColors()->At(kBlue)); color->SetRGB(0,0,0);
  }

  //Double_t ymin=-.09,ymax=1.09;
  Double_t ymin=0.0,ymax=1.09;
  Double_t eff1minpt=0.,eff1maxpt=649;  
  Double_t eff2minpt=0.,eff2maxpt=399;
  Double_t eff3minpt=0.,eff3maxpt=249;
  Double_t eff4minpt=0.,eff4maxpt=149;
  if (! plotL1){
    //eff1maxpt=749; eff2maxpt=499; eff3maxpt=349; eff4maxpt=249;
    //eff1maxpt=449; eff2maxpt=399; eff3maxpt=349; eff4maxpt=249;
    eff1maxpt=399; eff2maxpt=349; eff3maxpt=299; eff4maxpt=249;
  }

  //  Double_t hltcut1=250, hltcut2=150, hltcut3=110, hltcut4=60;
  Double_t hltcut1=200, hltcut2=150, hltcut3=110, hltcut4=60;
  Double_t l1cut1=150, l1cut2=100, l1cut3=70, l1cut4=30;
  Double_t ylmin=0.75, ylmax=ymax, yeff=0.95;
  if (! plotL1){
    yeff=0.5; ylmin=yeff-0.25; ylmax=yeff; 
  }
  Double_t ylminl1=0.25, ylmaxl1=0.5;
  Double_t yarr1=0.5, yarr2=yarr1+0.2;
  Double_t ytxt=ylmin-0.075;

  Double_t ytxtl1=ylminl1-0.075;


  //  TString hltlbl1="HLT > 250",hltlbl2="HLT > 150",hltlbl3="HLT > 110",hltlbl4="HLT > 60";
  TString hltlbl1="HLT > 200",hltlbl2="HLT > 150",hltlbl3="HLT > 110",hltlbl4="HLT > 60";
  TString l1lbl1= "L1 > 150" ,l1lbl2= "L1 > 100" ,l1lbl3= "L1 > 70"  ,l1lbl4 ="L1 > 30";

  //hnum2->Scale(10.);
  TH1F *eff1 = (TH1F*)hnum1->Clone(); eff1->SetName("eff1");
  TH1F *eff2 = (TH1F*)hnum1->Clone(); eff2->SetName("eff2");
  TH1F *eff3 = (TH1F*)hnum1->Clone(); eff3->SetName("eff3");
  TH1F *eff4 = (TH1F*)hnum1->Clone(); eff4->SetName("eff4");
  eff1->Divide(hnum1,hden,1.,1.,"B");
  eff2->Divide(hnum2,hden,1.,1.,"B");
  eff3->Divide(hnum3,hden,1.,1.,"B");
  eff4->Divide(hnum4,hden,1.,1.,"B");
  
  eff1->GetXaxis()->SetRangeUser(eff1minpt,eff1maxpt);
  eff2->GetXaxis()->SetRangeUser(eff2minpt,eff2maxpt);
  eff3->GetXaxis()->SetRangeUser(eff3minpt,eff3maxpt);
  eff4->GetXaxis()->SetRangeUser(eff4minpt,eff4maxpt);

  eff1->GetYaxis()->SetRangeUser(ymin,ymax);
  eff2->GetYaxis()->SetRangeUser(ymin,ymax);
  eff3->GetYaxis()->SetRangeUser(ymin,ymax);
  eff4->GetYaxis()->SetRangeUser(ymin,ymax);

  SetupCanvas();

  TH1F *hists[]={eff4,eff3,eff2,eff1};
  TPad *pads[]={pad1,pad2,pad3,pad4};
  TLine *l1[]={0,0,0,0};
  TLine *l2[]={0,0,0,0};
  TLine *l3[]={0,0,0,0};
  TLine *l4[]={0,0,0,0};

  TMarker *m1[]={0,0,0,0};

  TLine *arrow[]={0,0,0,0};
  Double_t xlineh[]={hltcut4,hltcut3,hltcut2,hltcut1};
  Double_t xlinel[]={l1cut4,l1cut3,l1cut2,l1cut1};
  Double_t ptmax[]={eff4maxpt,eff3maxpt,eff2maxpt,eff1maxpt};

  TPaveText *lbls[]={0,0,0,0};
  TString hltlbls[]={hltlbl4,hltlbl3,hltlbl2,hltlbl1};
  TString l1lbls[]={l1lbl4,l1lbl3,l1lbl2,l1lbl1};
  int nhist=4;

  Style_t ltype=3; Width_t lwid=2;
  Double_t xp1=0.6, xp2=xp1+0.2, yp1=0.35, yp2=yp1+0.2;
  //if (!plotL1) {yp1=0.55; yp2=yp1+0.1;}
  if (!plotL1) {yp1=0.55; yp2=yp1+0.15;}
  Float_t arrsiz=0.02;

  TLatex *t = new TLatex(); Double_t xtxt;

  for (int i=0;i<nhist;i++){
    pads[i]->cd();
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.05);
    if (i == 0 || i == 2)
      hists[i]->GetYaxis()->SetTitleOffset(0.8);
    else
      hists[i]->GetYaxis()->SetTitleOffset(1.8);

    hists[i]->Draw("hist");
    hists[i]->Draw("e,same");

    hists[i]->GetYaxis()->SetTitle("Efficiency");
    hists[i]->GetXaxis()->SetTitle(xaxis);

    l1[i] = new TLine(xlineh[i],ylmin,xlineh[i],ylmax);
    l1[i]->SetLineStyle(ltype); l1[i]->SetLineWidth(lwid); ; l1[i]->SetLineColor(kMagenta);
    l1[i]->Draw();

    if (plotL1) {
      l2[i] = new TLine(0,yeff,ptmax[i],yeff);
      l2[i]->SetLineStyle(ltype); l2[i]->SetLineWidth(lwid); ; l2[i]->SetLineColor(kMagenta);
      l2[i]->Draw();
    }

    if (!plotL1){
      m1[i] = new TMarker(xlineh[i],yeff,20);
      m1[i]->SetMarkerColor(kMagenta);
      m1[i]->Draw();
    }

    if (DrawL1Cut){
      l3[i] = new TLine(xlinel[i],ylminl1,xlinel[i],ylmaxl1);
      l3[i]->SetLineStyle(ltype); l3[i]->SetLineWidth(lwid); ; l3[i]->SetLineColor(kMagenta);
      l3[i]->Draw();

      l4[i] = new TLine(0,0.5,ptmax[i],0.5);
      l4[i]->SetLineStyle(ltype); l4[i]->SetLineWidth(lwid); ; l4[i]->SetLineColor(kMagenta);
      //l4[i]->Draw();
      m1[i] = new TMarker(xlinel[i],0.5,20);
      m1[i]->SetMarkerColor(kMagenta);
      m1[i]->Draw();

      xtxt=xlinel[i]; 
      t->SetTextColor(kMagenta);
      t->DrawLatex(xtxt,ytxtl1,"L1 Cut");
      
    }
    //arrow[i] = new TArrow(xline[i],yarr1,xline[i],yarr2,arrsiz,"->-");
    //arrow[i]->SetLineStyle(1); arrow[i]->SetLineWidth(1);
    //arrow[i]->Draw();

    xtxt=xlineh[i]; 
    t->SetTextColor(kMagenta);
    t->DrawLatex(xtxt,ytxt,"HLT Cut");

    
    lbls[i] = new TPaveText(xp1,yp1,xp2,yp2,"NDC");
    //    if (plotL1)
    TText *t1=lbls[i]->AddText(l1lbls[i]);
    TText *t1=lbls[i]->AddText(hltlbls[i]);
    lbls[i]->Draw();
  }

  c1->Update();

  bool saveIt=false;
  saveIt=true;
  if (saveIt){
    TString psfile="hltturnOns" ;
    if (plotL1) psfile="l1turnOns" ;
    psfile=psfile + ".ps";
    cout << "Writing postscript file: " << psfile << endl;
    c1->Print(psfile);
  }

  return;  
}
