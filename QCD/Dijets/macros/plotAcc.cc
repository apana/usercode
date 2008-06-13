#include "rootfuncs.h"

void plotAcc(){

  gROOT->Reset();
  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetFuncWidth(2);

  Double_t etamin=-3.0, etamax=3.0;
  etamin=-2.;etamax=2.;

  //gDirectory->ls();
  TH2F *htmp = (TH2F*)gDirectory->Get("h2");
  if ( htmp != 0 ) delete htmp;
  //gDirectory->ls();

  TH2F *h2 = new TH2F("h2","xxx ",100,etamin,etamax,100,etamin,etamax);
  h2->GetXaxis()->SetTitle("#eta_{1}");
  h2->GetYaxis()->SetTitle("#eta_{2}");

  h2->Draw();

  bool drawBox=true;

  if (drawBox){
    Double_t xb1,xb2,xb3,xb4;
    Double_t yb1,yb2,yb3,yb4;
    if (fabs(etamax-2.0)<0.01){
      //xb1=-1; xb2=-2; xb3=1; xb4=2;
      //yb1=-2; yb2=-1; yb3=2; yb4=1;
      xb1=1; xb2=-2; xb3=-1; xb4=2;
      yb1=-2; yb2=1; yb3=2; yb4=-1;
    }else if (fabs(etamax-3.0)<0.01){
      xb1=0; xb2=-3; xb3=0; xb4=3;
      yb1=-3; yb2=0; yb3=3; yb4=0;
    }else{
      cout << "Problem with eta limits" << endl;
      return;
    }

    Double_t x[4] = {xb1,xb2,xb3,xb4};
    Double_t y[4] = {yb1,yb2,yb3,yb4};

    TPolyLine *pline = new TPolyLine(4,x,y);
    Color_t boxColor=40;
    pline->SetFillColor(boxColor);
    pline->SetLineColor(boxColor);
    pline->SetLineWidth(1);
    pline->Draw("f");
    pline->Draw();

  }

  const int nmax = 11;

  for (int i=0; i != nmax; i++){
    ostringstream ch_i("");
    ch_i<< i;
    TString fname="f" + ch_i.str();
    TF1 *pfuncs = new TF1(fname,"pol1",etamin,etamax);
    pfuncs->SetParameter( 0, -5+i);
    pfuncs->SetParameter( 1, 1);
    pfuncs->SetLineStyle(1);
    if ( i != 5) pfuncs->Draw("same");
  }

  for (int i=0; i != nmax; i++){
    ostringstream ch_i("");
    ch_i << i;
    TString fname="n" + ch_i.str();
    TF1* nfuncs = new TF1(fname,"pol1",etamin,etamax);
    nfuncs->SetParameter( 0, -5+i);
    nfuncs->SetParameter( 1, -1);
    nfuncs->SetLineStyle(2);
    if ( i != 5) nfuncs->Draw("same");
  }

  TLatex l1,l2,l3,l4,l5,l6;
  Int_t align=21; Double_t tsize=0.04;
  l1.SetTextAlign(align);  l1.SetTextSize(tsize);
  l2.SetTextAlign(align);  l2.SetTextSize(tsize);
  l3.SetTextAlign(align);  l3.SetTextSize(tsize);
  l4.SetTextAlign(align);  l4.SetTextSize(tsize);

  float etab=1.5,etastar=1.5;
  ostringstream ch_bmin(""), ch_bmax("");
  ostringstream ch_smin(""), ch_smax("");
  Float_t xl1,xl2,xl3,xl4;
  Float_t yl1,yl2,yl3,yl4;
  if (fabs(etamax-2.0)<0.01){
    // etastar=0.5;
    etab=0.5;
    ch_bmax << etab;
    ch_bmin << -etab;
    ch_smax << etastar;
    ch_smin << -etastar;

    xl1=-0.62; yl1=-0.62;
    xl2=.37; yl2=.37;
    xl3=-1.35; yl3=1.35;
    xl4=1.65; yl4=-1.65;

  }else if (fabs(etamax-3.0)<0.01){
    etastar=1.5;
    ch_bmax << etab;
    ch_bmin << -etab;
    ch_smax << etastar;
    ch_smin << -etastar;

    xl1=-1.35; yl1=-1.35;
    xl2=1.3; yl2=1.3;
    xl3=-1.25; yl3=1.25;
    xl4=1.45; yl4=-1.45;
  }

  TString label1="#eta_{boost}= " + ch_bmin.str();
  TString label2="#eta_{boost}= " + ch_bmax.str();
  TString label3="#eta^{*}= " + ch_smin.str();
  TString label4="#eta^{*}= " + ch_smax.str();
  l1.SetTextAngle(-45); l1.DrawLatex(xl1, yl1, label1);
  l2.SetTextAngle(-45); l2.DrawLatex(xl2, yl2, label2);
  l3.SetTextAngle(45); l3.DrawLatex(xl3, yl3, label3);
  l4.SetTextAngle(45); l4.DrawLatex(xl4, yl4, label4);

  ostringstream ch_eta("");
  ch_eta << etamax;
  TString psfile="DiJet_Acceptance_etamax" + ch_eta.str() + ".eps";
  cout << "Writing postscript file: " << psfile << endl;
  c1->Print(psfile);

  return;

}
