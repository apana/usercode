#include "rootfuncs.h"

double invMass(double *xx, double *pars);

double invMass(double *xx, double *pars)
{

  double x=xx[0];
  double et=pars[0];

  double value=2.*pow(et,2)*(cosh(log(x))+1);
  double fitval=value;

  return sqrt(fitval);
}

void DrawBox(double *x, double* y, Color_t boxColor=40){

  //Double_t x[4] = {xx[],.-3.,.0,3.,};
  //Double_t y[4] = {-3, 0.,  3,0.,};
  //TPolyLine *pline = new TPolyLine(4,x,y);
  TBox *box1 = new TBox(x[0],y[0],x[1],y[1]);
  TBox *box2 = new TBox(x[0],y[0],x[1],y[1]);

  box1->SetFillColor(boxColor);
  box1->SetFillStyle(3002);
  box1->SetLineColor(boxColor);
  box1->SetLineWidth(1);
  box1->Draw();

  box2->SetFillStyle(0);
  box2->SetLineColor(boxColor);
  box2->Draw();

}

void DrawMassLimits(double x, double y, const char* text="xxx"){

  TLatex *l = new TLatex();
  Int_t align=12; Double_t tsize=0.04;
  l->SetTextAlign(align);  l->SetTextSize(tsize), l->SetTextAngle(0);
  l->DrawLatex(x, y, text);

}

void plotMvsChi_CMS(){

  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetFuncWidth(2);
  gStyle->SetTitleOffset(1.35,"y");
 
  TF1 *f1 = new TF1("f1",invMass,1.,30.,1);
  TF1 *f2 = new TF1("f2",invMass,1.,30.,1);
  TF1 *f3 = new TF1("f3",invMass,1.,30.,1);
  TF1 *f4 = new TF1("f4",invMass,1.,30.,1);
  TF1 *f5 = new TF1("f5",invMass,1.,30.,1);
  TF1 *f6 = new TF1("f6",invMass,1.,30.,1);
  
  Color_t col1=kRed, col2=kBlue, col3=kDarkGreen, col4=38, col5=48, col6=58;

  TLatex l1,l2,l3,l4,l5,l6;
  Int_t align=12; Double_t tsize=0.04;
  l1.SetTextAlign(align);  l1.SetTextSize(tsize);
  l2.SetTextAlign(align);  l2.SetTextSize(tsize);
  l3.SetTextAlign(align);  l3.SetTextSize(tsize);
  l4.SetTextAlign(align);  l4.SetTextSize(tsize);
  l5.SetTextAlign(align);  l5.SetTextSize(tsize);
  l6.SetTextAlign(align);  l6.SetTextSize(tsize);

  int nEt=6;

  //double et1=375.,et2=310.,et3=246., et4=190., et5=125, et6=60;  // CMS org
  double et1=400.,et2=310.,et3=250., et4=190., et5=125, et6=80; 

  ostringstream ch_et1(""),ch_et2(""),ch_et3(""),ch_et4(""),ch_et5(""),ch_et6("");
  ch_et1 << et1; ch_et2 << et2; ch_et3 << et3; ch_et4 << et4; ch_et5 << et5; ch_et6 << et6;
  TString ch_label1="E_{T}=" + ch_et1.str() + " GeV";
  TString ch_label2="E_{T}=" + ch_et2.str() + " GeV";
  TString ch_label3="E_{T}=" + ch_et3.str() + " GeV";
  TString ch_label4="E_{T}=" + ch_et4.str() + " GeV";
  TString ch_label5="E_{T}=" + ch_et5.str() + " GeV";
  TString ch_label6="E_{T}=" + ch_et6.str() + " GeV";

  double xx,yy;

  double maxmass=2350;
  h2 = new TH2F("h2"," ",100,1.,30,100,50,maxmass);
  h2->Draw();

  f1->SetParameter(0,et1);
  f2->SetParameter(0,et2);
  f3->SetParameter(0,et3);
  f4->SetParameter(0,et4);
  f5->SetParameter(0,et5);
  f6->SetParameter(0,et6);

  h2->GetXaxis()->SetTitle("#chi");
  h2->GetYaxis()->SetTitle("Invariant Mass (GeV)");
  //f1->GetYaxis()->SetRangeUser(60,1100);

  f1->SetLineColor(col1);
  f2->SetLineColor(col2);
  f3->SetLineColor(col3);
  f4->SetLineColor(col4);
  f5->SetLineColor(col5);
  f6->SetLineColor(col6);

  Int_t lsty=5;
  f1->SetLineStyle(lsty);
  f2->SetLineStyle(lsty);
  f3->SetLineStyle(lsty);
  f4->SetLineStyle(lsty);
  f5->SetLineStyle(lsty);
  f6->SetLineStyle(lsty);

  //Double_t x[4] = {  1,   1,  11 , 11};
  //Double_t y[4] = {635,1100,1100, 635};

  Double_t dy=20;
  Double_t dy1=dy+30,dy2=dy+50,dy3=dy,dy4=dy,dy5=dy+20,dy6=dy+20;

  Double_t x[2] = {  1,   20};
  Double_t y[2] = {1880,maxmass};
  DrawBox(x,y, col1);
  l1.SetTextAngle(27); l1.DrawLatex(22, (y[1]+y[0])/2-dy1, ch_label1);
  xx=2. ; yy=(y[1]+y[0])/2;  DrawMassLimits(xx,yy,"M > 1880 GeV");

  Double_t x[2] = {  1,   20};
  Double_t y[2] = {1460, 1880};
  DrawBox(x,y, col2);
  l2.SetTextAngle(26); l2.DrawLatex(22, (y[1]+y[0])/2-dy2, ch_label2);
  xx=2. ; yy=(y[1]+y[0])/2;  DrawMassLimits(xx,yy,"1460 < M < 1880 GeV");

  Double_t x[2] = {  1,   20};
  Double_t y[2] = {1160, 1460};
  DrawBox(x,y, col3);
  l3.SetTextAngle(20); l3.DrawLatex(22, (y[1]+y[0])/2-dy3, ch_label3);
  xx=2. ; yy=(y[1]+y[0])/2;  DrawMassLimits(xx,yy,"1160 < M < 1460 GeV");

  Double_t x[2] = {  1,  20};
  Double_t y[2] = {900, 1160};
  DrawBox(x,y, col4);
  l4.SetTextAngle(15); l4.DrawLatex(22, (y[1]+y[0])/2-dy4, ch_label4);
  xx=2. ; yy=(y[1]+y[0])/2.;  DrawMassLimits(xx,yy,"900 < M < 1160 GeV");

  Double_t x[2] = {  1,  20};
  Double_t y[2] = {600, 900};
  DrawBox(x,y, col5);
  l5.SetTextAngle(10); l5.DrawLatex(22, (y[1]+y[0])/2-dy5, ch_label5);
  xx=2. ; yy=(y[1]+y[0])/2;  DrawMassLimits(xx,yy,"600 < M < 900 GeV");

  Double_t x[2] = {  1,  20};
  Double_t y[2] = {380, 600};
  DrawBox(x,y, col6);
  l6.SetTextAngle(5); l6.DrawLatex(22, (y[1]+y[0])/2-dy6, ch_label6);
  xx=2. ; yy=(y[1]+y[0])/2;  DrawMassLimits(xx,yy,"380 < M < 600 GeV");

  f1->Draw("same");
  f2->Draw("same");
  f3->Draw("same");
  f4->Draw("same");
  f5->Draw("same");
  f6->Draw("same");
  gPad->RedrawAxis();

  bool saveIt=false;
  if (saveIt){
    TString psfile="MassvsChi.eps";
    cout << "Writing postscript file: " << psfile << endl;
    c1->Print(psfile);

    TString pdffile="MassvsChi.pdf";
    cout << "Writing pdf file: " << pdffile << endl;
    c1->Print(pdffile);
  }

}
