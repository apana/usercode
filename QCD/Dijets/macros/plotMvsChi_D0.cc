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
void plotMvsChi_D0(){

  mySetup(0); // 0 -- linear plots, 1 -- log plots  
  gStyle->SetFuncWidth(2);
  TF1 *f1 = new TF1("f1",invMass,1.,30.,1);
  TF1 *f2 = new TF1("f2",invMass,1.,30.,1);
  TF1 *f3 = new TF1("f3",invMass,1.,30.,1);
  TF1 *f4 = new TF1("f4",invMass,1.,30.,1);
  
  Color_t col1=kRed, col2=kBlue, col3=kDarkGreen, col4=38;

  bool prl=true;

  int nEt=4;
  double et1=175.,et2=120.,et3=55., et4=0.;  // D0
  if (prl){
    nEt=4;
    et1=175; et2=120.; et3=90 ; et4=55;
  }

  //double et1=200.,et2=80.,et3=50.;  // D0

  h2 = new TH2F("h2"," ",100,1.,30,100,60,1100);
  h2->Draw();

  f1->SetParameter(0,et1);
  f2->SetParameter(0,et2);
  f3->SetParameter(0,et3);
  f4->SetParameter(0,et4);

  h2->GetXaxis()->SetTitle("#chi");
  h2->GetYaxis()->SetTitle("Invariant Mass (GeV)");
  //f1->GetYaxis()->SetRangeUser(60,1100);

  f1->SetLineColor(col1);
  f2->SetLineColor(col2);
  f3->SetLineColor(col3);
  f4->SetLineColor(col4);

  Int_t lsty=5;
  f1->SetLineStyle(lsty);
  f2->SetLineStyle(lsty);
  f3->SetLineStyle(lsty);
  f4->SetLineStyle(lsty);

  //Double_t x[4] = {  1,   1,  11 , 11};
  //Double_t y[4] = {635,1100,1100, 635};
  Double_t x[4] = {  1,   11};
  Double_t y[4] = {635,1100};
  DrawBox(x,y, col1);

  Double_t x[4] = {  1,   13};
  Double_t y[4] = {475, 635};
  DrawBox(x,y, col2);

  Double_t x[4] = {  1,   20};
  Double_t y[4] = {425, 475};
  DrawBox(x,y, col3);

  Double_t x[4] = {  1,  20};
  Double_t y[4] = {260, 425};
  DrawBox(x,y, col4);

  TLatex l1,l2,l3,l4;
  Int_t align=12; Double_t tsize=0.04;
  l1.SetTextAlign(align);  l1.SetTextSize(tsize), l1.SetTextAngle(25);
  l1.DrawLatex(22, 800, "E_{T}=175 GeV");

  l2.SetTextAlign(align);  l2.SetTextSize(tsize), l2.SetTextAngle(15);
  l2.DrawLatex(22, 550, "E_{T}=120 GeV");

  l3.SetTextAlign(align);  l3.SetTextSize(tsize), l3.SetTextAngle(10);
  l3.DrawLatex(22, 385, "E_{T}=90 GeV");

  l4.SetTextAlign(align);  l4.SetTextSize(tsize), l4.SetTextAngle(5);
  l4.DrawLatex(22, 220, "E_{T}=55 GeV");

  f1->Draw("same");
  f2->Draw("same");
  f3->Draw("same");
  if (prl) f4->Draw("same");
  gPad->RedrawAxis();
}
