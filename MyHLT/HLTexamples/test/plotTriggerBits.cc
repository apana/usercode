void plotTriggerBits(const int i1=1, const int i2=35){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy(0);

  TString rootfile="histo.root";
  cout << "\nRootfile: " << rootfile << "\n\n";
  TFile *_file0 = TFile::Open(rootfile);
  _file0->cd("plots");
  // gDirectory->ls();

  TString hname="TriggerResults";
  TH1F* h;
  gDirectory->GetObject(hname,h);


  TCanvas *c1 = new TCanvas("c1","Root Canvas 1",400, 52, 1100, 550);
  float bmargin=0.16, rmargin=0.07, lmargin=0.055, tmargin=0.02;

  c1->SetBottomMargin(bmargin);
  c1->SetTopMargin(tmargin);
  c1->SetRightMargin(rmargin);
  c1->SetLeftMargin(lmargin);
  h->SetFillColor(kYellow);

  const char * binlabel1=h->GetXaxis()->GetBinLabel(i1);
  const char * binlabel2=h->GetXaxis()->GetBinLabel(i2);

  cout << "\nPlotting trigger occupancy for bits " << i1 << "-" << i2 
       << " (" << binlabel1 << " - " << binlabel2 << ")" << "\n\n";

  h->GetXaxis()->SetRange(i1,i2);
  h->Draw();
}
