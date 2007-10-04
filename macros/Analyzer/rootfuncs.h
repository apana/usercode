TFile* OpenRootFile(const TString& rootfile) {

  cout << "Rootfile: " << rootfile << endl;

  TFile *file;
  if( gSystem->AccessPathName(rootfile) ){
    cout << endl << "File: " << rootfile << " not there!!!" << endl << endl;
    file=0;
  }
  else
    file = new TFile(rootfile);

  return file;
}

TFile* OpenDCacheFile(const TString& rootfile) {

  cout << "Rootfile: " << rootfile << endl;

  TFile *file;
  if( gSystem->AccessPathName(rootfile) ){
    cout << endl << "File: " << rootfile << " not there!!!" << endl << endl;
    file=0;
  }
  else
    file = new TDCacheFile(rootfile,"READ","Demo ROOT file with histograms",0);

  return file;
}

bool hExist(TFile *file,const TString& hname){

  bool retval=true;

  TKey *key = file->FindKey(hname);
  if (key ==0){
    cout << "!!Histogram " << hname << " does not exist!!" << endl;
    retval=false;
  }
  return retval;
}

TH1* GetHist(TFile* file,const TString& hname){

  TH1* h=0;
  cout << "Getting histogram: " << hname << endl;
  if (hExist(file,hname)){
    h=(TH1F*)file->Get(hname);
  }

  return h;
}

TH2* GetHist2D(TFile* file,const TString& hname){

  TH2* h;
  if (hExist(file,hname)){
    h=(TH2F*)file->Get(hname);
  }

  return h;
}
