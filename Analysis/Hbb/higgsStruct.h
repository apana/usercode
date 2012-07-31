struct FatHiggsInfo{
  Bool_t FatHiggsFlag;
  Float_t mass;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t filteredmass;
  Float_t filteredpt;
};


struct HiggsInfo{
  Float_t mass;
  Float_t pt;  
  Float_t eta; 
  Float_t phi; 
  Float_t dR;  
  Float_t dPhi;
  Float_t dEta;
};

struct VInfo{
  Float_t mass;
  Float_t pt;  
  Float_t eta; 
  Float_t phi; 
};

struct EventInfo{
  Int_t run;
  Int_t lumi;
  Int_t event;
  Bool_t json;
};

struct METInfo
{
  Float_t et; 
  Float_t sumet;   
  Float_t sig;
  Float_t phi;
};
