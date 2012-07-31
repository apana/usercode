struct FatHiggsInfo{
  Int_t FatHiggsFlag;
  Float_t mass;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t filteredmass;
  Float_t filteredpt;
  Float_t filteredeta;
  Float_t filteredphi;
};

struct HiggsInfo{
  int HiggsFlag;
  float mass;
  float pt;  
  float eta; 
  float phi; 
  float dR;  
  float dPhi;
  float dEta;
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

struct genParticleInfo
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
};
