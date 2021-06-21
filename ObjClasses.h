// ---- Classes -------------------------------
class ClassTrack {
 public:
  //Double_t E; //
  //Double_t Et; //
  Double_t Pt; //
  //Double_t Mass; //
  Double_t Eta; //
  Double_t Phi; //
  TLorentzVector FourVector; //
  Double_t DR; //
  Double_t DRleadtrk; //
  Double_t z0; //
  Double_t d0; //
  Double_t BDTscore; //
  //Int_t ntracks; //
  Int_t istruth; //
  //Bool_t mtoLeadTau; //
  Int_t charge; //
  Int_t nPiHits; // 
  Int_t nSiHoles; //
  //Int_t mtoOfflTau;
  //Int_t mtoL1Tau;
  //Int_t mtoL1Tauv2;
  Int_t i_pt; //
  Int_t i_bdt; //
  Int_t i_aux; //
  //Int_t index_DR;
  //Int_t contained;
  bool operator==( const ClassTrack& other ) const {
    if ( (Phi == other.Phi) && (Eta == other.Eta) ) {
      return true;
    }
    return false;
  }
  
  ClassTrack ();
  ClassTrack (TLorentzVector);
  ~ClassTrack ();
};

class ClassTau {
 public:
  //Double_t E;
  //Double_t Et;
  Double_t Pt;
  //Double_t Mass;
  Double_t Eta;
  Double_t Phi;
  Double_t PtCalib;
  Int_t charge;
  TLorentzVector FourVector;
  Int_t Nprong;
  Int_t Nneutral;
  Int_t Ntracks;
  Double_t HLT_Pt;
  Double_t Calo_Pt;
  Double_t Calo_Eta;
  Double_t Calo_Phi;
  Double_t Calo_Had;
  Double_t Calo_EM;

  std::vector<ClassTrack> offlTracks;
  
  //ClassTrack lead1stgTrack;
  //Int_t matched_lead1stgTrack;
  
  //ClassTrack close1stgTrack;
  //Int_t matched_close1stgTrack;
  
  //ClassTrack leadofflTrack;
  //Int_t matched_leadofflTrack;
  
  //ClassTrack closeofflTrack;
  //Int_t matched_closeofflTrack;
  
  //ClassTrack leadprecTrack;
  //Int_t matched_leadprecTrack;
  
  //ClassTrack closeprecTrack;
  //Int_t matched_closeprecTrack;

  bool hastruthtrk;
  //bool hasL1tau;
  //Int_t L1tau;
  //Int_t L1tausub;
  //Int_t matched_L1tau;
  //std::vector<ClassTrack> stg1Tracks_m;
  //std::vector<ClassTrack> stg1Tracks_c;
  //ClassTrack target1stg_minDR;
  //ClassTrack target1stg_minDz;
  
  bool operator==( const ClassTau& other ) const {
    if ( (Phi == other.Phi) && (Eta == other.Eta) ) {
      return true;
    }
    return false;
  }

  ClassTau ();
  ClassTau (TLorentzVector);
  ~ClassTau ();
};



// Constructors ---------------------------------------------- 
ClassTau::ClassTau () {};

ClassTau::ClassTau ( TLorentzVector FourVec ) {
  //E = FourVec.E();
  Pt = FourVec.Pt();
  Eta = FourVec.Eta();
  Phi = FourVec.Phi();
  //Et = FourVec.Et();
  //Mass = FourVec.M();
  FourVector = FourVec;
  //matched_lead1stgTrack= -1;
  //matched_close1stgTrack= -1;
  //matched_leadofflTrack= -1;
  //matched_closeofflTrack= -1;
  //matched_leadprecTrack= -1;
  //matched_closeprecTrack= -1;
  Ntracks = -1;
  //HLT_Pt = -2;
  Calo_Pt = -2;
  Calo_Eta = -99;
  Calo_Phi = -99;
  Calo_Had = -2;
  Calo_EM = -2;
  Nprong = -1;
  Nneutral = -1;

  hastruthtrk = false;
  //hasL1tau = false;
  //matched_L1tau=0;
  //L1tau=-1;
  //L1tausub=-1;
  PtCalib=-1.;
  charge=-99;
};

ClassTau::~ClassTau () {};

ClassTrack::ClassTrack () {};

ClassTrack::ClassTrack ( TLorentzVector FourVec ) {
  //E = FourVec.E();
  Pt = FourVec.Pt();
  Eta = FourVec.Eta();
  Phi = FourVec.Phi();
  //Et = FourVec.Et();
  //Mass = FourVec.M();
  FourVector = FourVec;
  DR=-999;
  DRleadtrk=-999;
  //ntracks=-1;
  //contained = -1;
  //charge=0;
  nPiHits=-1;
  nSiHoles=-1;
  i_pt = -1;
  i_bdt = -1;
  //mtoOfflTau=0;
  //mtoLeadTau=false;
  istruth = -1;
  //BDTscore = -999;
};

ClassTrack::~ClassTrack () {};

// ---- Comparators ---------------------------------------
 
struct ComparePt {
  bool operator()(TLorentzVector a, TLorentzVector b) {
    return (a.Pt() > b.Pt());
  }
  bool operator()(ClassTau a, ClassTau b) {
    return (a.Pt > b.Pt);
  }
    bool operator()(ClassTrack a, ClassTrack b) {
    return (a.Pt > b.Pt);
  }
};

struct CompareEta {
  bool operator()(TLorentzVector a, TLorentzVector b) {
    return (a.Eta() < b.Eta());
  }
  bool operator()(ClassTau a, ClassTau b) {
    return (a.Eta < b.Eta);
  }
    bool operator()(ClassTrack a, ClassTrack b) {
    return (a.Eta < b.Eta);
  }
};

struct ComparePhi {
  bool operator()(TLorentzVector a, TLorentzVector b) {
    return (a.Phi() < b.Phi());
  }
  bool operator()(ClassTau a, ClassTau b) {
    return (a.Phi < b.Phi);
  }
  bool operator()(ClassTrack a, ClassTrack b) {
   return (a.Phi < b.Phi);
  }
};

struct CompareBDT {
  bool operator()(ClassTrack a, ClassTrack b) {
   return (a.BDTscore > b.BDTscore);
  }
};

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}