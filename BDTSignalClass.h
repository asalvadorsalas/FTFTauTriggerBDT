//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 12 15:12:54 2018 by ROOT version 6.12/06
// from TTree tree_vec/tree_vec
// found on file: ../FTF_Adrian/run/run_1dd77f8d986744aca3e9436ba497acd6/data-hist/sample.root
//////////////////////////////////////////////////////////

#ifndef BDTSignalClass_h
#define BDTSignalClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
// Header file for the classes stored in the TTree if any.
#include "vector"
#include <TLorentzVector.h>
#include "ObjClasses.h"

class BDTSignalClass {
public :
  
  TTree *tree;
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  TFile          *hfile;

  std::vector<ClassTau> OfflineTaus;

  // Declaration of leaf types
  Int_t           EventNumber;
  Int_t           passL1_Tau12IM;

  std::vector<float>   *Coretrack_pt;
  std::vector<float>   *Coretrack_eta;
  std::vector<float>   *Coretrack_phi;
  std::vector<float>   *Coretrack_z0;
  std::vector<float>   *tau_offl_pt;
  std::vector<float>   *tau_offl_ptFinalCalib;
  std::vector<int>     *tau_offl_isMediumRNN;
  std::vector<float>   *tau_offl_eta;
  std::vector<float>   *tau_offl_phi;
  std::vector<int>     *tau_offl_charge;
  std::vector<int>     *tau_truth_Nprong;
  std::vector<int>     *tau_truth_Nneutral;
  std::vector<int>     *tau_nTracks;
  std::vector<float>   *tau_pt;
  std::vector<float>   *tau_ptTrigCalo;
  std::vector<float>   *tau_etaTrigCalo;
  std::vector<float>   *tau_phiTrigCalo;
  std::vector<float>   *tau_ptHad;
  std::vector<float>   *tau_ptEM;
  std::vector<std::vector<float> > *offltrack_pt;
  std::vector<std::vector<float> > *offltrack_eta;
  std::vector<std::vector<float> > *offltrack_phi;
  std::vector<std::vector<float> > *offltrack_z0;
  std::vector<std::vector<float> > *offltrack_d0;
  std::vector<std::vector<int> > *offltrack_nPiHits;
  std::vector<std::vector<int> > *offltrack_nSiHoles;
  std::vector<std::vector<int> >   *offltrack_charge;
  std::vector<std::vector<int> > *offltrack_truthtype;

  // List of branches
  TBranch        *b_EventNumber;   //!
  TBranch        *b_passL1_Tau12IM;   //!

  TBranch        *b_Coretrack_pt; //!
  TBranch        *b_Coretrack_eta; //!
  TBranch        *b_Coretrack_phi; //!
  TBranch        *b_Coretrack_z0; //!
  TBranch        *b_tau_offl_pt; //!
  TBranch        *b_tau_offl_ptFinalCalib;   //!
  TBranch        *b_tau_offl_isMediumRNN; //!
  TBranch        *b_tau_offl_eta; //!
  TBranch        *b_tau_offl_phi; //!
  TBranch        *b_tau_offl_charge; //!
  TBranch        *b_tau_truth_Nprong; //!
  TBranch        *b_tau_truth_Nneutral; //!
  TBranch        *b_tau_nTracks; //!
  TBranch        *b_tau_pt; //!
  TBranch        *b_tau_ptTrigCalo; //!
  TBranch        *b_tau_etaTrigCalo; //!
  TBranch        *b_tau_phiTrigCalo; //!
  TBranch        *b_tau_ptHad; //!
  TBranch        *b_tau_ptEM; //!
  TBranch        *b_offltrack_pt; //!
  TBranch        *b_offltrack_eta; //!
  TBranch        *b_offltrack_phi; //!
  TBranch        *b_offltrack_z0; //!
  TBranch        *b_offltrack_d0; //!
  TBranch        *b_offltrack_nPiHits; //!
  TBranch        *b_offltrack_nSiHoles; //!
  TBranch        *b_offltrack_charge; //!
  TBranch        *b_offltrack_truthtype; //!


  std::vector<int>     *Coretrack_signalv0;
  TBranch        *b_Coretrack_signalv0;   //!
  std::vector<int>     *Coretrack_signalv1;
  TBranch        *b_Coretrack_signalv1;   //!

  BDTSignalClass(const char* inputfile, const char* outputfile);
  virtual ~BDTSignalClass();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     DeclareHistos();
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef BDTSignalClass_cxx

Int_t BDTSignalClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BDTSignalClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BDTSignalClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Coretrack_pt = 0;
   Coretrack_eta = 0;
   Coretrack_phi = 0;
   Coretrack_z0 = 0;
   tau_offl_pt = 0;
   tau_offl_ptFinalCalib = 0;
   tau_offl_isMediumRNN = 0;
   tau_offl_eta = 0;
   tau_offl_phi = 0;
   tau_offl_charge = 0;
   tau_truth_Nprong = 0;
   tau_truth_Nneutral = 0;
   tau_nTracks = 0;
   tau_pt = 0;
   tau_ptTrigCalo = 0;
   tau_etaTrigCalo = 0;
   tau_phiTrigCalo = 0;
   tau_ptHad = 0;
   tau_ptEM = 0;
   offltrack_pt = 0;
   offltrack_eta = 0;
   offltrack_phi = 0;
   offltrack_z0 = 0;
   offltrack_d0 = 0;
   offltrack_nPiHits = 0;
   offltrack_nSiHoles = 0;
   offltrack_charge = 0;
   offltrack_truthtype = 0;


   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("passL1_Tau12IM", &passL1_Tau12IM, &b_passL1_Tau12IM);

   fChain->SetBranchAddress("Coretrack_pt", &Coretrack_pt, &b_Coretrack_pt);
   fChain->SetBranchAddress("Coretrack_eta", &Coretrack_eta, &b_Coretrack_eta);
   fChain->SetBranchAddress("Coretrack_phi", &Coretrack_phi, &b_Coretrack_phi);
   fChain->SetBranchAddress("Coretrack_z0", &Coretrack_z0, &b_Coretrack_z0);
   fChain->SetBranchAddress("tau_offl_pt", &tau_offl_pt, &b_tau_offl_pt);
   fChain->SetBranchAddress("tau_offl_ptFinalCalib", &tau_offl_ptFinalCalib, &b_tau_offl_ptFinalCalib);
   fChain->SetBranchAddress("tau_offl_isMediumRNN", &tau_offl_isMediumRNN, &b_tau_offl_isMediumRNN);
   fChain->SetBranchAddress("tau_offl_eta", &tau_offl_eta, &b_tau_offl_eta);
   fChain->SetBranchAddress("tau_offl_phi", &tau_offl_phi, &b_tau_offl_phi);
   fChain->SetBranchAddress("tau_offl_charge", &tau_offl_charge, &b_tau_offl_charge);
   fChain->SetBranchAddress("tau_truth_Nprong", &tau_truth_Nprong, &b_tau_truth_Nprong);
   fChain->SetBranchAddress("tau_truth_Nneutral", &tau_truth_Nneutral, &b_tau_truth_Nneutral);
   fChain->SetBranchAddress("tau_nTracks", &tau_nTracks, &b_tau_nTracks);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_ptTrigCalo", &tau_ptTrigCalo, &b_tau_ptTrigCalo);
   fChain->SetBranchAddress("tau_etaTrigCalo", &tau_etaTrigCalo, &b_tau_etaTrigCalo);
   fChain->SetBranchAddress("tau_phiTrigCalo", &tau_phiTrigCalo, &b_tau_phiTrigCalo);
   fChain->SetBranchAddress("tau_ptHad", &tau_ptHad, &b_tau_ptHad);
   fChain->SetBranchAddress("tau_ptEM", &tau_ptEM, &b_tau_ptEM);
   fChain->SetBranchAddress("offltrack_pt", &offltrack_pt, &b_offltrack_pt);
   fChain->SetBranchAddress("offltrack_eta", &offltrack_eta, &b_offltrack_eta);
   fChain->SetBranchAddress("offltrack_phi", &offltrack_phi, &b_offltrack_phi);
   fChain->SetBranchAddress("offltrack_z0", &offltrack_z0, &b_offltrack_z0);
   fChain->SetBranchAddress("offltrack_d0", &offltrack_d0, &b_offltrack_d0);
   fChain->SetBranchAddress("offltrack_nPiHits", &offltrack_nPiHits, &b_offltrack_nPiHits);
   fChain->SetBranchAddress("offltrack_nSiHoles", &offltrack_nSiHoles, &b_offltrack_nSiHoles);
   fChain->SetBranchAddress("offltrack_charge", &offltrack_charge, &b_offltrack_charge);
   fChain->SetBranchAddress("offltrack_truthtype", &offltrack_truthtype, &b_offltrack_truthtype);

   Coretrack_signalv0 = 0;
   fChain->SetBranchAddress("Coretrack_signalv0", &Coretrack_signalv0, &b_Coretrack_signalv0);  
   Coretrack_signalv1 = 0;
   fChain->SetBranchAddress("Coretrack_signalv1", &Coretrack_signalv1, &b_Coretrack_signalv1);   
  
   Notify();
}

Bool_t BDTSignalClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BDTSignalClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BDTSignalClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif