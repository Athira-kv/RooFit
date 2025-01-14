//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  1 14:09:38 2024 by ROOT version 6.32.06
// from TTree tree2/tree2
// found on file: ../efficiency/w8_rich_phi.root
//////////////////////////////////////////////////////////

#ifndef fit_phi_h
#define fit_phi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

class fit_phi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Long64_t        Evt;
   Int_t           TriggerMask;
   TLorentzVector  *lv_beam;
   TLorentzVector  *lv_scat;
   Double_t        Q2;
   Double_t        xbj;
   Double_t        y;
   Double_t        Emiss;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   TLorentzVector  *lv_kp;
   TLorentzVector  *lv_km;
   TLorentzVector  *lv_phi;
   Double_t        pt1;
   Double_t        pt2;
   Double_t        alpha;
   Double_t        pp_x;
   Double_t        pp_y;
   Double_t        pp_mom;
   Double_t        pp_theta;
   Double_t        pm_x;
   Double_t        pm_y;
   Double_t        pm_mom;
   Double_t        pm_theta;
   Double_t        pm_lh[7];
   Double_t        pp_lh[7];
   Double_t        k_thr;
   Double_t        pi_thr;
   Double_t        p_thr;
   Int_t           Nout;
   Int_t           charge1;
   Int_t           charge2;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_TriggerMask;   //!
   TBranch        *b_lv_beam;   //!
   TBranch        *b_lv_scat;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_xbj;   //!
   TBranch        *b_y;   //!
   TBranch        *b_Emiss;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_lv_kp;   //!
   TBranch        *b_lv_km;   //!
   TBranch        *b_lv_phi;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_pp_x;   //!
   TBranch        *b_pp_y;   //!
   TBranch        *b_pp_mom;   //!
   TBranch        *b_pp_theta;   //!
   TBranch        *b_pm_x;   //!
   TBranch        *b_pm_y;   //!
   TBranch        *b_pm_mom;   //!
   TBranch        *b_pm_theta;   //!
   TBranch        *b_pm_lh;   //!
   TBranch        *b_pp_lh;   //!
   TBranch        *b_k_thr;   //!
   TBranch        *b_pi_thr;   //!
   TBranch        *b_p_thr;   //!
   TBranch        *b_Nout;   //!
   TBranch        *b_charge1;   //!
   TBranch        *b_charge2;   //!

   fit_phi(TTree *tree=0);
   virtual ~fit_phi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fit_phi_cxx
fit_phi::fit_phi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../w8_rich_phi.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../w8_rich_phi.root");
      }
      f->GetObject("tree2",tree);

   }
   Init(tree);
}

fit_phi::~fit_phi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fit_phi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fit_phi::LoadTree(Long64_t entry)
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

void fit_phi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   lv_beam = 0;
   lv_scat = 0;
   lv_kp = 0;
   lv_km = 0;
   lv_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("TriggerMask", &TriggerMask, &b_TriggerMask);
   fChain->SetBranchAddress("lv_beam", &lv_beam, &b_lv_beam);
   fChain->SetBranchAddress("lv_scat", &lv_scat, &b_lv_scat);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("xbj", &xbj, &b_xbj);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("Emiss", &Emiss, &b_Emiss);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("lv_kp", &lv_kp, &b_lv_kp);
   fChain->SetBranchAddress("lv_km", &lv_km, &b_lv_km);
   fChain->SetBranchAddress("lv_phi", &lv_phi, &b_lv_phi);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("alpha", &alpha, &b_alpha);
   fChain->SetBranchAddress("pp_x", &pp_x, &b_pp_x);
   fChain->SetBranchAddress("pp_y", &pp_y, &b_pp_y);
   fChain->SetBranchAddress("pp_mom", &pp_mom, &b_pp_mom);
   fChain->SetBranchAddress("pp_theta", &pp_theta, &b_pp_theta);
   fChain->SetBranchAddress("pm_x", &pm_x, &b_pm_x);
   fChain->SetBranchAddress("pm_y", &pm_y, &b_pm_y);
   fChain->SetBranchAddress("pm_mom", &pm_mom, &b_pm_mom);
   fChain->SetBranchAddress("pm_theta", &pm_theta, &b_pm_theta);
   fChain->SetBranchAddress("pm_lh", pm_lh, &b_pm_lh);
   fChain->SetBranchAddress("pp_lh", pp_lh, &b_pp_lh);
   fChain->SetBranchAddress("k_thr", &k_thr, &b_k_thr);
   fChain->SetBranchAddress("pi_thr", &pi_thr, &b_pi_thr);
   fChain->SetBranchAddress("p_thr", &p_thr, &b_p_thr);
   fChain->SetBranchAddress("Nout", &Nout, &b_Nout);
   fChain->SetBranchAddress("charge1", &charge1, &b_charge1);
   fChain->SetBranchAddress("charge2", &charge2, &b_charge2);
   Notify();
}

bool fit_phi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void fit_phi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fit_phi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef phi_roofit_mom_cxx
