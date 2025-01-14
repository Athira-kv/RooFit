//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 10 14:28:08 2024 by ROOT version 6.32.06
// from TTree tree/tree
// found on file: w8_k0_uDST.root
//////////////////////////////////////////////////////////

#ifndef fit_k_h
#define fit_k_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

class fit_k {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Evt;
   Int_t           TriggerMask;
   TLorentzVector  *lv_beam;
   TLorentzVector  *lv_scat;
   Double_t        Q2;
   Double_t        xbj;
   Double_t        y;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         v2x;
   Float_t         v2y;
   Float_t         v2z;
   Double_t        Chi2;
   TLorentzVector  *lv_p;
   TLorentzVector  *lv_pi;
   TLorentzVector  *lv_lambda;
   TLorentzVector  *lv_pip;
   TLorentzVector  *lv_pim;
   TLorentzVector  *lv_k0;
   Double_t        pt1;
   Double_t        pt2;
   Double_t        alpha;
   Double_t        pp_x;
   Double_t        pp_y;
   Double_t        pp_mom;
   Double_t        pp_ch;
   Double_t        pp_theta;
   Double_t        pm_x;
   Double_t        pm_y;
   Double_t        pm_mom;
   Double_t        pm_ch;
   Double_t        pm_theta;
   Double_t        pm_lh[7];
   Double_t        pp_lh[7];
   Double_t        k_thr;
   Double_t        pi_thr;
   Double_t        p_thr;
   Double_t        photons_p;
   Double_t        photons_m;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_TriggerMask;   //!
   TBranch        *b_lv_beam;   //!
   TBranch        *b_lv_scat;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_xbj;   //!
   TBranch        *b_y;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_v2x;   //!
   TBranch        *b_v2y;   //!
   TBranch        *b_v2z;   //!
   TBranch        *b_Chi2;   //!
   TBranch        *b_lv_p;   //!
   TBranch        *b_lv_pi;   //!
   TBranch        *b_lv_lambda;   //!
   TBranch        *b_lv_pip;   //!
   TBranch        *b_lv_pim;   //!
   TBranch        *b_lv_k0;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_pp_x;   //!
   TBranch        *b_pp_y;   //!
   TBranch        *b_pp_mom;   //!
   TBranch        *b_pp_ch;   //!
   TBranch        *b_pp_theta;   //!
   TBranch        *b_pm_x;   //!
   TBranch        *b_pm_y;   //!
   TBranch        *b_pm_mom;   //!
   TBranch        *b_pm_ch;   //!
   TBranch        *b_pm_theta;   //!
   TBranch        *b_pm_lh;   //!
   TBranch        *b_pp_lh;   //!
   TBranch        *b_k_thr;   //!
   TBranch        *b_pi_thr;   //!
   TBranch        *b_p_thr;   //!
   TBranch        *b_photons_p;   //!
   TBranch        *b_photons_m;   //!

   fit_k(TTree *tree=0);
   virtual ~fit_k();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fit_k_cxx
fit_k::fit_k(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../w8_k0_uDST.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../w8_k0_uDST.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

fit_k::~fit_k()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fit_k::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fit_k::LoadTree(Long64_t entry)
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

void fit_k::Init(TTree *tree)
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
   lv_p = 0;
   lv_pi = 0;
   lv_lambda = 0;
   lv_pip = 0;
   lv_pim = 0;
   lv_k0 = 0;
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
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("v2x", &v2x, &b_v2x);
   fChain->SetBranchAddress("v2y", &v2y, &b_v2y);
   fChain->SetBranchAddress("v2z", &v2z, &b_v2z);
   fChain->SetBranchAddress("Chi2", &Chi2, &b_Chi2);
   fChain->SetBranchAddress("lv_p", &lv_p, &b_lv_p);
   fChain->SetBranchAddress("lv_pi", &lv_pi, &b_lv_pi);
   fChain->SetBranchAddress("lv_lambda", &lv_lambda, &b_lv_lambda);
   fChain->SetBranchAddress("lv_pip", &lv_pip, &b_lv_pip);
   fChain->SetBranchAddress("lv_pim", &lv_pim, &b_lv_pim);
   fChain->SetBranchAddress("lv_k0", &lv_k0, &b_lv_k0);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("alpha", &alpha, &b_alpha);
   fChain->SetBranchAddress("pp_x", &pp_x, &b_pp_x);
   fChain->SetBranchAddress("pp_y", &pp_y, &b_pp_y);
   fChain->SetBranchAddress("pp_mom", &pp_mom, &b_pp_mom);
   fChain->SetBranchAddress("pp_ch", &pp_ch, &b_pp_ch);
   fChain->SetBranchAddress("pp_theta", &pp_theta, &b_pp_theta);
   fChain->SetBranchAddress("pm_x", &pm_x, &b_pm_x);
   fChain->SetBranchAddress("pm_y", &pm_y, &b_pm_y);
   fChain->SetBranchAddress("pm_mom", &pm_mom, &b_pm_mom);
   fChain->SetBranchAddress("pm_ch", &pm_ch, &b_pm_ch);
   fChain->SetBranchAddress("pm_theta", &pm_theta, &b_pm_theta);
   fChain->SetBranchAddress("pm_lh", pm_lh, &b_pm_lh);
   fChain->SetBranchAddress("pp_lh", pp_lh, &b_pp_lh);
   fChain->SetBranchAddress("k_thr", &k_thr, &b_k_thr);
   fChain->SetBranchAddress("pi_thr", &pi_thr, &b_pi_thr);
   fChain->SetBranchAddress("p_thr", &p_thr, &b_p_thr);
   fChain->SetBranchAddress("photons_p", &photons_p, &b_photons_p);
   fChain->SetBranchAddress("photons_m", &photons_m, &b_photons_m);
   Notify();
}

bool fit_k::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void fit_k::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fit_k::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef roofit_ana_cxx
