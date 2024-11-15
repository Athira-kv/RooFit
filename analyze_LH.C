#define analyze_LH_cxx
#include "analyze_LH.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
double m_pi = 0.139570;


void analyze_LH::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH2D* h_arment = new TH2D("arment","arment", 2000, -1, 1, 400, 0, 0.4);
   TH1D* h_invmass_k = new TH1D("kmass","kmass",200,0.4,0.6);

   TH2I *h2 = new TH2I("h2", "Negative vs Positive Particle Identification;Negative Particle;Positive Particle", 4, -0.5, 3.5, 4, -0.5, 3.5);

   TH2I *pho_mom = new TH2I("h3", "photons;Positive Particle ID;Momentum ", 4, -0.5, 3.5, 350, 10, 80);

   TH2D* h4 = new TH2D("h4", "Difference vs photons; LH(#pi) - LH (K); #photons",100, -1, 3, 150, 0, 150);
   TH2D* h5 = new TH2D("h5", "Ratio Vs Difference with # photons ; LH ratio ; LH difference",150, 0,150,100,-1,3);
   TH1F* invMassHist[4][4][13];

   //momentum bins: Momentum p(GeV/c) = (10,11,12,13,15,17,19,22,25,27,30,35,40,50)

   std::vector<float> mombins = {10.0 ,11.0 ,12.0 ,13.0 ,15.0 ,17.0 ,19.0 ,22.0 ,25.0 ,27.0 ,30.0 ,35.0 ,40.0 ,50.0 };

   
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double n_rich = sqrt((m_pi*m_pi/(pi_thr*pi_thr)) + 1);
      if(pp_mom < pi_thr || pm_mom < pi_thr){ continue ; }


      if( pt1 < 0.10 || pt2 < 0.10){ continue ;}

      if(pp_mom < 10. || pm_mom < 10.){ continue ; }
      //if( pt1 < 0.1 | pt2 < 0.1 ){ continue ; }

      float theta_pp = pp_ch;
      float theta_pm = pm_ch;

      if(theta_pp == 0 || theta_pp == -1){ continue ; }
      if(theta_pm == 0 || theta_pm == -1){ continue ; }

      double k_mass = lv_k0->M();
      
      //      if( abs(k_mass - 0.497) > 0.020 ) { continue ; }
      if( lv_k0->Pt() < 0.040 ) { continue ; }

      double lh_values_p[4] = {pp_lh[0], pp_lh[1], pp_lh[2], pp_lh[5]};
      double lh_values_m[4] = {pm_lh[0], pm_lh[1], pm_lh[2], pm_lh[5]};

      float p_maxvalue = -999.0;
      int p_index = -10;
      float m_maxvalue = -999.0;
      int m_index = -10;

      bool proton = true;
      bool pion = true;

      for(int m = 0; m < 4; m++) {
	if(lh_values_p[m] > p_maxvalue) {
	  p_maxvalue = lh_values_p[m];
	  p_index = m;
        }
        if(lh_values_m[m] > m_maxvalue) {
	  m_maxvalue = lh_values_m[m];
	  m_index = m;
        }
      }

       size_t pBinIndex = 0;
    while (pBinIndex < mombins.size() - 1 && !(pp_mom >= mombins[pBinIndex] && pp_mom < mombins[pBinIndex + 1])) {
        ++pBinIndex;
    }

    if( m_index == 0  ){

      h4->Fill((lh_values_p[0]-lh_values_p[1]), photons_p);
      h5->Fill((lh_values_p[0]/lh_values_p[1]), (lh_values_p[0]-lh_values_p[1]), photons_p);
      pho_mom->Fill(p_index,pp_mom,photons_p);
      //     invMassHist[m_index][p_index][pBinIndex]->Fill(photons_p);
    }
    
      
    h2->Fill(m_index, p_index);
      h_invmass_k->Fill(k_mass);
      h_arment->Fill(alpha, pt1);

      
   }


    TFile *outputFile = new TFile("w8_k0_LH.root", "RECREATE");

   h_arment->Write();
   h_invmass_k->Write();
   h2->Write();
   pho_mom->Write();
   h4->Write();
   h5->Write();
  
   outputFile->Close();
   
}
