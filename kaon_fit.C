#define kaon_fit_cxx
#include "kaon_fit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream> 

double m_pi = 0.139570;
double m_k = 0.493677;

TProfile2D* logplot(const double xMin, const double yMin, const double xMax, const double yMax, const char* name = "hist") {

    const int nBinsX = 100;
    const int nBinsY = 100;

    double* binEdgesX = new double[nBinsX + 1];
    double* binEdgesY = new double[nBinsY + 1];

    for (int i = 0; i <= nBinsX; ++i) {
        binEdgesX[i] = xMin * TMath::Power(xMax/xMin, static_cast<double>(i)/nBinsX);
    }

    for (int i = 0; i <= nBinsY; ++i) {
        binEdgesY[i] = yMin * TMath::Power(yMax/yMin, static_cast<double>(i)/nBinsY);
    }

    TProfile2D* hist = new TProfile2D(name," ", nBinsX, binEdgesX, nBinsY, binEdgesY);

    return hist;
}

void kaon_fit::Loop()
{

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH2D* h_arment = new TH2D("arment","arment", 2000, -1, 1, 400, 0, 0.4);
   TH1D* h_invmass_k = new TH1D("kmass","kmass",200,0.4,0.6);
   TH2D* h_x_y = new TH2D("xy", "xy", 300, -150, 150, 300, -150, 150);
   TH1D* h_mom = new TH1D("mom", "mom", 100, 0, 10);

   TProfile2D* h_0_1 = logplot(1e-3, 1e-3, 4, 4, "h_0_1");
   TProfile2D* h_0_2 = logplot(1e-3, 1e-3, 4, 4, "h_0_2");
   TProfile2D* h_0_3 = logplot(1e-3, 1e-3, 4, 4, "h_0_3");

   
   TProfile2D* h3D = logplot(1e-3, 1, 4, 60.0, "h3D");

    const int nBinsX = 100;
   const int nBinsY = 100;
   double xMin = 1e-5;  double xMax = 1.;
   double yMin = 1e-3;
   double yMax = 100.;

   double* binEdgesX = new double[nBinsX + 1];
   for (int i = 0; i <= nBinsX; ++i) {
     binEdgesX[i] = xMin * TMath::Power(xMax/xMin, static_cast<double>(i)/nBinsX);
   }

   double* binEdgesY = new double[nBinsY + 1];
   for (int i = 0; i <= nBinsY; ++i) {
     binEdgesY[i] = yMin * TMath::Power(yMax/yMin, static_cast<double>(i)/nBinsY);
   }

   TH2F* hist = new TH2F("hist", "", nBinsX, binEdgesX, nBinsY, binEdgesY);


   std::vector<float> mombins = {2.0,5.0, 10.0 ,11.0 ,12.0 ,13.0 ,15.0 ,17.0 ,19.0 ,22.0 ,25.0 ,27.0 ,30.0 ,35.0 ,40.0 ,50.0 };
   size_t numMomentumBins = mombins.size();
   

   std::vector<float> thetabins = {0, 0.01, 0.04, 0.12, 0.3};
   size_t numthetabins = thetabins.size();
   
   std::vector<std::vector<std::vector<TH1F*>>> invMassHist(6,std::vector<std::vector<TH1F*>>(numMomentumBins, std::vector<TH1F*>(numthetabins, nullptr)  ) ) ;

   std::vector<std::vector<TH2D*>> IDhistograms(numMomentumBins, std::vector<TH2D*>(numthetabins, nullptr));


   for (int p = 0; p < 6; ++p) {
     for (size_t i = 0; i < numMomentumBins; ++i) {
       double min_p = mombins.at(i);
       for (size_t j = 0; j < numthetabins; ++j) {

	 double min_theta = thetabins.at(j);
	 TString histName = Form("invMass_%d_mom_%.f_theta_%.2f", p, min_p, min_theta);
	 TString histTitle = Form("Invariant Mass for p_index = %d, mom %.f, theta %.2f",
				  p, min_p, min_theta);
	 invMassHist[p][i][j] = new TH1F(histName, histTitle, 200, 0.4, 0.6);	 
       }
     }
   }
   
   for (size_t i = 0; i < numMomentumBins; ++i) {
     double min_p = mombins.at(i);
     for (size_t j = 0; j < numthetabins; ++j) {
       
       double min_theta = thetabins.at(j);
       TString histName = Form("mom_%.f_theta_%.2f", min_p, min_theta);
       TString histTitle = Form(" mom %.f, theta %.2f", min_p, min_theta);
       IDhistograms[i][j] = new TH2D(histName, histTitle, 4, -1.5, 2.5, 4, -1.5, 2.5);
     }
   }


   
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double n_rich = sqrt((m_pi*m_pi/(pi_thr*pi_thr)) + 1);
      if(pp_mom < pi_thr || pm_mom < pi_thr){ continue ; }


      if(pp_mom < 2. || pm_mom < 2.){ continue ; }
    
      float theta_pp = pm_ch;
      float theta_pm = pp_ch;

      if(theta_pp == 0 || theta_pp == -1){ continue ; }
      if(theta_pm == 0 || theta_pm == -1){ continue ; }

      // x and y cordinates at the RICH entrance should lead to radius > 5 cm

      if((pp_x*pp_x + pp_y*pp_y ) <= 25 ){ continue; }
      if((pm_x*pm_x + pm_y*pm_y ) <= 25 ){ continue; }
      
      
      double k_mass = lv_k0->M();

      double lh_values_p[4] = {pp_lh[0], pp_lh[1], pp_lh[2], pp_lh[5]};
      double lh_values_m[4] = {pm_lh[0], pm_lh[1], pm_lh[2], pm_lh[5]};

      int p_index = -1; // Unknown by default
      int m_index = -1; // Unknown by default

      // Determine p_index
      if (pp_lh[0] / pp_lh[1] > 1.0 && pp_lh[0] / pp_lh[2] > 1.0 && pp_lh[0] / pp_lh[5] > 1.0) {
        p_index = 0;
      } else if (pp_lh[1] / pp_lh[0] > 1.08 && pp_lh[1] / pp_lh[2] > 1.0 && pp_lh[1] / pp_lh[5] > 1.24){
        p_index = 1;
      }
      
      if(pp_mom > p_thr + 5.0){
	if(pp_lh[2] / pp_lh[0] > 1.0 && pp_lh[2] / pp_lh[1] > 1.0 && pp_lh[2] / pp_lh[5] > 1.0) {
	  p_index = 2;
	}
      }

      if(pp_mom < p_thr - 5. && pp_mom > k_thr){
	if(pp_lh[0] / pp_lh[5] < 2.0 && pp_lh[1] / pp_lh[5] < 3.0) {
          p_index = 2;
        }
      }

      if( pp_mom > p_thr - 5. && pp_mom < p_thr + 5.){
	if(pp_lh[2] / pp_lh[0] > 1.0 && pp_lh[2] / pp_lh[1] > 1.0 && pp_lh[2] / pp_lh[5] > 1.0 && pp_lh[0] / pp_lh[5] < 2.0 && pp_lh[1] / pp_lh[5] < 3.0){
	  p_index = 2;
	}
      }
      
      if(pp_lh[3]/pp_lh[0] > 1.8 ){p_index = -11; }
      
      // Determine m_index
      if (pm_lh[0] / pm_lh[1] > 1.0 && pm_lh[0] / pm_lh[2] > 1.0 && pm_lh[0] / pm_lh[5] > 1.0) {
        m_index = 0;
      } else if (pm_lh[1] / pm_lh[0] > 1.08 && pm_lh[1] / pm_lh[2] > 1.0 && pm_lh[1] / pm_lh[5] > 1.24){
        m_index = 1;
      }
      
      if(pm_mom > p_thr + 5.0){
        if(pm_lh[2] / pm_lh[0] > 1.0 && pm_lh[2] / pm_lh[1] > 1.0 && pm_lh[2] / pm_lh[5] > 1.0) {
          m_index = 2;
        }
      }
      
      if(pm_mom < p_thr - 5. && pm_mom > k_thr){
        if(pm_lh[0] / pm_lh[5] < 2.0 && pm_lh[1] / pm_lh[5] < 3.0) {
          m_index = 2;
        }
      }
      
      if( pm_mom > p_thr - 5. && pm_mom < p_thr + 5.){
        if(pm_lh[2] / pm_lh[0] > 1.0 && pm_lh[2] / pm_lh[1] > 1.0 && pm_lh[2] / pm_lh[5] > 1.0 && pm_lh[0] / pm_lh[5] < 2.0 && pm_lh[1] / pm_lh[5] < 3.0){
          m_index = 2;
        }
      }
      
      if(pm_lh[3]/pm_lh[0] > 1.8 ){ m_index = -11; }

    
      size_t pBinIndex = 0;
      while (pBinIndex < mombins.size() - 1 && !(pp_mom >= mombins[pBinIndex] && pp_mom < mombins[pBinIndex + 1])) {
        ++pBinIndex;
      }
      
      size_t thetaBinIndex = 0;
      while (thetaBinIndex < thetabins.size() - 1 && !(pp_theta >= thetabins[thetaBinIndex] && pp_theta < thetabins[thetaBinIndex + 1])) {
	++thetaBinIndex;
      }

      
      
      IDhistograms[pBinIndex][thetaBinIndex]->Fill(m_index, p_index);
      
      if(m_index == 0){

	//if(invMassHist[4][pBinIndex][thetaBinIndex]){ std::cout<<" true for m=0 "<<std::endl;}
	h_arment->Fill(alpha, pt1);
	h_invmass_k->Fill(k_mass);
	h_x_y->Fill(pp_x, pp_y);
	h_0_1->Fill(pp_lh[1],pp_lh[0],pp_mom);
	h_0_2->Fill(pp_lh[2], pp_lh[0], pp_mom);
	h_0_3->Fill(pp_lh[3], pp_lh[0], pp_mom);
	h3D->Fill(pp_lh[0], pp_mom, photons_p);

        invMassHist[4][pBinIndex][thetaBinIndex]->Fill(k_mass);
	hist->Fill(xbj, Q2);
      }
      
      if(m_index == 0 && p_index >= 0 ){

	//if(invMassHist[p_index][pBinIndex][thetaBinIndex]){ std::cout<<" true for m = 0 & p != -10 "<<std::endl;}
	//else{std::cout<<" NOT true for m = 0 & p != -10 "<<std::endl;}
	
        invMassHist[p_index][pBinIndex][thetaBinIndex]->Fill(k_mass);
	
      }
      if(m_index == 0 && p_index == -1 ){
	//if( invMassHist[5][pBinIndex][thetaBinIndex]){ std::cout<<" true for m = 0 & p = -10 "<<std::endl ; }
	//else{std::cout<<" NOT true for m = 0 & p == -10 "<<std::endl;}
	invMassHist[5][pBinIndex][thetaBinIndex]->Fill(k_mass);     

      }

      if(m_index == 0 && p_index == -11 ){
	h_mom->Fill(pp_mom);
      }
      
   }


   TFile *outputFile = new TFile("k0_input_pr.root", "RECREATE");

   h_arment->Write();
   h_invmass_k->Write();
   hist->Write();
   h_mom->Write();
   h_x_y->Write();
   hist->Write();
   h_0_1->Write();
   h_0_2->Write();
   h_0_3->Write();
   h3D->Write();

     for (int p = 0; p < 6; ++p) {
     for (size_t i = 0; i < numMomentumBins; ++i) {
       for (size_t j = 0; j < numthetabins; ++j) {
	 invMassHist[p][i][j]->Write();
       }
     }
     }

     for (size_t i = 0; i < numMomentumBins; ++i) {
       for (size_t j = 0; j < numthetabins; ++j) {
         IDhistograms[i][j]->Write();
       }
     }
   
     outputFile->Close();

     

     // RooFitResult* fit_results[6][numMomentumBins][numthetabins] = {nullptr};

     std::vector<std::vector<std::vector<RooFitResult*>>> fit_results(
    6,   std::vector<std::vector<RooFitResult*>>(numMomentumBins, 
        std::vector<RooFitResult*>(numthetabins, nullptr) ));
     
   std::ofstream covdata("kmDST_pr_covmatrix.txt");

   std::ofstream covvalues("kmDST_pr_covvalues.txt");

   covvalues<<"Category\tMomentum\t angle\tN_as\tN_pi_s\tN_k_s\tN_p_s\tN_u_s\tN_unk_s\tindex1\tindex2\tvalue\tQ2\n";
   
   std::ofstream outputtext("kmDST_ang_signal_background_results_pr_simfit.txt");

   outputtext << "Category \t Momentum Bin \t Angle \t N_as \t\t N_pi_s \t N_k_s \t\t N_p_s \t\t N_u_s \t\t N_unk_s\n";
  
   TCanvas* c1 = new TCanvas("c1", "Fit Results", 800, 600);
   c1->Print("kmDST_ang_momfit_results_pr_sim.pdf[");

   for (int p_index = 0; p_index < 6; ++p_index) { 

     for (int pBinIndex = 0; pBinIndex < mombins.size() - 1; ++pBinIndex) {

       for(int thetaindex = 0; thetaindex < thetabins.size() - 1; thetaindex ++ ){ 
	 
	 RooRealVar mass("mass", "Invariant Mass of Pion Pairs", 0.45, 0.55);  // in GeV
	 
	 RooRealVar mean1("mean1", "Kaon Mass Peak 1", 0.49, 0.4, 0.5);  // First Gaussian (main peak)
	 RooRealVar sigma1("sigma1", "Width 1", 0.003, 0.001, 0.01);
	 RooGaussian gauss1("gauss1", "First Gaussian", mass, mean1, sigma1);

	 RooRealVar mean2("mean2", "Kaon Mass Peak 2", 0.49, 0.4, 0.5);  // Second Gaussian (wider tail)
	 RooRealVar sigma2("sigma2", "Width 2", .01, 0.01, 0.05);
	 RooGaussian gauss2("gauss2", "Second Gaussian", mass, mean2, sigma2);

	 RooRealVar frac("frac", "Fraction of First Gaussian", 0.5, 0.0, 1.0);

	 RooAddPdf signal("signal", "Double Gaussian Signal", RooArgList(gauss1, gauss2), RooArgList(frac));

	 RooRealVar a0("a0", "Coefficient of 1st term", 0.0);
	 RooRealVar a1("a1", "Coefficient of 2nd term", 0.0, -0.5, 0.5);
	 RooRealVar a2("a2", "Coefficient of 3rd term", 0.0, -0.5, 0.5);
	 RooRealVar a3("a3", "Coefficient of 4rd term", 0.0);
	 RooChebychev background("background", "Background Polynomial", mass, RooArgSet(a0, a1, a2, a3));


	 // extra part for error using covariant matrix

	 // Number of signal and background events for each category
	 int ent = invMassHist[4][pBinIndex][thetaindex]->GetEntries();
	 RooRealVar N_a_b("N_a_b", "N_a_b", 0.05*ent, 0., 1.05*ent);

	 ent = invMassHist[0][pBinIndex][thetaindex]->GetEntries();
	 RooRealVar N_pi_s("N_pi_s","N_pi_s",0.9*ent, 0.,1.05*ent) ;
	 RooRealVar N_pi_b("N_pi_b","N_pi_b",0.1*ent, 0.,1.05*ent) ;

	 ent = invMassHist[1][pBinIndex][thetaindex]->GetEntries();
	 RooRealVar N_k_s("N_k_s","N_k_s",0.3*ent,0.,1.05*ent) ;
	 RooRealVar N_k_b("N_k_b","N_k_b",0.7*ent,0.,1.05*ent) ;

	 ent = invMassHist[2][pBinIndex][thetaindex]->GetEntries();
	 RooRealVar N_p_s("N_p_s","N_p_s",0.2*ent,0.,1.05*ent) ;
	 RooRealVar N_p_b("N_p_b","N_p_b",0.8*ent,0.,1.05*ent) ;

	 ent = invMassHist[3][pBinIndex][thetaindex]->GetEntries();
	 RooRealVar N_u_s("N_u_s","N_u_s",0.8*ent,0.,1.05*ent) ;
	 RooRealVar N_u_b("N_u_b","N_u_b",0.05*ent,0.,1.05*ent) ;

	 ent = invMassHist[5][pBinIndex][thetaindex]->GetEntries();
         RooRealVar N_unk_s("N_unk_s","N_unk_s",0.8*ent,0.,1.05*ent) ;
         RooRealVar N_unk_b("N_unk_b","N_unk_b",0.05*ent,0.,1.05*ent) ;
	 
	 RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s + N_unk_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s,N_unk_s));

       RooAddPdf model_all("model_all","model_all",RooArgList(signal,background),RooArgList(N_a_s,N_a_b)) ;
       RooAddPdf model_pi("model_pi","model_pi",RooArgList(signal, background),RooArgList(N_pi_s,N_pi_b)) ; // checked with bgnpi
       RooAddPdf model_k("model_k","model_k",RooArgList(signal, background),RooArgList(N_k_s,N_k_b)) ; // checked with bgnk
       RooAddPdf model_p("model_p","model_p",RooArgList(signal, background),RooArgList(N_p_s,N_p_b)) ;
       RooAddPdf model_bkg("model_bkg","model_bkg",RooArgList(signal, background),RooArgList(N_u_s,N_u_b)) ;
       RooAddPdf model_unk("model_unk","model_unk",RooArgList(signal, background),RooArgList(N_unk_s,N_unk_b)) ;

       
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////       

       
       //RooRealVar N_signal("N_signal", "Signal Events", 1000000, 0, 1000000);
       //RooRealVar N_background("N_background", "Background Events", 1000, 0, 10000);

       // Combine signal and background for individual fitting
       //RooAddPdf model("model", "Signal + Background", RooArgList(signal, background), RooArgList(N_signal, N_background));


       RooCategory category("category", "Category");
       category.defineType("Pi_as_Pi");
       category.defineType("Pi_as_K");
       category.defineType("Pi_as_Proton");
       category.defineType("Pi_as_bkg");
       category.defineType("Pi_not_identified");
       category.defineType("All_Events");

       RooDataHist data_pi_as_pi("data_pi_as_pi", "Pi+ as Pi+", mass, RooFit::Import(*invMassHist[0][pBinIndex][thetaindex]));
       RooDataHist data_pi_as_k("data_pi_as_k", "Pi+ as K+", mass, RooFit::Import(*invMassHist[1][pBinIndex][thetaindex]));
       RooDataHist data_pi_as_proton("data_pi_as_proton", "Pi+ as Proton", mass, RooFit::Import(*invMassHist[2][pBinIndex][thetaindex]));
       RooDataHist data_pi_as_bkg("data_pi_as_bkg", "Pi+ as bkg", mass, RooFit::Import(*invMassHist[3][pBinIndex][thetaindex]));
       RooDataHist data_pi_not_identified("data_pi_not_identified", "Pi+ not identified", mass, RooFit::Import(*invMassHist[5][pBinIndex][thetaindex]));
       RooDataHist data_all_events("data_all_events", "All Events", mass, RooFit::Import(*invMassHist[4][pBinIndex][thetaindex]));

       // Combine the datasets using the category
       RooDataHist combinedData("combinedData", "Combined Data", RooArgSet(mass), RooFit::Index(category),
				RooFit::Import("Pi_as_Pi", data_pi_as_pi),
				RooFit::Import("Pi_as_K", data_pi_as_k),
				RooFit::Import("Pi_as_Proton", data_pi_as_proton),
				RooFit::Import("Pi_as_bkg", data_pi_as_bkg),
				RooFit::Import("Pi_not_identified", data_pi_not_identified),
				RooFit::Import("All_Events", data_all_events));

       // Define a simultaneous model using the category
       RooSimultaneous simPdf("simPdf", "Simultaneous Fit", category);
       simPdf.addPdf(model_pi, "Pi_as_Pi");
       simPdf.addPdf(model_k, "Pi_as_K");
       simPdf.addPdf(model_p, "Pi_as_Proton");
       simPdf.addPdf(model_bkg, "Pi_as_bkg");
       simPdf.addPdf(model_unk, "Pi_not_identified");
       simPdf.addPdf(model_all, "All_Events");

       RooFormulaVar restriction("restriction", "100000 * (TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s + N_unk_s) / N_a_s) > 1e-4)",
				 RooArgSet(N_a_s, N_pi_s, N_k_s, N_p_s, N_u_s, N_unk_s));

       int iter = 0; 
       do {
	 // Adjust parameters slightly if reattempting the fit
	 if (iter > 0) {
	   N_pi_s.setVal(1.01 * N_pi_s.getVal());
	   N_k_s.setVal(0.98 * N_k_s.getVal());
	   N_p_s.setVal(0.98 * N_p_s.getVal());
	   N_u_s.setVal(0.98 * N_u_s.getVal());
	   N_unk_s.setVal(0.98 * N_unk_s.getVal());
	 }


       
       // Create a simple negative log-likelihood for a PDF (assuming simPdf and combData are defined)
       RooAbsReal* nll = simPdf.createNLL(combinedData, RooFit::Extended(true));
       
       // Add the restriction to the NLL to penalize mismatches in counts
       RooAddition nll_r("nll_r", "nll_r", RooArgSet(*nll, restriction));
       
       // Set up the minimizer and run the fit
       RooMinimizer minimizer(*nll); // Create a RooMinimizer instance
       minimizer.setPrintLevel(0); // Set print level for minimizer
       minimizer.migrad(); // Run the minimization
       minimizer.hesse();             // Compute covariance matrix initially
       minimizer.simplex();           // Run the Simplex algorithm for initial optimization
       int ret = minimizer.migrad();  // Main minimization step using MIGRAD
       // Retry MIGRAD if initial attempt fails
       if (ret) {
	 ret = minimizer.migrad();
	 if (!ret) printf("=== Successful retry\n");
       }

       // Improve the fit further and re-evaluate covariance matrix
       if (!ret) minimizer.improve();
       minimizer.hesse();

       // Save fit results for this iteration
       fit_results[p_index][pBinIndex][thetaindex] = minimizer.save();  // Store fit results for current conditions
       
       iter++;
       
       } while (fit_results[p_index][pBinIndex][thetaindex]->covQual() != 3 && iter <= 100);
       
       printf("=> %d iterations -> Covariance Quality: %d\n", iter, fit_results[p_index][pBinIndex][thetaindex]->covQual());

       std::unordered_set<int> relindex = {2, 4, 6, 8, 10};  // Use unordered_set for fast lookup
       
       if (fit_results[p_index][pBinIndex][thetaindex] != nullptr) { 
	 covdata << " Covariant Matrix for ["<<p_index<<","<<pBinIndex<<", "<<thetaindex<<"] \n";
	 TMatrixDSym covMatrix = fit_results[p_index][pBinIndex][thetaindex]->covarianceMatrix();

	 for (int i = 0; i < covMatrix.GetNrows(); ++i) {
	   for (int j = 0; j < covMatrix.GetNcols(); ++j) {
	     covdata<<std::setw(12)<<covMatrix(i,j)<<" " ;
	     // Check if both i and j are in relindex and if i >= j
	     if (relindex.find(i) != relindex.end() && relindex.find(j) != relindex.end() && (i >= j)) {
	      double value = covMatrix(i, j);
	      covvalues<<p_index<<"\t"<<(mombins.at(pBinIndex) + mombins.at(pBinIndex + 1))/2<<"\t"<<(thetabins.at(thetaindex) + thetabins.at(thetaindex +1))/2<<"\t"<<N_a_s.getVal() <<"\t\t"<<N_pi_s.getVal()<<"\t\t"<<N_k_s.getVal()<<"\t\t"<<N_p_s.getVal()<<"\t\t"<<N_u_s.getVal()<<"\t\t"<<N_unk_s.getVal()<<"\t\t"<<i << "\t" << j << "\t" << value << "\t" <<Q2<<"\n";
	     }
	   }

	   covdata<< std::endl;
	 }

	 covdata<< std::endl;
       }
       
     



     
       outputtext << p_index<<"\t\t"<<mombins.at(pBinIndex) << "\t\t"<<thetabins.at(thetaindex)<<"\t\t"<<N_a_s.getVal() <<"\t\t"<<N_pi_s.getVal()<<"\t\t"<<N_k_s.getVal()<<"\t\t"<<N_p_s.getVal()<<"\t\t"<<N_u_s.getVal()<<"\t\t"<<N_unk_s.getVal()<<"\n";
       
     
       // List of category names and titles
       std::vector<std::string> categories = {"Pi_as_Pi", "Pi_as_K", "Pi_as_Proton", "Pi_as_bkg", "Pi_not_identified", "All_Events"};
       std::vector<std::string> titles = {"Pi+ as Pi+", "Pi+ as K+", "Pi+ as Proton", "Pi+ as bkg", "Pi+ not identified", "All Events"};
       
       RooPlot* frame = mass.frame();
       combinedData.plotOn(frame, RooFit::Cut(Form("category==category::%s", categories[p_index].c_str())));

       // Plot the overall model (signal + background)
       simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData),
		     RooFit::Name("Total Fit"));

       // Plot the signal component only
       simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData),
		     RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::Name("Signal"));

	   // Plot the background component only
       simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData),
		     RooFit::Components("background"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed), RooFit::Name("Background"));
       simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData));
       
       frame->SetTitle(titles[p_index].c_str());
       frame->Draw();
       
       TLatex latex;
       latex.SetNDC();
       latex.SetTextSize(0.03);
       latex.DrawLatex(0.6, 0.85, Form("Chi2/ndf: %.2f", frame->chiSquare()));
       latex.DrawLatex(0.6, 0.8, Form(" Momentum >  %.1f", mombins.at(pBinIndex)));
       latex.DrawLatex(0.6, 0.75, Form(" Theta > %.2f", thetabins.at(thetaindex)));
       
       c1->Print("kmDST_ang_momfit_results_pr_sim.pdf");

       delete frame;
       }
     }
   }
  
   c1->Print("kmDST_ang_momfit_results_pr_sim.pdf]"); 
   delete c1;
   outputtext.close();
   covdata.close();
   covvalues.close();
}
