#define fit_k_cxx
#include "fit_k.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream> // Include fstream for file I/O

double m_pi = 0.139570;

void fit_k::Loop()
{

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH2D* h_arment = new TH2D("arment","arment", 2000, -1, 1, 400, 0, 0.4);
   TH1D* h_invmass_k = new TH1D("kmass","kmass",200,0.4,0.6);

   //momentum bins: Momentum p(GeV/c) = (10,11,12,13,15,17,19,22,25,27,30,35,40,50)

   std::vector<float> mombins = {10.0 ,11.0 ,12.0 ,13.0 ,15.0 ,17.0 ,19.0 ,22.0 ,25.0 ,27.0 ,30.0 ,35.0 ,40.0 ,50.0 };
   size_t numMomentumBins = mombins.size();
   
   std::vector<std::vector<TH1F*>> invMassHist(5, std::vector<TH1F*>(numMomentumBins - 1, nullptr));
   
   //TH1F* invMassHist[5][13];
   for (int p = 0; p < 5; ++p) {
     for (size_t i = 0; i < mombins.size(); ++i) {
       double min_p = mombins.at(i);
       TString histName = Form("invMass_%d_for_mom_%.2f", p, min_p);
       invMassHist[p][i] = new TH1F(histName, Form("Invariant Mass for p_index = %d for mom %.2f ", p, min_p), 200, 0.4, 0.6);
     }
   }

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
	invMassHist[p_index][pBinIndex]->Fill(k_mass);
	invMassHist[4][pBinIndex]->Fill(k_mass);
	h_invmass_k->Fill(k_mass);
	h_arment->Fill(alpha, pt1);

      }


      
   }

   TFile *outputFile = new TFile("k0_momfit_input.root", "RECREATE");

   h_arment->Write();
   h_invmass_k->Write();

   for (int p = 0; p < 5; ++p) {
     for(int i = 0; i < 13; i++){
	
       invMassHist[p][i]->Write();
     }
   }
    
   outputFile->Close();
   
   RooFitResult* fit_results[5][13] = {nullptr};
   
   std::ofstream covdata("k_covmatrix.txt");

   std::ofstream covvalues("k_covvalues.txt");
   covvalues<<"Category\tMomentum\tN_as\tN_pi_s\tN_k_s\tN_p_s\tN_u_s\tindex1\tindex2\tvalue\n";
   
   std::ofstream outputtext("k_signal_background_results_simfit.txt");

   outputtext << "Category \t Momentum Bin \t N_as \t\t N_pi_s \t N_k_s \t\t N_p_s \t\t N_u_s \n";
  
   TCanvas* c1 = new TCanvas("c1", "Fit Results", 800, 600);
   c1->Print("k_momfit_results_sim.pdf[");

   for (int p_index = 0; p_index < 5; ++p_index) { 

     for (int pBinIndex = 0; pBinIndex < mombins.size() - 1; ++pBinIndex) {
      
       
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
       int ent = invMassHist[4][pBinIndex]->GetEntries();
       RooRealVar N_a_b("N_a_b", "N_a_b", 0.05*ent, 0., 1.05*ent);

       ent = invMassHist[0][pBinIndex]->GetEntries();
       RooRealVar N_pi_s("N_pi_s","N_pi_s",0.9*ent, 0.,1.05*ent) ;
       RooRealVar N_pi_b("N_pi_b","N_pi_b",0.1*ent, 0.,1.05*ent) ;

       ent = invMassHist[1][pBinIndex]->GetEntries();
       RooRealVar N_k_s("N_k_s","N_k_s",0.3*ent,0.,1.05*ent) ;
       RooRealVar N_k_b("N_k_b","N_k_b",0.7*ent,0.,1.05*ent) ;

       ent = invMassHist[2][pBinIndex]->GetEntries();
       RooRealVar N_p_s("N_p_s","N_p_s",0.2*ent,0.,1.05*ent) ;
       RooRealVar N_p_b("N_p_b","N_p_b",0.8*ent,0.,1.05*ent) ;

       ent = invMassHist[3][pBinIndex]->GetEntries();
       RooRealVar N_u_s("N_u_s","N_u_s",0.8*ent,0.,1.05*ent) ;
       RooRealVar N_u_b("N_u_b","N_u_b",0.05*ent,0.,1.05*ent) ;

       RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));

       RooAddPdf model_all("model_all","model_all",RooArgList(signal,background),RooArgList(N_a_s,N_a_b)) ;
       RooAddPdf model_pi("model_pi","model_pi",RooArgList(signal, background),RooArgList(N_pi_s,N_pi_b)) ; // checked with bgnpi
       RooAddPdf model_k("model_k","model_k",RooArgList(signal, background),RooArgList(N_k_s,N_k_b)) ; // checked with bgnk
       RooAddPdf model_p("model_p","model_p",RooArgList(signal, background),RooArgList(N_p_s,N_p_b)) ;
       RooAddPdf model_unk("model_unk","model_unk",RooArgList(signal, background),RooArgList(N_u_s,N_u_b)) ;

       
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////       

       
       //RooRealVar N_signal("N_signal", "Signal Events", 1000000, 0, 1000000);
       //RooRealVar N_background("N_background", "Background Events", 1000, 0, 10000);

       // Combine signal and background for individual fitting
       //RooAddPdf model("model", "Signal + Background", RooArgList(signal, background), RooArgList(N_signal, N_background));


       RooCategory category("category", "Category");
       category.defineType("Pi_as_Pi");
       category.defineType("Pi_as_K");
       category.defineType("Pi_as_Proton");
       category.defineType("Pi_not_identified");
       category.defineType("All_Events");

       RooDataHist data_pi_as_pi("data_pi_as_pi", "Pi+ as Pi+", mass, RooFit::Import(*invMassHist[0][pBinIndex]));
       RooDataHist data_pi_as_k("data_pi_as_k", "Pi+ as K+", mass, RooFit::Import(*invMassHist[1][pBinIndex]));
       RooDataHist data_pi_as_proton("data_pi_as_proton", "Pi+ as Proton", mass, RooFit::Import(*invMassHist[2][pBinIndex]));
       RooDataHist data_pi_not_identified("data_pi_not_identified", "Pi+ not identified", mass, RooFit::Import(*invMassHist[3][pBinIndex]));
       RooDataHist data_all_events("data_all_events", "All Events", mass, RooFit::Import(*invMassHist[4][pBinIndex]));

       // Combine the datasets using the category
       RooDataHist combinedData("combinedData", "Combined Data", RooArgSet(mass), RooFit::Index(category),
				RooFit::Import("Pi_as_Pi", data_pi_as_pi),
				RooFit::Import("Pi_as_K", data_pi_as_k),
				RooFit::Import("Pi_as_Proton", data_pi_as_proton),
				RooFit::Import("Pi_not_identified", data_pi_not_identified),
				RooFit::Import("All_Events", data_all_events));

       // Define a simultaneous model using the category
       RooSimultaneous simPdf("simPdf", "Simultaneous Fit", category);
       simPdf.addPdf(model_pi, "Pi_as_Pi");
       simPdf.addPdf(model_k, "Pi_as_K");
       simPdf.addPdf(model_p, "Pi_as_Proton");
       simPdf.addPdf(model_unk, "Pi_not_identified");
       simPdf.addPdf(model_all, "All_Events");

       RooFormulaVar restriction("restriction", "100000 * (TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s) / N_a_s) > 1e-4)",
				  RooArgSet(N_a_s, N_pi_s, N_k_s, N_p_s, N_u_s));

       int iter = 0; 
       do {
	 // Adjust parameters slightly if reattempting the fit
	 if (iter > 0) {
	   N_pi_s.setVal(1.01 * N_pi_s.getVal());
	   N_k_s.setVal(0.98 * N_k_s.getVal());
	   N_p_s.setVal(0.98 * N_p_s.getVal());
	   N_u_s.setVal(0.98 * N_u_s.getVal());
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
       fit_results[p_index][pBinIndex] = minimizer.save();  // Store fit results for current conditions
       
       iter++;
       
       } while (fit_results[p_index][pBinIndex]->covQual() != 3 && iter <= 100);
       
       printf("=> %d iterations -> Covariance Quality: %d\n", iter, fit_results[p_index][pBinIndex]->covQual());

       std::unordered_set<int> relindex = {2, 4, 6, 8};  // Use unordered_set for fast lookup
       
       if (fit_results[p_index][pBinIndex] != nullptr) { 
	 covdata << " Covariant Matrix for ["<<p_index<<","<<pBinIndex<<"] \n";
	 TMatrixDSym covMatrix = fit_results[p_index][pBinIndex]->covarianceMatrix();

	 for (int i = 0; i < covMatrix.GetNrows(); ++i) {
	   for (int j = 0; j < covMatrix.GetNcols(); ++j) {
	     covdata<<std::setw(12)<<covMatrix(i,j)<<" " ;
	     // Check if both i and j are in relindex and if i >= j
	     if (relindex.find(i) != relindex.end() && relindex.find(j) != relindex.end() && (i >= j)) {
	      double value = covMatrix(i, j);
	      covvalues<<p_index<<"\t"<<mombins.at(pBinIndex)<<"\t"<<N_a_s.getVal() <<"\t\t"<<N_pi_s.getVal()<<"\t\t"<<N_k_s.getVal()<<"\t\t"<<N_p_s.getVal()<<"\t\t"<<N_u_s.getVal()<<"\t\t"<<i << "\t" << j << "\t" << value << "\n";
	     }
	   }

	   covdata<< std::endl;
	 }

	 covdata<< std::endl;
       }
       
	// Save the fit result
       //RooFitResult* result = minimizer.save(); // Save the fit result
       
       // Fit the model to the data
       //RooFitResult* fitResult = simPdf.fitTo(combinedData, RooFit::Save());

       // Write the results to the text file
       /* outputtext <<" "  << "\t\t" << mombins.at(pBinIndex) << "\t\t" 
	  << N_signal.getVal() << "\t\t" << N_signal.getError() << "\t\t" 
	  << N_background.getVal() << "\t\t" << N_background.getError() << "\n";
       */

       outputtext << p_index<<"\t\t"<<mombins.at(pBinIndex) << "\t\t"<<N_a_s.getVal() <<"\t\t"<<N_pi_s.getVal()<<"\t\t"<<N_k_s.getVal()<<"\t\t"<<N_p_s.getVal()<<"\t\t"<<N_u_s.getVal()<<"\n";
       

       // List of category names and titles
       std::vector<std::string> categories = {"Pi_as_Pi", "Pi_as_K", "Pi_as_Proton", "Pi_not_identified", "All_Events"};
       std::vector<std::string> titles = {"Pi+ as Pi+", "Pi+ as K+", "Pi+ as Proton", "Pi+ not identified", "All Events"};

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
       latex.DrawLatex(0.6, 0.8, Form(" %.1f < Momentum < %.1f", mombins.at(pBinIndex), mombins.at(pBinIndex + 1)));

       c1->Print("k_momfit_results_sim.pdf");

       delete frame;
     }
   }
   
   c1->Print("k_momfit_results_sim.pdf]"); 
   delete c1;
   outputtext.close();
   covdata.close();
   covvalues.close();
}

