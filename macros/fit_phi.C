#define fit_phi_cxx
#include "fit_phi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream> // Include fstream for file I/O

double m_pi = 0.139570;

void fit_phi::Loop()
{

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH2D* h_arment = new TH2D("arment","arment", 2000, -1, 1, 400, 0, 0.4);
   
   TH1D* h_invmass_phi = new TH1D("phimass","phimass",200,0.98,1.11);

   //momentum bins: Momentum p(GeV/c) = (10,11,12,13,15,17,19,22,25,27,30,35,40,50)

   std::vector<float> mombins = {10.0 ,11.0 ,12.0 ,13.0 ,15.0 ,17.0 ,19.0 ,22.0 ,25.0 ,27.0 ,30.0 ,35.0 ,40.0 ,50.0 };
   size_t numMomentumBins = mombins.size();
   
   std::vector<std::vector<TH1F*>> invMassHist(5, std::vector<TH1F*>(numMomentumBins - 1, nullptr));
   
   //TH1F* invMassHist[5][13];
   for (int p = 0; p < 5; ++p) {
     for (size_t i = 0; i < mombins.size(); ++i) {
       double min_p = mombins.at(i);
       TString histName = Form("invMass_%d_for_mom_%.2f", p, min_p);
       invMassHist[p][i] = new TH1F(histName, Form("Invariant Mass for p_index = %d for mom %.2f ", p, min_p), 200, 0.98, 1.11);
     }
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double n_rich = sqrt((m_pi*m_pi/(pi_thr*pi_thr)) + 1);
      if(pp_mom < k_thr || pm_mom < k_thr){ continue ; }


      // if( pt1 < 0.10 || pt2 < 0.10){ continue ;}

      //if(pp_mom < 10. || pm_mom < 10.){ continue ; }
      //if( pt1 < 0.1 | pt2 < 0.1 ){ continue ; }

      float theta_pp = pp_theta;
      float theta_pm = pm_theta;

      if(theta_pp == 0 || theta_pp == -1){ continue ; }
      if(theta_pm == 0 || theta_pm == -1){ continue ; }

      double phi_mass = lv_phi->M();
      

      //      if( abs(k_mass - 0.497) > 0.020 ) { continue ; }
      //if( lv_k0->Pt() < 0.040 ) { continue ; }

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
	
	/* if(lh_values_m[m] > m_maxvalue) {
          m_maxvalue = lh_values_m[m];
          m_index = m;
	  }*/
      }

      //new selection for indices for negative particle

      // Find m_maxvalue and m_index with new conditions
      for(int m = 0; m < 4; m++) {
	bool condition_met = true;
	
	// Check conditions for j = 0 and 2 , i.e, for pion and proton
	for(int j : {0, 2}) {
	  if (lh_values_m[m] / lh_values_m[j] <= 1.08) {
            condition_met = false;
            break;
	  }
	}
	
	// Check condition for j = 3 , i.e, background
	if (lh_values_m[m] / lh_values_m[3] <= 1.24) {
	  condition_met = false;
	}

	//new condition that (LH pi - LH bkg) > 0.1

	if( abs(lh_values_m[0] - lh_values_m[3]) < 0.1){
	  condition_met = false;
	}
	
	// Update m_maxvalue and m_index if conditions are met and value is greater
	if(condition_met && lh_values_m[m] > m_maxvalue) {
	  m_maxvalue = lh_values_m[m];
	  m_index = m;
	}
      }

      
   
      size_t pBinIndex = 0;
      while (pBinIndex < mombins.size() - 1 && !(pp_mom >= mombins[pBinIndex] && pp_mom < mombins[pBinIndex + 1])) {
        ++pBinIndex;
      }

      if( m_index == 1  ){
	invMassHist[p_index][pBinIndex]->Fill(phi_mass);
        invMassHist[4][pBinIndex]->Fill(phi_mass);
	h_arment->Fill(alpha, pt1);
	h_invmass_phi->Fill(phi_mass);
      }

   }

   TFile *outputFile = new TFile("phi_momfit_input_lh3.root", "RECREATE");

   h_arment->Write();
   h_invmass_phi->Write();

   for (int p = 0; p < 5; ++p) {
     for(int i = 0; i < 13; i++){
	
       invMassHist[p][i]->Write();
     }
   }
    
   outputFile->Close();
   
   RooFitResult* fit_results[5][13] = {nullptr};
   
   std::ofstream covdata("phi_covmatrix_lh3.txt");

   std::ofstream covvalues("phi_covvalues_lh3.txt");
   covvalues<<"Category\tMomentum\tN_as\tN_pi_s\tN_k_s\tN_p_s\tN_u_s\tindex1\tindex2\tvalue\n";
   
   std::ofstream outputtext("phi_signal_background_results_simfit_lh3.txt");

   outputtext << "Category \t Momentum Bin \t N_as \t\t N_pi_s \t N_k_s \t\t N_p_s \t\t N_u_s \n";
  
   TCanvas* c1 = new TCanvas("c1", "Fit Results", 800, 600);
   c1->Print("phi_momfit_results_sim_lh3.pdf[");

   for (int p_index = 0; p_index < 5; ++p_index) { 

     for (int pBinIndex = 0; pBinIndex < mombins.size() - 1; ++pBinIndex) {

       //Signal
       RooRealVar x("x","M",0.995,1.042,"GeV");

       RooRealVar mass("mass","mass #phi",1.01945,"GeV") ;
       RooRealVar width("width","width #phi",0.00426,"GeV");
       RooRealVar mean("mean","mean",0,-.04,.04,"GeV") ;
       RooRealVar sigma("sigma","sigma",0.0031,0.0025,0.005,"GeV");
       RooRealVar d1("d1","d1",1.);
       RooRealVar d2("d2","d2",0.4937);

       //RooRelBreitWigner sb("sb","relBW",x,mass,width,d1,d2,d2,d1);
       RooBreitWigner sb("sb", "relBW", x, mass, width);
       RooGaussian sg("sg","gauss",x,mean,sigma);


       RooFFTConvPdf sig("sig","sig",x,sb,sg);
       //sig.setShift(0,0);

       //Background
       RooRealVar thr("thr","threshold",0.987354,0.985,0.99,"GeV");
       RooRealVar thr2("thr2","threshold",0.987354,0.985,0.99,"GeV");
       RooRealVar thr3("thr3","threshold",0.987354,0.985,0.99,"GeV");
       RooRealVar thr4("thr4","threshold",0.987354,0.985,0.99,"GeV");
       RooRealVar thr5("thr5","threshold",0.987354,0.985,0.99,"GeV");

       RooRealVar n("n","n",  0.4,  0.,   1.5) ;
       RooRealVar a("a","a",  4.,   0.,    100.) ;
       RooGenericPdf bgn1("bgn1","x<thr ? 0 :TMath::Power(x - thr,n)*TMath::Exp(-a*(x - thr))",RooArgList(x,thr,n,a));

       RooRealVar n2("n2","n2",  0.4,  0.,   1.5) ;
       RooRealVar a2("a2","a2",  4.,   0.,    100.) ;
       RooGenericPdf bgn2("bgn2","x<thr2 ? 0 :TMath::Power(x - thr2,n2)*TMath::Exp(-a2*(x - thr2))",RooArgList(x,thr2,n2,a2));

       RooRealVar n3("n3","n3",  0.4,  0.,   1.5) ;
       RooRealVar a3("a3","a3",  4.,   0.,    100.) ;
       RooGenericPdf bgn3("bgn3","x<thr3 ? 0 :TMath::Power(x - thr3,n3)*TMath::Exp(-a3*(x - thr3))",RooArgList(x,thr3,n3,a3));

       RooRealVar n4("n4","n4",  0.4,  0.,   1.5) ;
       RooRealVar a4("a4","a4",  4.,   0.,    100.) ;
       RooGenericPdf bgn4("bgn4","x<thr4 ? 0 :TMath::Power(x - thr4,n4)*TMath::Exp(-a4*(x - thr4))",RooArgList(x,thr4,n4,a4));

       RooRealVar n5("n5","n5",  0.4,  0.,   1.5) ;
       RooRealVar a5("a5","a5",  4.,   0.,    100.) ;
       RooGenericPdf bgn5("bgn5","x<thr5 ? 0 :TMath::Power(x - thr5,n5)*TMath::Exp(-a5*(x - thr5))",RooArgList(x,thr5,n5,a5));
       
       // extra part for error using covariant matrix

       // Number of signal and background events for each category
       int ent = invMassHist[4][pBinIndex]->GetEntries();
       RooRealVar N_a_b("N_a_b", "N_a_b", 0.2*ent, 0., 1.05*ent);

       ent = invMassHist[0][pBinIndex]->GetEntries();
       RooRealVar N_pi_s("N_pi_s","N_pi_s",0.12*ent,0.,1.05*ent) ;
       RooRealVar N_pi_b("N_pi_b","N_pi_b",0.5*ent,0.,1.05*ent) ;
       RooRealVar N_pi_b2("N_pi_b2","N_pi_b2",0.38*ent,0.,1.05*ent) ;

       ent = invMassHist[1][pBinIndex]->GetEntries();
       RooRealVar N_k_s("N_k_s","N_k_s",0.87*ent,0.,1.05*ent) ;
       RooRealVar N_k_b("N_k_b","N_k_b",0.13*ent,0.,1.05*ent) ;

       ent = invMassHist[2][pBinIndex]->GetEntries();
       RooRealVar N_p_s("N_p_s","N_p_s",0.87*ent,0.,1.05*ent) ;
       RooRealVar N_p_b("N_p_b","N_p_b",0.13*ent,0.,1.05*ent) ;

       ent = invMassHist[3][pBinIndex]->GetEntries();
       RooRealVar N_u_s("N_u_s","N_u_s",0.77*ent,0.,1.05*ent) ;
       RooRealVar N_u_b("N_u_b","N_u_b",0.23*ent,0.,1.05*ent) ;
       RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));

       RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgn1),RooArgList(N_a_s,N_a_b)) ;
       RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgn2),RooArgList(N_pi_s,N_pi_b)) ;
       RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgn3),RooArgList(N_k_s,N_k_b)) ;
       RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgn1),RooArgList(N_p_s,N_p_b)) ;
       RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgn1),RooArgList(N_u_s,N_u_b)) ;

       x.setBins(30);

       
       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////       

       
       //RooRealVar N_signal("N_signal", "Signal Events", 1000000, 0, 1000000);
       //RooRealVar N_background("N_background", "Background Events", 1000, 0, 10000);

       // Combine signal and background for individual fitting
       //RooAddPdf model("model", "Signal + Background", RooArgList(signal, background), RooArgList(N_signal, N_background));


       RooCategory category("category", "Category");
       category.defineType("K_as_Pi");
       category.defineType("K_as_K");
       category.defineType("K_as_Proton");
       category.defineType("K_not_identified");
       category.defineType("All_Events");

       RooDataHist data_k_as_pi("data_k_as_pi", "K+ as Pi+", x, RooFit::Import(*invMassHist[0][pBinIndex]));
       RooDataHist data_k_as_k("data_k_as_k", "K+ as K+", x, RooFit::Import(*invMassHist[1][pBinIndex]));
       RooDataHist data_k_as_proton("data_k_as_proton", "K+ as Proton", x, RooFit::Import(*invMassHist[2][pBinIndex]));
       RooDataHist data_k_not_identified("data_k_not_identified", "K+ not identified", x, RooFit::Import(*invMassHist[3][pBinIndex]));
       RooDataHist data_all_events("data_all_events", "All Events", x, RooFit::Import(*invMassHist[4][pBinIndex]));

       // Combine the datasets using the category
       RooDataHist combinedData("combinedData", "Combined Data", RooArgSet(x), RooFit::Index(category),
				RooFit::Import("K_as_Pi", data_k_as_pi),
				RooFit::Import("K_as_K", data_k_as_k),
				RooFit::Import("K_as_Proton", data_k_as_proton),
				RooFit::Import("K_not_identified", data_k_not_identified),
				RooFit::Import("All_Events", data_all_events));

       // Define a simultaneous model using the category
       RooSimultaneous simPdf("simPdf", "Simultaneous Fit", category);
       simPdf.addPdf(model_pi, "K_as_Pi");
       simPdf.addPdf(model_k, "K_as_K");
       simPdf.addPdf(model_p, "K_as_Proton");
       simPdf.addPdf(model_unk, "K_not_identified");
       simPdf.addPdf(model_all, "All_Events");

       	RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));

		
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
       
       } while (fit_results[p_index][pBinIndex]->covQual() != 3 && iter <= 5);
       
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
	       covvalues<<p_index<<"\t"<<mombins.at(pBinIndex)<<"\t"<<N_a_s.getVal() <<"\t\t"<<N_pi_s.getVal()<<"\t\t"<<N_k_s.getVal \
		 ()<<"\t\t"<<N_p_s.getVal()<<"\t\t"<<N_u_s.getVal()<<"\t\t"<<i << "\t" << j << "\t" << value << "\n";
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
       std::vector<std::string> categories = {"K_as_Pi", "K_as_K", "K_as_Proton", "K_not_identified", "All_Events"};
       std::vector<std::string> titles = {"K+ as Pi+", "K+ as K+", "K+ as Proton", "K+ not identified", "All Events"};
       std::map<std::string, std::string> backgroundComponents = {
	 {"K_as_Pi", "bgn2"},
	 {"K_as_K", "bgn3"},
	 {"K_as_Proton", "bgn1"},
	 {"K_not_identified", "bgn1"},
	 {"All_Events", "bgn1"}
       };
     
       //       for (size_t i = 0; i < categories.size(); ++i) {
	 RooPlot* frame = x.frame();  
	 combinedData.plotOn(frame, RooFit::Cut(Form("category==category::%s", categories[p_index].c_str())));

	 // Plot the overall model (signal + background)
	 simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData),
		       RooFit::Name("Total Fit"));

	 // Plot the signal component only
	 simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData),
		       RooFit::Components("sig"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::Name("Signal"));

	 // Plot the specific background component based on the category
	 simPdf.plotOn(frame, RooFit::Slice(category, categories[p_index].c_str()), RooFit::ProjWData(category, combinedData),
		       RooFit::Components(backgroundComponents[categories[p_index]].c_str()), RooFit::LineColor(kGreen),
		       RooFit::LineStyle(kDashed), RooFit::Name("Background"));

	 frame->SetTitle(titles[p_index].c_str());
	 frame->Draw();

	 TLatex latex;
	 latex.SetNDC();
	 latex.SetTextSize(0.03);
	 latex.DrawLatex(0.6, 0.85, Form("Chi2/ndf: %.2f", frame->chiSquare()));
	 latex.DrawLatex(0.6, 0.8, Form(" %.1f < Momentum < %.1f", mombins.at(pBinIndex), mombins.at(pBinIndex + 1)));

	 c1->Print("phi_momfit_results_sim_lh3.pdf");
	 delete frame;
	 //}
      
     }
   }
   
   c1->Print("phi_momfit_results_sim_lh3.pdf]"); 
   delete c1;
   outputtext.close();
   covdata.close();
   covvalues.close();
}

