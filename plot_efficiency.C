#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

void partials(const std::vector<double>& Ni, const std::vector<double>& Na, std::vector<double>& dNi) {
  if (Ni.size() != Na.size()) {
    throw std::invalid_argument("Vectors Ni and Na must have the same size.");
  }
  dNi.resize(Ni.size());
  
  for (size_t i = 0; i < Ni.size(); ++i) {
    dNi[i] = (Na[i] - Ni[i]) / (Na[i] * Na[i]);
  }
}

double var_term(const std::vector<double>& dNi, const std::vector<double>& v) {  
    double a = 0.0;
    for (size_t i = 0; i < dNi.size(); ++i) {
        a += std::pow(dNi[i], 2) * v[i];
    }
    return a;
}

double covar_term(const std::vector<double>& dNi, const std::vector<std::vector<double>>& cv) {
    double a = 0.0;

    for (size_t i = 0; i < dNi.size() - 1; ++i) {
        for (size_t j = i + 1; j < dNi.size(); ++j) {
            a += dNi[i] * dNi[j] * cv[i][j];
        }
    }
    return a;
}

void create_covariance_matrix(const std::vector<size_t>& index1, const std::vector<size_t>& index2, const std::vector<double>& value, 
                              std::vector<std::vector<double>>& cv, size_t max_index) {

  cv.resize(max_index + 1, std::vector<double>(max_index + 1, 0.0));

  if (index1.size() != index2.size() || index1.size() != value.size()) {
    throw std::invalid_argument("Vectors index1, index2, and value must have the same size.");
  }

  for (size_t i = 0; i < index1.size(); ++i) {
    if (index1[i] <= max_index && index2[i] <= max_index) {
      cv[index1[i]][index2[i]] = value[i];
      cv[index2[i]][index1[i]] = value[i]; // Ensure symmetry
    }
  }
}


void plot_efficiency() {

  std::ifstream infile("k_covvalues.txt");

    int category;
    double momentum, N_as, N_pi_s, N_k_s, N_p_s, N_u_s;
    int index1, index2;
    double value;
    
    std::map<int, std::vector<std::pair<double, double>>> efficiencies_by_category;

    std::string line;
    std::getline(infile, line);


    while (std::getline(infile, line)) {
        std::stringstream ss(line);

	if (line.empty()) continue;

	ss >> category >> momentum >> N_as >> N_pi_s >> N_k_s >> N_p_s >> N_u_s >> index1 >> index2 >> value;

	double efficiency = 0.0;
	if(category == 0){
	  efficiency = N_pi_s/N_as;
	}
	else if (category == 1){
	  efficiency = N_k_s/N_as;
	}
	else if (category == 2){
	  efficiency = N_p_s/N_as;
	}
	else if ( category == 3){
	  efficiency = N_u_s/N_as;
	}
	else if (category == 4){
	  efficiency = 1.0 ;
	}
	efficiencies_by_category[category].push_back(std::make_pair(momentum, efficiency));
    }

	
    // Close the input file
    infile.close();

    // Print efficiencies for debugging
    std::cout << "Efficiency values:" << std::endl;
    for (auto& entry : efficiencies_by_category) {
        int cat = entry.first;
        std::vector<std::pair<double, double>> mom_efficiencies = entry.second;
        std::cout << "Category " << cat << ":" << std::endl;
        for (const auto& pair : mom_efficiencies) {
            std::cout << "Momentum: " << pair.first << " | Efficiency: " << pair.second << std::endl;
        }
    }

    std::vector<std::vector<double>> cv;
    size_t max_index = 3;
    
    create_covariance_matrix("covariance_data.txt", cv, max_index);
        
    for (size_t i = 0; i <= max_index; ++i) {
      for (size_t j = 0; j <= max_index; ++j) {
	std::cout << cv[i][j] << " ";
      }
      std::cout << std::endl;
    }




    // Create a canvas for plotting
    TCanvas *canvas = new TCanvas("canvas", "Efficiency vs Momentum", 800, 600);

    // Create a legend for the plot
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Plot the efficiency vs momentum for each category
    for (auto& entry : efficiencies_by_category) {
        int cat = entry.first;
        std::vector<std::pair<double, double>> mom_efficiencies = entry.second;
        // Create vectors to hold momenta and efficiencies for plotting
        std::vector<double> momenta_vals;
        std::vector<double> efficiencies_vals;

        for (const auto& pair : mom_efficiencies) {
            momenta_vals.push_back(pair.first);
            efficiencies_vals.push_back(pair.second);
        }

        // Create a graph for the current category
        int npoints = momenta_vals.size();
        TGraph *graph = new TGraph(npoints, &momenta_vals[2], &efficiencies_vals[2]);

        graph->SetTitle(Form("Category %d", cat));
        graph->SetMarkerStyle(20 + cat);  // Different marker for each category
        //graph->SetMarkerColor(cat + 3);   // Different color for each category
        graph->SetLineColor(cat + 3);     // Same color for lines (if connected)

        // Draw the graph on the canvas
        if (cat == 0) {
            graph->Draw("AP");  // First graph, draw with axes
        } else {
            graph->Draw("P");   // Subsequent graphs, only plot points
        }

        // Add the category label to the legend
        legend->AddEntry(graph, Form("Category %d", cat), "p");

        // Clean up memory (graph will be automatically deleted at the end of the function)
    }

    // Add the legend to the canvas
    legend->Draw();

    // Show the canvas
    canvas->Update();
    canvas->Draw();
}    

