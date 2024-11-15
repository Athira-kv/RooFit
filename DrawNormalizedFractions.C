void DrawNormalizedFractions(TH2I* hist) {
    // Get total number of entries
    double totalEntries = hist->GetEntries();
    
    if (totalEntries == 0) {
        std::cout << "No entries in the histogram." << std::endl;
        return;
    }
    
    // Create a new TH2F with the same bins to store the normalized values
    TH2F* histFrac = new TH2F("histFrac", "Histogram with Fractions", 
                              hist->GetNbinsX(), 0, 4,
                              hist->GetNbinsY(), 0, 4);
    
    // Loop over all bins to normalize bin content by total entries
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= hist->GetNbinsY(); ++j) {
            double binContent = hist->GetBinContent(i, j);
            double fraction = binContent / totalEntries;
            histFrac->SetBinContent(i, j, fraction);  // Set as fraction in TH2F
        }
    }

    // Set text display format to show fractions with two decimal places
    gStyle->SetPaintTextFormat("0.4f");
    histFrac->SetMarkerSize(1.5);  // Adjust text size if needed

    // Draw the new histogram with the TEXT option
    histFrac->Draw("TEXT COLZ");
}
