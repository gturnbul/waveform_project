void treechi2(){

    // Open the treebuilding root file
    TFile *file = TFile::Open("treebuild.root", "READ");
    if (!file || file->IsZombie()) {
        cout << "Error: Could not open input file!" << endl;
        return;
    }

    // Access the treebuild tree
    TTree *tree = (TTree*)file->Get("treebuild");
    if (!tree) {
        cout << "Error: Could not find tree in file!" << endl;
        file->Close();
        return;
    }

    //create canva to draw histograms
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    // Create histograms for each branch
    TH1F *hist1 = new TH1F("hist1", "Chi2/ndf of Electron", 100, 0, 100);
    TH1F *hist2 = new TH1F("hist2", "Chi2 Comparison of Fixed and Set A parameter", 100, 0, 100);

    // Set histogram colors
    hist1->SetLineColor(kBlue);   // Red for e_chi2ndf
    hist1->SetStats(0);  // Turn off the statistics box
    hist1->GetXaxis()->SetTitle("Chi2/ndf");
    hist1->GetYaxis()->SetTitle("Counts");
    hist2->SetLineColor(kRed);  // Blue for e_chi2ndf_set
    hist2->SetStats(0);  // Turn off the statistics box
    hist2->GetXaxis()->SetTitle("Chi2/ndf");
    hist2->GetYaxis()->SetTitle("Counts");

    // Fill histograms with the respective branches
    tree->Draw("e_chi2ndf>>hist1", "e_event>-0.1");
    // tree->Draw("e_chi2ndf_set>>hist2", "e_event>-0.1", "same");

    // Draw histograms on the same canvas
    hist1->Draw();
    // hist2->Draw("same");

    // // Add a legend to differentiate
    TLegend *legend = new TLegend(0.7,0.7,0.9,0.9); // Adjust position as needed
    legend->AddEntry(hist1, "Start Parameter Set", "l");
    legend->AddEntry(hist2, "Amplitude Fixed", "l");
    legend->Draw("same");

    // update the canvas
    c1->Update();
    // Save the canvas as a .png file
    c1->SaveAs("treecollect_electron.png");
    

    // Close the file
    file->Close();
}
