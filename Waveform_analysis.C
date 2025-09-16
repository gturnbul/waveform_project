#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <memory>
#include <string>
#include <cmath>         // For math functions (with _USE_MATH_DEFINES)
#include <sys/stat.h>    // For file system operations
#include <sys/types.h>   // For file system types
#include <cassert>       // For debugging assertions
#include <thread>        // For sleeping

// ROOT Core Libraries
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TKey.h"
#include "TParameter.h"
#include <TString.h>      // For ROOT string handling
#include <TObjArray.h>    // For working with ROOT arrays

// ROOT Graphical Libraries
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <TPaletteAxis.h>
#include <TText.h>
#include <TLatex.h>
#include <TLegend.h>
#include "TMultiGraph.h"

// ROOT Graphing and Histogram Libraries
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraphErrors.h>

// ROOT Random Number Generator
#include <TRandom3.h>

// ROOT Fitting Libraries
#include <TF1.h>

// Additional ROOT Utilities
#include <TPaveText.h>
#include <TArrayF.h>
#include <unordered_map>
#include <TVirtualFFT.h>
#include <complex>
#include <TComplex.h>
#include <map>
#include <vector>
#include <cmath>
#include <THStack.h>

#include "TSpectrum.h"

#include <cstdlib> // for rand()
// Macro Definition
#define _USE_MATH_DEFINES

using namespace std;

const short EMPTY_BIN = -1;
//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Functions //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Function to load FEB offsets from CSV
std::map<int, std::vector<double>> load_feb_offsets(const std::string& filename) {
    std::map<int, std::vector<double>> feb_offsets;
    std::ifstream infile(filename);
    std::string line;

    // skip header
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        int om, bin;
        double data_mean, fit_mean, sigma, chi2ndf, baseline;
        char comma;

        ss >> om >> comma >> bin >> comma >> data_mean >> comma
           >> fit_mean >> comma >> sigma >> comma >> chi2ndf >> comma >> baseline;

        feb_offsets[om].resize(1024, 0.0);  // ensure 1024 bins per OM
        feb_offsets[om][bin] = fit_mean;    // use fitted mean as offset
    }
    return feb_offsets;
}

// Function to load EOW offsets from CSV
std::map<int, std::vector<double>> load_eow_offsets(const std::string& filename) {
    std::map<int, std::vector<double>> eow_offsets;
    std::ifstream infile(filename);
    std::string line;

    // skip header
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        int om, bin;
        double offset_value, sigma, chi2ndf, baseline;
        char comma;

        ss >> om >> comma >> bin >> comma >> offset_value >> comma
           >> sigma >> comma >> chi2ndf >> comma >> baseline;

        if (bin >= 976 && bin <= 1023) {
            eow_offsets[om].resize(48, 0.0);    // only bins 976–1023
            eow_offsets[om][bin - 976] = offset_value;
        }
    }
    return eow_offsets;
}

// Forward declare reorder functions
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024;
}

std::vector<short> reorder_waveform(const std::vector<short>& waveform, int fcr) {
    std::vector<short> reordered(1024, EMPTY_BIN);  // initialize all bins as empty

    int max_bins = std::min((int)waveform.size(), 976);  // only use bins 0 to 975

    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);
        reordered[reordered_index] = waveform[i];
    }

    // bins 976 to 1023 remain EMPTY_BIN, indicating gaps

    return reordered;
}

std::vector<short> inverse_reorder_waveform(const std::vector<short>& reordered, int fcr) {
    int max_bins = 1024;  // original waveform length
    
    std::vector<short> original(max_bins, EMPTY_BIN);
    
    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);
        // Safety check: reordered_index should be within bounds
        if (reordered_index >= 0 && reordered_index < (int)reordered.size()) {
            original[i] = reordered[reordered_index];
        } else {
            // handle unexpected index (optional)
            original[i] = EMPTY_BIN;
        }
    }
    
    return original;
}
//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Main ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
int main() {
    // Load CSV offsets
    auto feb_offsets = load_feb_offsets("feb_offsets_1143_width.csv");
    auto eow_offsets = load_eow_offsets("eow_spike_offsets_1143_width.csv");

    // Open input ROOT file (run 1112)
    TFile *infile = new TFile("snemo_run-1112_udd.root", "READ");
    TTree *tree = (TTree*) infile->Get("SimData");

    // Branches
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *calo_wall = nullptr;
    std::vector<int> *calo_side = nullptr;
    std::vector<int> *calo_column = nullptr;
    std::vector<int> *calo_row = nullptr;
    std::vector<int> *calo_type = nullptr;
    std::vector<int> *fcr = nullptr;

    int calo_nohits = 0;

    tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
    tree->SetBranchAddress("digicalo.waveform", &wave);
    tree->SetBranchAddress("digicalo.wall", &calo_wall);
    tree->SetBranchAddress("digicalo.side", &calo_side);
    tree->SetBranchAddress("digicalo.column", &calo_column);
    tree->SetBranchAddress("digicalo.row", &calo_row);
    tree->SetBranchAddress("digicalo.type", &calo_type);
    tree->SetBranchAddress("digicalo.fcr", &fcr);

    // Output file
    TFile *outfile = new TFile("snemo_run-1112_corrected.root", "RECREATE");
    TTree *outtree = tree->CloneTree(0); // clone structure

    int image_count = 0; // put this before the event loop

    // Loop over entries
    int nentries = tree->GetEntries();
     std::cout << "Total entries in tree: " << nentries << std::endl;

    for (int i = 0; i < nentries; i++) {
        tree->GetEntry(i);

       for (int j = 0; j < calo_nohits; ++j) {
            int om_num = calo_side->at(j) * 260 + calo_column->at(j) * 13 + calo_row->at(j);

            // Step 0: read waveform and save a true copy from the tree
            std::vector<short> original_waveform = wave->at(j);  // original waveform from TTree

            // Step 1: create a copy for EOW correction
            std::vector<short> eow_corrected = original_waveform;

            // Step 2: apply EOW offsets to last 48 bins
            if (eow_offsets.count(om_num)) {
                for (int bin = 976; bin <= 1023; ++bin) {
                    if (eow_corrected[bin] != EMPTY_BIN) {
                        eow_corrected[bin] -= eow_offsets[om_num][bin - 976];
                    }
                }
            }

            // Step 2: reorder spike-corrected waveform
            std::vector<short> reordered = reorder_waveform(eow_corrected, fcr->at(j));

            // Step 3: apply FEB offsets
            if (feb_offsets.count(om_num)) {
                for (int bin = 0; bin < 1024; ++bin) {
                    if (reordered[bin] != EMPTY_BIN) {
                        reordered[bin] -= feb_offsets[om_num][bin];
                    }
                }
            }

            // Step 4: inverse reorder back
            std::vector<short> corrected = inverse_reorder_waveform(reordered, fcr->at(j));

            // Replace waveform in tree
            wave->at(j) = corrected;

            // --- optional plotting (only first 10) ---
            if (image_count < 10) {
                TCanvas *c = new TCanvas("c", "Waveform Comparison", 800, 800);
                c->Divide(1, 2); // top=before, bottom=after

                // X values
                std::vector<double> xbins(original_waveform.size());
                for (size_t b = 0; b < xbins.size(); b++) xbins[b] = b;

                // Graph BEFORE
                std::vector<double> y_before(original_waveform.begin(), original_waveform.end());
                TGraph *g_before = new TGraph(xbins.size(), xbins.data(), y_before.data());

                // Graph AFTER
                std::vector<double> y_after(corrected.begin(), corrected.end());
                TGraph *g_after = new TGraph(xbins.size(), xbins.data(), y_after.data());

                // Top pad: before
                c->cd(1);
                g_before->SetTitle("Waveform BEFORE Corrections;Bin;ADC");
                g_before->SetLineColor(kRed);
                g_before->Draw("AL");

                // Bottom pad: after
                c->cd(2);
                g_after->SetTitle("Waveform AFTER Corrections;Bin;ADC");
                g_after->SetLineColor(kBlue);
                g_after->Draw("AL");

                // Save
                c->SaveAs(Form("waveform_split_%d.png", image_count));

                delete c;
                image_count++;
            }
            // // --- inside your loop ---
            // if (image_count < 10) {
            //     TCanvas *c = new TCanvas("c", "EOW Correction Check", 800, 800);
            //     c->Divide(1, 2); // top=original, bottom=EOW corrected

            //     std::vector<double> xbins(1024);
            //     for (int b = 0; b < 1024; ++b) xbins[b] = b;

            //     // Top: true original
            //     std::vector<double> y_orig(1024);
            //     for (int b = 0; b < 1024; ++b) y_orig[b] = original_waveform[b];
            //     TGraph *g_orig = new TGraph(1024, xbins.data(), y_orig.data());
            //     c->cd(1);
            //     g_orig->SetTitle("Waveform BEFORE EOW Correction;Bin;ADC");
            //     g_orig->SetLineColor(kRed);
            //     g_orig->Draw("AL");

            //     // Bottom: after EOW correction
            //     std::vector<double> y_eow(1024);
            //     for (int b = 0; b < 1024; ++b) y_eow[b] = eow_corrected[b];
            //     TGraph *g_eow = new TGraph(1024, xbins.data(), y_eow.data());
            //     c->cd(2);
            //     g_eow->SetTitle("Waveform AFTER EOW Correction;Bin;ADC");
            //     g_eow->SetLineColor(kBlue);
            //     g_eow->Draw("AL");

            //     c->SaveAs(Form("waveform_EOWcheck_%d.png", image_count));
            //     delete c;
            //     ++image_count;
            // }

        }

        // Write this entry to the output tree
        outtree->Fill();

        // Stop outer loop once we've produced 10 images (optional: keeps runtime small)
        if (image_count >= 10) break;
    }

    // Finish up
    outfile->cd();
    outtree->Write();
    outfile->Close();
    infile->Close();

    return 0;
}
