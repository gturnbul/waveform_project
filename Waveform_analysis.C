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

    // Open input ROOT file
    TFile *infile = TFile::Open("filtered_gamma_fcr.root", "READ");
    if (!infile || infile->IsZombie()) { std::cerr << "Cannot open input file!\n"; return 1; }

    TTree *tree = (TTree*) infile->Get("Result_tree");
    if (!tree) { std::cerr << "Cannot find Result_tree!\n"; return 1; }

    // Branches
    std::vector<int> *om_number = nullptr;
    std::vector<int> *is_electron = nullptr;
    std::vector<double> *energy = nullptr;
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *charge = nullptr;
    std::vector<int> *amplitude = nullptr;
    std::vector<unsigned short> *calo_fcr = nullptr;

    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("is_electron", &is_electron);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("wave", &wave);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("amplitude", &amplitude);
    tree->SetBranchAddress("calo_fcr", &calo_fcr);


    // Step 0: Create output file and tree
    TFile *newFile = new TFile("waveform_analysis_1112.root", "RECREATE");
    TTree *newTree = new TTree("wave", "wave");

    // Define branches
    int event = -1, electron = 0, om_num = 0, fcr_num = 0;
    double baseline = 0, wave_min_value = 0, wave_min_time = 0;
    double pulse_onset = 0, pulse_end = 0, calc_amplitude = 0;
    double chi2 = 0, ndf = 0, chi2ndf = 0;
    double fit_a = 0, fit_p = 0, fit_r = 0, fit_d = 0, fit_s = 0;
    double out_energy = 0, out_charge = 0, out_amplitude = 0;

    newTree->Branch("event", &event);
    newTree->Branch("electron", &electron);
    newTree->Branch("om_number", &om_num);
    newTree->Branch("fcr_num", &fcr_num);
    newTree->Branch("baseline", &baseline);
    newTree->Branch("wave_min_value", &wave_min_value);
    newTree->Branch("wave_min_time", &wave_min_time);
    newTree->Branch("pulse_onset", &pulse_onset);
    newTree->Branch("pulse_end", &pulse_end);
    newTree->Branch("calc_amplitude", &calc_amplitude);
    newTree->Branch("chi2", &chi2);
    newTree->Branch("ndf", &ndf);
    newTree->Branch("chi2ndf", &chi2ndf);
    newTree->Branch("fit_a", &fit_a);
    newTree->Branch("fit_p", &fit_p);
    newTree->Branch("fit_r", &fit_r);
    newTree->Branch("fit_d", &fit_d);
    newTree->Branch("fit_s", &fit_s);
    newTree->Branch("energy", &out_energy);
    newTree->Branch("charge", &out_charge);
    newTree->Branch("amplitude", &out_amplitude);


    int image_count = 0;

    // Loop over entries
    int nentries = tree->GetEntries();
     std::cout << "Total entries in tree: " << nentries << std::endl;

     for (int i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        size_t nhits = om_number->size();
        for (size_t j = 0; j < nhits; ++j) {
            int om_num = om_number->at(j);
            std::vector<short> original_waveform = wave->at(j);

            // --- EOW offsets ---
            std::vector<short> eow_corrected = original_waveform;
            if (eow_offsets.count(om_num)) {
                for (int bin = 976; bin <= 1023; ++bin) {
                    if (eow_corrected[bin] != EMPTY_BIN)
                        eow_corrected[bin] -= eow_offsets[om_num][bin - 976];
                }
            }

            // --- Reorder ---
            int fcr_index = static_cast<int>(calo_fcr->at(j));
            std::vector<short> reordered = reorder_waveform(eow_corrected, fcr_index);

            // --- FEB offsets ---
            if (feb_offsets.count(om_num)) {
                for (int bin = 0; bin < 1024; ++bin) {
                    if (reordered[bin] != EMPTY_BIN)
                        reordered[bin] -= feb_offsets[om_num][bin];
                }
            }

            // --- Inverse reorder ---
            std::vector<short> corrected = inverse_reorder_waveform(reordered, fcr_index);
            wave->at(j) = corrected;

            // --- Plot first 10 waveforms ---
            if (image_count < 10) {
                TCanvas *c = new TCanvas("c", "Waveform Comparison", 800, 800);
                c->Divide(1, 2);

                std::vector<double> xbins(original_waveform.size());
                for (size_t b = 0; b < xbins.size(); ++b) xbins[b] = b;

                TGraph *g_before = new TGraph(xbins.size(), xbins.data(), std::vector<double>(original_waveform.begin(), original_waveform.end()).data());
                TGraph *g_after = new TGraph(xbins.size(), xbins.data(), std::vector<double>(corrected.begin(), corrected.end()).data());

                c->cd(1); g_before->SetTitle("Waveform BEFORE Corrections;Bin;ADC"); g_before->SetLineColor(kRed); g_before->Draw("AL");
                c->cd(2); g_after->SetTitle("Waveform AFTER Corrections;Bin;ADC"); g_after->SetLineColor(kBlue); g_after->Draw("AL");

                c->SaveAs(Form("waveform_split_%d.png", image_count));
                delete c;
                image_count++;
            }
        }

        //fill the tree
        newTree->Fill();

    }
    newFile->cd();
    newTree->Write();
    newFile->Close();
    infile->Close();

    return 0;
}
