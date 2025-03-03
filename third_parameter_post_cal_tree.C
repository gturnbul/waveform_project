// Standard C++ libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <utility> // For std::pair
#include <cmath>         // For math functions (with _USE_MATH_DEFINES)
#include <cassert>       // For debugging assertions
#include <numeric>
#include <set>

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

// ROOT Graphing and Histogram Libraries
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPie.h>
#include <TPieSlice.h>

// ROOT Random Number Generator
#include <TRandom3.h>

// ROOT Fitting Libraries
#include <TF1.h>

// Additional ROOT Utilities
#include <TPaveText.h>
#include <TArrayF.h>

// Boost for hash_map
#include <boost/functional/hash.hpp> // For boost::hash

// Macro Definition
#define _USE_MATH_DEFINES

using namespace std;

void waveform_tree() {
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 0: Open the existing ROOT file with the original tree               ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TFile *file = TFile::Open("waveform_gill_calibration.root", "READ");
    if (!file || file->IsZombie()) {
        cout << "Error: Could not open input file!" << endl;
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 1: Access the original tree                                         ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TTree *tree = (TTree*)file->Get("treebuild");
    if (!tree) {
        cout << "Error: Could not find tree in file!" << endl;
        file->Close();
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 2: Set up variables to read the branches of the original tree       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<int> *om_number = nullptr;
    std::vector<int> *is_electron = nullptr;
    std::vector<double> *feb_energy = nullptr;
    std::vector<int> *feb_charge = nullptr;
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *feb_amplitude = nullptr;
    std::vector<int> *calo_fcr = nullptr;

    // Set the branch addresses
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("is_electron", &is_electron);
    tree->SetBranchAddress("energy", &feb_energy);
    tree->SetBranchAddress("charge", &feb_charge);
    tree->SetBranchAddress("wave", &wave);
    tree->SetBranchAddress("amplitude", &feb_amplitude);
    tree->SetBranchAddress("calo_fcr", &calo_fcr);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 3: Create a new ROOT file to store the new tree and new tree        ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TFile *newFile = new TFile("waveform_tree_gill.root", "RECREATE");
    TTree *newTree = new TTree("treebuild", "treebuild");

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 4: Set variables to hold the new branches of the new tree           ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    int event = -1;
    int electron = 0;
    int om_num = 0;
    double energy = 0;
    int charge = 0;
    int amplitude = 0;
    int fcr_num = 0;

    double baseline = 0;
    double wave_min_value = 0;
    double wave_min_time = 0;
    double pulse_onset = 0;
    double pulse_end = 0;
    double calc_amplitude = 0;
    double chi2 = 0;
    double ndf = 0;
    double chi2ndf = 0;

    double fit_a = 0;
    double fit_p = 0;
    double fit_r = 0;
    double fit_d = 0;
    double fit_s = 0;


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 5: Set the branch addresses                                         ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    newTree->Branch("event", &event);
    newTree->Branch("electron", &electron);
    newTree->Branch("om_number", &om_num);
    newTree->Branch("energy", &energy);
    newTree->Branch("charge", &charge);
    newTree->Branch("amplitude", &amplitude);
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


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 6: Fill the tree with data                                          ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    int nEntries = tree->GetEntries();  // Get the actual number of entries
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Loop over all the elements in the vectors
        for (size_t j = 0; j < is_electron->size(); j++) {        
            
            baseline = 0;
            int FCR_index = (*calo_fcr)[j];  // Get the memory cell start value
            fcr_num = FCR_index;

           ////////////// Calculate BASELINE    //////////////////////////////////////////////////////
            for (int k = 0; k < 96; k++) {
                baseline += wave->at(j).at(k);  // Access wave for event i
            }
            baseline = baseline / 96;  // Calculate average baseline

            double sum_sq_diff = 0;
            for (int k = 0; k < 96; k++) {
                sum_sq_diff += pow(wave->at(j).at(k) - baseline, 2);  // Calculate variance
            }
            double std_dev = sqrt(sum_sq_diff / 95);  // Standard deviation 

            ////////////// Calculate WAVE MINIMUM    /////////////////////////////////////////////////
            wave_min_value = wave->at(j).at(0);  // Initialize with the first value of the waveform
            wave_min_time = 0;  // Initialize with the first time bin (0-based)

            for (int k = 0; k < wave->at(j).size(); k++) {  // Start from 1 because we already initialized with the first value
                double current_value = wave->at(j).at(k);

                // If the current value is smaller than the current minimum, update the minimum
                if (current_value < wave_min_value) {
                    wave_min_value = current_value;
                    wave_min_time = k;  // Store the bin index of the minimum value
                }
            }

            ////////////// Calculate PULSE ONSET    /////////////////////////////////////////////////
            int count_decrease = 0;
            pulse_onset = 0;

            for(int k = 1; k < wave->at(j).size(); k++) {  // Start from 1 because we already initialized with the first value
                double current_value = wave->at(j).at(k);
                double previous_value = wave->at(j).at(k-1);

                // If the current value is smaller than the previous value, increment the counter
                if (current_value < previous_value) {
                    count_decrease++;
                }
                else {
                    count_decrease = 0;
                }
                // check there are 5 decrease in a row
                if (count_decrease == 5) {
                    pulse_onset = k-4;
                    break;
                }
            }

            ////////////// Calculate PULSE END     /////////////////////////////////////////////////  
            int count_stabilise = 0;
            pulse_end = 0;

            for (int k = wave_min_time; k < wave->at(j).size(); k++) {  // Start after wave_min_time
                double current_value = wave->at(j).at(k);

                // Check if the current value is close to the baseline
                if (std::abs(current_value - baseline) < std_dev) {  // Within std dev of baseline
                    count_stabilise++;
                } else {
                    count_stabilise = 0;  // Reset counter if the condition is not met
                }

                // Check if 5 consecutive points fall within the baseline range
                if (count_stabilise == 5) {
                    pulse_end = k - 4;  // First point in the 5-point sequence
                    break;
                }
            }
            
            ////////////// Calculate AMPLITUDE     /////////////////////////////////////////////////
            double calc_amp_time = 0;  // Initialize the time for when threshold is crossed
            calc_amplitude = baseline - wave_min_value;

            ////////////// Calculate AMP Threshold    /////////////////////////////////////////////// 
            double target_amplitude = calc_amplitude * 0.333;  // select a threshold limit of this difference
            double threshold = baseline - target_amplitude;

            for (int k = wave_min_time; k < wave->at(j).size(); k++) {
                double current_value = wave->at(j).at(k);  // Value of the waveform at current bin

                // Check if the waveform rises above the threshold (crosses back up toward baseline)
                if (current_value >= threshold) {
                    calc_amp_time = k;  // Record the time when the threshold is crossed
                    break;  // Exit the loop once the threshold is crossed
                }
            }
            
            /////////////////////////////////////////////////////////////////////////////////////////
            ////////////// Calculate fit with HISTOGRAM /////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////
            TH1D *h_wave = new TH1D(Form("h_wave_%d", j), Form("h_wave_%d", j), wave->at(j).size(), 0, wave->at(j).size());
            for (int k=1; k < wave->at(j).size(); k++) {
                h_wave->SetBinContent(k, wave->at(j).at(k-1));
                h_wave->SetBinError(k, std_dev);
            }

            double start = pulse_onset - 200;
            double end = pulse_end + 200;
            ////////////// Declare Fit Function     /////////////////////////////////////////////////     
            TF1* fitFunc = new TF1("fitFunc", Form("(1/(1-exp(-[4]*(x-[1]))))*[0]*(exp(-(x-[1])/[2]) - exp(-(x-[1])/[3])) + %f", baseline), 0, 1000);
            fitFunc->SetParameters(calc_amplitude*-1.21, wave_min_time-6, 42, 3, 0.6);       

            ////////////// Perform Fit   /////////////////////////////////////////////////////////// 
            h_wave->Fit("fitFunc", "RQ"); //fit the range, quietly

            ////////////// Find Fit Parameters   /////////////////////////////////////////////////// 
            fit_a = fitFunc->GetParameter(0);
            fit_p = fitFunc->GetParameter(1);
            fit_r = fitFunc->GetParameter(2);
            fit_d = fitFunc->GetParameter(3);
            fit_s = fitFunc->GetParameter(4);

            chi2 = fitFunc->GetChisquare();
            ndf = fitFunc->GetNDF();
            if (ndf > 0) {
                chi2ndf = chi2 / ndf;
            } else {
                chi2ndf = 0; // Set to 0 if no valid degrees of freedom
            }

            // /////////////////////////////////////////////////////////////////////////////////////////
            // ////////////// CHI2NDF VARIATIONS           /////////////////////////////////////////////
            // /////////////////////////////////////////////////////////////////////////////////////////                     
        

            /////////////////////////////////////////////////////////////////////////////////////////
            ////////////// ORIGINGAL VALUES             /////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////  
            om_num = (*om_number)[j];
            energy = (*feb_energy)[j];
            charge = (*feb_charge)[j];
            amplitude = (*feb_amplitude)[j];
            fcr_num = (*calo_fcr)[j];
            event = i;  // Optionally track the event number
            electron = (*is_electron)[j];

            // Fill the new tree with data
            newTree->Fill();
            delete h_wave;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 7: Write the tree to the file                                       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    newTree->Write();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 8: Close the file to save everything                                ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    newFile->Close();
    file->Close();
}

// Necessary for the makefile      /////////////////////////////////////////////////////////////////
int main() {
    waveform_tree();
    return 0;
}
