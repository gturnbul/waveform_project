#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
#include <TPaletteAxis.h>
#include <TArrayF.h>

using namespace std;

void treebuilding(){

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 0: Open the existing ROOT file with the original tree               ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    TFile *file = TFile::Open("filtered_gamma_fcr.root", "READ");  // input file path
    if (!file || file->IsZombie()) {
        cout << "Error: Could not open input file!" << endl;
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 1: Access the original tree                                         ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
   
    TTree *tree = (TTree*)file->Get("Result_tree");
    if (!tree) {
        cout << "Error: Could not find tree in file!" << endl;
        file->Close();
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 2: Set up variables to read the branches of the original tree       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

     // Define variables to hold the waveform data, electron type, and optical module number
    std::vector<int> *om_number = nullptr;
    std::vector<int> *is_electron = nullptr;
    std::vector<double> *energy = nullptr;
    std::vector<int> *charge = nullptr;
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *amplitude = nullptr;
    std::vector<int> *calo_fcr = nullptr;

    // Set the branch addresses
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("is_electron", &is_electron);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("wave", &wave);
    tree->SetBranchAddress("amplitude", &amplitude);
    tree->SetBranchAddress("calo_fcr", &calo_fcr);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 3: Create a new ROOT file to store the new tree and new tree        ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Create a new ROOT file to store the new tree
    TFile *newFile = new TFile("treebuild_fcr.root", "RECREATE");    
    // create new tree
    TTree *newTree = new TTree("treebuild", "treebuild");

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 4: Set variables to hold the new branches of the new tree           ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    int e_event = -1;
    int e_om_number = 0;
    double e_energy = 0;
    int e_charge = 0;
    int e_amplitude = 0;
    int e_calo_fcr = 0;
    double e_baseline = 0;
    double e_wave_peak_value = 0;
    int e_wave_peak_time = 0;
    double e_calc_amplitude = 0;
    double e_pulse_onset = 0;
    double e_chi2 = 0;
    double e_ndf = 0;
    double e_chi2ndf = 0;
    double e_fit_a_factor = 0;
    double e_fit_pulse_onset = 0;
    double e_fit_r_factor = 0;
    double e_fit_d_factor = 0;
    double e_fit_s_factor = 0;

    int g_event = -1;
    int g_om_number = 0;
    double g_energy = 0;
    int g_charge = 0;
    int g_amplitude = 0;
    int g_calo_fcr = 0;
    double g_baseline = 0;
    double g_wave_peak_value = 0;
    int g_wave_peak_time = 0;
    double g_calc_amplitude = 0;
    double g_pulse_onset = 0;
    double g_chi2 = 0;
    double g_ndf = 0;
    double g_chi2ndf = 0;
    double g_fit_a_factor = 0;
    double g_fit_pulse_onset = 0;
    double g_fit_r_factor = 0;
    double g_fit_d_factor = 0;
    double g_fit_s_factor = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 5: Set the branch addresses                                         ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    newTree->Branch("e_event", &e_event);
    newTree->Branch("e_om_number", &e_om_number);
    newTree->Branch("e_energy", &e_energy);
    newTree->Branch("e_charge", &e_charge);
    newTree->Branch("e_amplitude", &e_amplitude);
    newTree->Branch("e_calo_fcr", &e_calo_fcr);
    newTree->Branch("e_baseline", &e_baseline);
    newTree->Branch("e_wave_peak_value", &e_wave_peak_value);
    newTree->Branch("e_wave_peak_time", &e_wave_peak_time);
    newTree->Branch("e_calc_amplitude", &e_calc_amplitude);
    newTree->Branch("e_pulse_onset", &e_pulse_onset);
    newTree->Branch("e_chi2", &e_chi2);
    newTree->Branch("e_ndf", &e_ndf);
    newTree->Branch("e_chi2ndf", &e_chi2ndf);
    newTree->Branch("e_fit_a_factor", &e_fit_a_factor);
    newTree->Branch("e_fit_pulse_onset", &e_fit_pulse_onset);
    newTree->Branch("e_fit_r_factor", &e_fit_r_factor);
    newTree->Branch("e_fit_d_factor", &e_fit_d_factor);
    newTree->Branch("e_fit_s_factor", &e_fit_s_factor);

    newTree->Branch("g_event", &g_event);
    newTree->Branch("g_om_number", &g_om_number);
    newTree->Branch("g_energy", &g_energy);
    newTree->Branch("g_charge", &g_charge);
    newTree->Branch("g_amplitude", &g_amplitude);
    newTree->Branch("g_calo_fcr", &g_calo_fcr);
    newTree->Branch("g_baseline", &g_baseline);
    newTree->Branch("g_wave_peak_value", &g_wave_peak_value);
    newTree->Branch("g_wave_peak_time", &g_wave_peak_time);
    newTree->Branch("g_calc_amplitude", &g_calc_amplitude);
    newTree->Branch("g_pulse_onset", &g_pulse_onset);
    newTree->Branch("g_chi2", &g_chi2);
    newTree->Branch("g_ndf", &g_ndf);
    newTree->Branch("g_chi2ndf", &g_chi2ndf);
    newTree->Branch("g_fit_a_factor", &g_fit_a_factor);
    newTree->Branch("g_fit_pulse_onset", &g_fit_pulse_onset);
    newTree->Branch("g_fit_r_factor", &g_fit_r_factor);
    newTree->Branch("g_fit_d_factor", &g_fit_d_factor);
    newTree->Branch("g_fit_s_factor", &g_fit_s_factor);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 6: Fill the tree with data                                          ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Determine the number of entries in the tree
    Long64_t nEntries = tree->GetEntries();

    // Set waveform length
    int ticks = 1024;
    

    for (Long64_t j = 0; j < nEntries; j++) {
        tree->GetEntry(j);

        double e_std_dev = 0;
        e_baseline =0;
        double g_std_dev = 0;
        g_baseline = 0;
        // Loop over elements in the `electron` vector for the current event
        for (size_t i = 0; i < is_electron->size(); i++) {
            if (is_electron->at(i) == 1) { // Electron histogram

               int FCR_start_index = calo_fcr->at(i); // Get the memory cell start value

                // Save the FCR_start_index to the new tree
                e_calo_fcr= FCR_start_index;

                /// i need to know what the offsets per cell of each om are so that I can reduce the offset error on each waveform???

                // Calculate the baseline and standard deviation of the waveform
                for (int k = 0; k < 96; k++) {
                    e_baseline += wave->at(i).at(k);
                }
                e_baseline = e_baseline / 96;
               

                double e_sum_sq_diff = 0;
                for (int k = 0; k < 96; k++) {
                    e_sum_sq_diff += pow(wave->at(i).at(k) - e_baseline, 2);
                }

                double e_std_dev = sqrt(e_sum_sq_diff / 95);
              

                // Find the minimum value of the waveform
                e_wave_peak_value = wave->at(i).at(0);  // Initialise with the first value of the waveform
                e_wave_peak_time = 0;  // Initialize with the first time bin (0-based)

                for (int k = 1; k < wave->at(i).size(); k++) {  // Start from 1 because we already initialized with the first value
                    double current_value = wave->at(i).at(k);

                    // If the current value is smaller than the current minimum, update the minimum
                    if (current_value < e_wave_peak_value) {
                        e_wave_peak_value = current_value;
                        e_wave_peak_time = k;  // Store the bin index of the minimum value
                    }
                }
                
                // Find the pulse onset of the waveform
                e_pulse_onset = wave->at(i).at(0);  // Initialise with the first value of the waveform
                int e_count_decrease = 0;

                for (int k = 1; k < wave->at(i).size(); k++) {  // Start from 1 because we already initialized with the first value
                    double current_value = wave->at(i).at(k);
                    double previous_value = wave->at(i).at(k-1);

                    // If the current value is smaller than the previous value, increment the counter
                    if (current_value < previous_value) {
                        e_count_decrease++;
                    }
                    else {
                        e_count_decrease = 0;
                    }
                    // check there are 5 decrease in a row
                    if (e_count_decrease == 5) {
                        e_pulse_onset = k-4;
                        break;
                    }
                }              

                // Calculate the amplitude of the waveform
                e_calc_amplitude = e_wave_peak_value - e_baseline;  

            
                // Declare the 1D histogram to hold the waveform data
                TH1D *e_h_wave = new TH1D(Form("e_h_wave_%d", j), Form("e_h_wave_%d", j), wave->at(i).size(), 0, wave->at(i).size());
                for (int k=1; k < wave->at(i).size(); k++) {
                    e_h_wave->SetBinContent(k, wave->at(i).at(k-1));
                    e_h_wave->SetBinError(k, e_std_dev);
                }

                // Declare the fit function for the waveform, embedding e_baseline as a constant value
                TF1* e_fitFunc = new TF1("e_fitFunc", Form("(1/(1-exp(-[4]*(x-[1]))))*[0]*(exp(-(x-[1])/[2]) - exp(-(x-[1])/[3])) + %f", e_baseline), 200, 900);
                e_fitFunc->SetParameters(e_calc_amplitude, e_wave_peak_time, 18, 10, 5);

                // Perform fit
                e_h_wave->Fit("e_fitFunc", "RQ"); //fit the range, quietly

                // Get the fit parameters
                e_chi2 = e_fitFunc->GetChisquare();
                e_ndf = e_fitFunc->GetNDF();
                e_chi2ndf = e_chi2 / e_ndf;
                e_fit_a_factor = e_fitFunc->GetParameter(0);
                e_fit_pulse_onset = e_fitFunc->GetParameter(1);
                e_fit_r_factor = e_fitFunc->GetParameter(2);
                e_fit_d_factor = e_fitFunc->GetParameter(3);
                e_fit_s_factor = e_fitFunc->GetParameter(4);
                

                // Get the other parameters of the waveform from the original file
                e_event = j;
                e_om_number = om_number->at(i);
                e_energy = energy->at(i);
                e_charge = charge->at(i);
                e_amplitude = amplitude->at(i);

                // Fill the tree with the data
                newTree->Fill();

                // Clean up
                delete e_h_wave;  // Prevent memory leaks by deleting the histogram after saving
            

                // Reset the variables for the next loop iteration
                e_event = -1;
            


            } 
            else if (is_electron->at(i) == 0) { // Gamma histogram

                int FCR_start_index = calo_fcr->at(i); // Get the memory cell start value

                // Save the FCR_start_index to the new tree
                g_calo_fcr= FCR_start_index;

                // Calculate the baseline and standard deviation of the waveform
                for (int k = 0; k < 96; k++) {
                    g_baseline += wave->at(i).at(k);
                }
                g_baseline = g_baseline / 96;
               

                double g_sum_sq_diff = 0;
                for (int k = 0; k < 96; k++) {
                    g_sum_sq_diff += pow(wave->at(i).at(k) - g_baseline, 2);
                }

                double g_std_dev = sqrt(g_sum_sq_diff / 95);
           

                // Find the minimum value of the waveform
                g_wave_peak_value = wave->at(i).at(0);  // Initialize with the first value of the waveform
                g_wave_peak_time = 0;  // Initialize with the first time bin (0-based)

                for (int k = 1; k < wave->at(i).size(); k++) {  // Start from 1 because we already initialized with the first value
                    double current_value = wave->at(i).at(k);

                    // If the current value is smaller than the current minimum, update the minimum
                    if (current_value < g_wave_peak_value) {
                        g_wave_peak_value = current_value;
                        g_wave_peak_time = k;  // Store the bin index of the minimum value
                    }
                } 
                
                // Find the pulse onset of the waveform
                g_pulse_onset = wave->at(i).at(0);  // Initialise with the first value of the waveform
                int g_count_decrease = 0;

                for (int k = 1; k < wave->at(i).size(); k++) {  // Start from 1 because we already initialized with the first value
                    double current_value = wave->at(i).at(k);
                    double previous_value = wave->at(i).at(k-1);

                    // If the current value is smaller than the previous value, increment the counter
                    if (current_value < previous_value) {
                        g_count_decrease++;
                    }
                    else {
                        g_count_decrease = 0;
                    }
                    // check there are 5 decrease in a row
                    if (g_count_decrease == 5) {
                        g_pulse_onset = k-4;
                        break;
                    }
                }   
             
                // Calculate the amplitude of the waveform
                g_calc_amplitude = g_wave_peak_value - g_baseline;  

                // Declare the 1D histogram to hold the waveform data
                TH1D *g_h_wave = new TH1D(Form("g_h_wave_%d", j), Form("g_h_wave_%d", j), wave->at(i).size(), 0, wave->at(i).size());
                for (int k=1; k < wave->at(i).size(); k++) {
                    g_h_wave->SetBinContent(k, wave->at(i).at(k-1));
                    g_h_wave->SetBinError(k, g_std_dev);
                }

                // Declare the fit function for the waveform, embedding g_baseline as a constant value
                TF1* g_fitFunc = new TF1("g_fitFunc", Form("(1/(1-exp(-[4]*(x-[1]))))*[0]*(exp(-(x-[1])/[2]) - exp(-(x-[1])/[3])) + %f", g_baseline), 200, 900);
                g_fitFunc->SetParameters(g_calc_amplitude, g_wave_peak_time, 18, 10, 5);

                // Perform fit
                g_h_wave->Fit("g_fitFunc", "RQ"); //fit the range, quietly

                // Get the fit parameters
                g_chi2 = g_fitFunc->GetChisquare();
                g_ndf = g_fitFunc->GetNDF();
                g_chi2ndf = g_chi2 / g_ndf;
                g_fit_a_factor = g_fitFunc->GetParameter(0);
                g_fit_pulse_onset = g_fitFunc->GetParameter(1);
                g_fit_r_factor = g_fitFunc->GetParameter(2);
                g_fit_d_factor = g_fitFunc->GetParameter(3);
                g_fit_s_factor = g_fitFunc->GetParameter(4);
                
                // Get the other parameters of the waveform from the original file

                g_event = j;
                g_om_number = om_number->at(i);
                g_energy = energy->at(i);
                g_charge = charge->at(i);
                g_amplitude = amplitude->at(i);

                // Fill the tree with the data
                newTree->Fill();

                // Clean up
                delete g_h_wave;  // Prevent memory leaks by deleting the histogram after saving
                
                // Reset the variables for the next loop iteration
                g_event = -1;
              

            }
          
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 7: Write the tree to the file                                       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    newTree->Write();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 8: Close the file to save everything                                ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    file->Close();
}

// Necessary for the makefile                                                   ////////////////////
int main() {
    treebuilding();
    return 0;
}
    treebuilding();
    return 0;
}
