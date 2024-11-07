#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

using namespace std;

void treebuilding(){

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 0: Open the existing ROOT file with the original tree               ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    TFile *file = TFile::Open("filtered_gamma_e.root", "READ");  // input file path
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
    std::vector<int> *electron = nullptr;
    std::vector<double> *energy = nullptr;
    std::vector<int> *charge = nullptr;
    std::vector<std::vector<short>> *wave = nullptr;

    // Set the branch addresses
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("electron", &electron);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("wave", &wave);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 3: Create a new ROOT file to store the new tree and new tree        ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Create a new ROOT file to store the new tree
    TFile *newFile = new TFile("treebuild.root", "RECREATE");    
    // create new tree
    TTree *newTree = new TTree("treebuild", "treebuild");

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 4: Set variables to hold the new branches of the new tree           ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    int e_event = -1;
    int e_om_number = 0;
    double e_energy = 0;
    int e_charge = 0;
    double e_baseline = 0;
    double e_wave_peak_value = 0;
    int e_wave_peak_time = 0;
    double e_calc_amplitude = 0;
    double e_chi2 = 0;
    double e_ndf = 0;
    double e_chi2ndf = 0;
    double e_fit_a_factor = 0;
    double e_fit_pulse_onset = 0;
    double e_fit_r_factor = 0;
    double e_fit_d_factor = 0;

    int g_event = -1;
    int g_om_number = 0;
    double g_energy = 0;
    int g_charge = 0;
    double g_baseline = 0;
    double g_wave_peak_value = 0;
    int g_wave_peak_time = 0;
    double g_calc_amplitude = 0;
    double g_chi2 = 0;
    double g_ndf = 0;
    double g_chi2ndf = 0;
    double g_fit_a_factor = 0;
    double g_fit_pulse_onset = 0;
    double g_fit_r_factor = 0;
    double g_fit_d_factor = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 5: Set the branch addresses                                         ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    newTree->Branch("e_event", &e_event);
    newTree->Branch("e_om_number", &e_om_number);
    newTree->Branch("e_energy", &e_energy);
    newTree->Branch("e_charge", &e_charge);
    newTree->Branch("e_baseline", &e_baseline);
    newTree->Branch("e_wave_peak_value", &e_wave_peak_value);
    newTree->Branch("e_wave_peak_time", &e_wave_peak_time);
    newTree->Branch("e_calc_amplitude", &e_calc_amplitude);
    newTree->Branch("e_chi2", &e_chi2);
    newTree->Branch("e_ndf", &e_ndf);
    newTree->Branch("e_chi2ndf", &e_chi2ndf);
    newTree->Branch("e_fit_a_factor", &e_fit_a_factor);
    newTree->Branch("e_fit_pulse_onset", &e_fit_pulse_onset);
    newTree->Branch("e_fit_r_factor", &e_fit_r_factor);
    newTree->Branch("e_fit_d_factor", &e_fit_d_factor);

    newTree->Branch("g_event", &g_event);
    newTree->Branch("g_om_number", &g_om_number);
    newTree->Branch("g_energy", &g_energy);
    newTree->Branch("g_charge", &g_charge);
    newTree->Branch("g_baseline", &g_baseline);
    newTree->Branch("g_wave_peak_value", &g_wave_peak_value);
    newTree->Branch("g_wave_peak_time", &g_wave_peak_time);
    newTree->Branch("g_calc_amplitude", &g_calc_amplitude);
    newTree->Branch("g_chi2", &g_chi2);
    newTree->Branch("g_ndf", &g_ndf);
    newTree->Branch("g_chi2ndf", &g_chi2ndf);
    newTree->Branch("g_fit_a_factor", &g_fit_a_factor);
    newTree->Branch("g_fit_pulse_onset", &g_fit_pulse_onset);
    newTree->Branch("g_fit_r_factor", &g_fit_r_factor);
    newTree->Branch("g_fit_d_factor", &g_fit_d_factor);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 6: Fill the tree with data                                          ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Determine the number of entries in the tree
    Long64_t nEntries = 100; // tree->GetEntries();

    // Set waveform length
    int ticks = 1024;
    
    // Define histograms for electrons and gammas with Event number on the x-axis and Waveform index on the y-axis
    TH2F *hist_electron = new TH2F("hist_electron", "Electron: Event vs Waveform",nEntries, 0, nEntries, ticks, 0, ticks);
    TH2F *hist_gamma = new TH2F("hist_gamma", "Gamma: Event vs Waveform",nEntries, 0, nEntries, ticks, 0, ticks);

    for (Long64_t j = 0; j < nEntries; j++) {
        tree->GetEntry(j);

        // Loop over elements in the `electron` vector for the current event
        for (size_t i = 0; i < electron->size(); i++) {
            if (electron->at(i) == 1) { // Electron histogram

                e_event = j;
                e_om_number = om_number->at(i);
                e_energy = energy->at(i);
                e_charge = charge->at(i);

                // Fill the tree with the data
                newTree->Fill();

                // Reset the variables for the next loop iteration
                e_event = -1;
            

            } 
            else if (electron->at(i) == 0) { // Gamma histogram

                g_event = j;
                g_om_number = om_number->at(i);
                g_energy = energy->at(i);
                g_charge = charge->at(i);

                // Fill the tree with the data
                newTree->Fill();

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
