#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>   // ROOT file handling
#include <TTree.h>   // ROOT tree structure

void prefixtest() {

    // Initialize vectors to hold the data with 8 entries each (4 for e_ and 4 for g_)
    std::vector<double> energy = {10.5, 20.2, 30.5, 40.3, 15.7, 25.4, 35.1, 45.6};  // e_ and g_ respectively
    std::vector<double> charge = {1.1, 2.2, 3.3, 4.4, 1.5, 2.6, 3.7, 4.8};    // e_ and g_ respectively
    std::vector<int> om_number = {101, 102, 103, 104, 105, 106, 107, 108};    // e_ and g_ respectively
    std::vector<int> event = {201, 201, 202, 202, 203, 203, 204, 204};        // e_ and g_ respectively

    // Vector to store the electron indicator (true for e_, false for g_)
    std::vector<bool> electron_indicator;

    // Open the input file to read the boolean indicators
    std::ifstream input_file("prefix_e_vector.txt");

    // Read the 1D vector from the file (1 or 0)
    int value;
    while (input_file >> value) {
        if (value == 1 || value == 0) {
            electron_indicator.push_back(value);
        } else {
            std::cerr << "Invalid value in input file. Only 1 and 0 are allowed." << std::endl;
            return;
        }
    }

    // Ensure the size of electron_indicator matches the size of other data vectors
    if (electron_indicator.size() != energy.size()) {
        std::cerr << "Mismatch in data size between electron_indicator and other data vectors." << std::endl;
        return;
    }

    // Open the ROOT file where the data will be saved
    TFile* newFile = new TFile("prefixtest.root", "RECREATE");

    // Create a TTree to store the data
    TTree* newTree = new TTree("Result_tree", "Tree with electron and gamma data");

    // Declare variables to hold data for each branch
    double energy_var; 
    double charge_var;
    int om_number_var;
    double event_var;

    // Create branches for electron data (e_)
    newTree->Branch("e_energy", &energy_var, "e_energy/D");
    newTree->Branch("e_charge", &charge_var, "e_charge/D");
    newTree->Branch("e_om_number", &om_number_var, "e_om_number/I");
    newTree->Branch("e_event", &event_var, "e_event/I");

    // Create branches for gamma data (g_)
    newTree->Branch("g_energy", &energy_var, "g_energy/D");
    newTree->Branch("g_charge", &charge_var, "g_charge/D");
    newTree->Branch("g_om_number", &om_number_var, "g_om_number/I");
    newTree->Branch("g_event", &event_var, "g_event/I");

    // Loop through the data and fill the TTree
    for (int i = 0; i < electron_indicator.size(); i++) {
        // Based on the prefix (e_ or g_), assign data to the corresponding variables
        if (electron_indicator[i]==1) {
            // For electron (e_)
            energy_var = energy[i];
            charge_var = charge[i];
            om_number_var = om_number[i];
            event_var = event[i];

            // Fill the branches for e_
            newTree->Fill();
        } else if (electron_indicator[i]==0) {
            // For gamma (g_)
            energy_var = energy[i];
            charge_var = charge[i];
            om_number_var = om_number[i];
            event_var = event[i];

            // Fill the branches for g_
            newTree->Fill();
        }
    }

    // Write the tree to the file
    newTree->Write();

    // Close the ROOT file
    newFile->Close();
}
============================ included for later, but not needed now =============================
// Main function to run the test
int main() {
    prefixtest();
    return 0;
}
