#include <vector>
#include <string>
#include <iostream>
#include <TFile.h>
#include <TTree.h>

void prefixtest() {
    // Vector where true represents "e_" prefix and false represents "g_" prefix
    std::vector<bool> electron_indicator = {true, false};

    // Corresponding values for each attribute
    std::vector<double> energy = {10.5, 20.2};  // e_ and g_ respectively
    std::vector<double> charge = {1.1, 2.2};    // e_ and g_ respectively
    std::vector<int> om_number = {101, 102};    // e_ and g_ respectively
    std::vector<int> event = {201, 202};        // e_ and g_ respectively

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////       Create a new ROOT file to store the new tree and new tree        ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TFile *file = new TFile("prefixtest.root", "RECREATE");

    // Create a TTree to store the data
    TTree *tree = new TTree("data", "Data with electron and gamma variables");

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////       Set variables to hold the new branches of the new tree           ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    double energy_val;
    double charge_val;
    int om_number_val;
    int event_val;
    std::string prefix;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////       Set the branch addresses                                         ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    tree->Branch("prefix", &prefix);
    tree->Branch("energy", &energy_val);
    tree->Branch("charge", &charge_val);
    tree->Branch("om_number", &om_number_val);
    tree->Branch("event", &event_val);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////       Fill the tree with data                                          ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
   
    // Loop over the data and fill the tree
    for (int i = 0; i < electron_indicator.size(); i++) {
        // Determine the prefix based on the value in electron_indicator
        prefix = electron_indicator[i] ? "e_" : "g_";

        // Get corresponding values
        energy_val = energy[i];
        charge_val = charge[i];
        om_number_val = om_number[i];
        event_val = event[i];

        // Fill the tree with the data
        tree->Fill();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////       Write the tree to the file                                       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    tree->Write();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////       Close the file to save everything                                ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    file->Close();
}
============================ included for later, but not needed now =============================
// Main function to run the test
int main() {
    prefixtest();
    return 0;
}
