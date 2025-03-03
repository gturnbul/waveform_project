/*
Second code to run. This uses the offsets calculated in first_read_from_udd.C 
and reads them in from a csv file. It then produces the tree again, but with 
the calibrated waveform.
*/

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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// FUNCTIONS /////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Get memory cell from bin number and fcr
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024; // Now ranges from 0 to 1023
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to load CSV data into calibration_map
void open_csv_file(const std::string& filename, 
                   std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int>>>& calibration_map) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << "Error: Could not open CSV file!" << std::endl;
        return;
    }

    std::string line;
    int lineNumber = 0;

    // Skip the header line
    std::getline(file, line);

    // Read the file line by line
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string value;
        int om_index = 0, memIndex = 0;
        double calibConst = 0.0;
        int columnIndex = 0;

        while (std::getline(ss, value, ',')) {
            if (columnIndex == 0) {  // Column A: OM Number
                om_index = std::stoi(value);
            }
            if (columnIndex == 1) {  // Column B: Memory Index
                memIndex = std::stoi(value);
            }
            if (columnIndex == 3) {  // Column D: Calibration Constant
                calibConst = std::stod(value);
            }
            columnIndex++;
        }

        // Store in map using (OM Number, Memory Index) as key
        calibration_map[{om_index, memIndex}] = calibConst;
    }

    std::cout << "CSV file loaded successfully into calibration map!" << std::endl;
    file.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double process_wave_bin(int wave_bin, short wave_value, int fcr, int om_number,
                       const std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int>>>& calibration_map) {
    int memory_cell = get_cell(wave_bin, fcr); // Get memory cell

    auto key = std::make_pair(om_number, memory_cell);

    // Look for the calibration constant in the map
    auto it = calibration_map.find(key);
    if (it == calibration_map.end()) {
        std::cerr << "Warning: Calibration constant not found for OM " << om_number
                  << " and memory cell " << memory_cell << ". Using uncalibrated value." << std::endl;
        return wave_value;
    }

    // Apply calibration
    return wave_value - it->second;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

void waveform_calibration() {
    // Define the map to store calibration constants
    std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int>>> calibration_map;

    // Load CSV into the map
    open_csv_file("calibration_output_noiseRed.csv", calibration_map);

    if (calibration_map.empty()) {
        std::cout << "Error: Failed to load calibration data from CSV!" << std::endl;
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 0: Open the existing ROOT file with the original tree               ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    TFile *file = TFile::Open("filtered_gamma_fcr.root", "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Could not open input file!" << std::endl;
        return;
    }

    // Get tree
    TTree *tree = (TTree*)file->Get("Result_tree");
    if (!tree) {
        std::cout << "Error: Could not find tree in file!" << std::endl;
        file->Close();
        return;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 1: Set up variables to read the branches of the original tree       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<int> *om_number = nullptr;
    std::vector<int> *is_electron = nullptr;
    std::vector<double> *feb_energy = nullptr;
    std::vector<int> *feb_charge = nullptr;
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *feb_amplitude = nullptr;
    std::vector<int> *calo_fcr = nullptr;

    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("is_electron", &is_electron);
    tree->SetBranchAddress("energy", &feb_energy);
    tree->SetBranchAddress("charge", &feb_charge);
    tree->SetBranchAddress("wave", &wave);
    tree->SetBranchAddress("amplitude", &feb_amplitude);
    tree->SetBranchAddress("calo_fcr", &calo_fcr);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 2: Create a new ROOT file to store the new tree and new tree        ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    TFile *newFile = new TFile("waveform_gill_calibration.root", "RECREATE");
    TTree *newTree = new TTree("treebuild", "treebuild");

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 3: Set variables to hold the new branches of the new tree           ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<int> cal_om_number, cal_electron, cal_feb_charge, cal_feb_amplitude, cal_calo_fcr;
    std::vector<double> cal_feb_energy;
    std::vector<std::vector<short>> cal_wave;

    newTree->Branch("om_number", &cal_om_number);
    newTree->Branch("is_electron", &cal_electron);
    newTree->Branch("energy", &cal_feb_energy);
    newTree->Branch("charge", &cal_feb_charge);
    newTree->Branch("wave", &cal_wave);
    newTree->Branch("amplitude", &cal_feb_amplitude);
    newTree->Branch("calo_fcr", &cal_calo_fcr);

    // Global variable to track missing calibration constants
    std::set<std::pair<int, int>> missing_calibrations;
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 4: Fill the tree with data                                          ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Check for OM numbers greater than 520
        bool skipEvent = false;
        for (size_t j = 0; j < om_number->size(); j++) {
            if (om_number->at(j) > 520) {
                skipEvent = true;
                break;
            }
        }

        if (skipEvent) {
            continue;  // Skip this event if any OM number > 520
        }

        cal_om_number = *om_number;
        cal_electron = *is_electron;
        cal_feb_energy = *feb_energy;
        cal_feb_charge = *feb_charge;
        cal_feb_amplitude = *feb_amplitude;
        cal_calo_fcr = *calo_fcr;

        cal_wave.clear();

        // process waveforms
        for (size_t j = 0; j < wave->size(); j++) {
            std::vector<short> calibrated_wave;
            for (size_t k = 0; k < wave->at(j).size(); k++) {
                calibrated_wave.push_back(process_wave_bin(k, wave->at(j)[k], calo_fcr->at(j), om_number->at(j), calibration_map));
            }
            cal_wave.push_back(calibrated_wave);
        }
        ///////////////////

        newTree->Fill();
    }
    if (!missing_calibrations.empty()) {
        std::cout << "Summary of missing calibration constants:" << std::endl;
        for (const auto& entry : missing_calibrations) {
            std::cout << "OM " << entry.first << ", Memory Cell " << entry.second << " missing calibration." << std::endl;
        }
    } else {
        std::cout << "All calibration constants found successfully!" << std::endl;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Step 5: Write the tree to the file                                       ////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    newTree->Write();
    newFile->Close();
    file->Close();
}


// Necessary for the makefile      /////////////////////////////////////////////////////////////////
int main() {
    waveform_calibration();
    return 0;
}

