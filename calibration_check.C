/*
a way to compare two csv files of offsets, to check they align.
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <utility> // For std::pair
#include <boost/functional/hash.hpp>
#include <TFile.h>
#include <TTree.h>

// Function to load CSV data into a calibration map
void open_csv_file(const std::string& filename,
                   std::unordered_map<std::pair<int, int>, std::pair<double, double>, boost::hash<std::pair<int, int>>>& calibration_map,
                   int om_col, int mem_col, int calib_col, int redchi_col, bool shift_memory = false) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << "Error: Could not open CSV file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string value;
        int temp_om = 0, temp_mem = 0;
        double temp_calib = 0.0, temp_red_chi2 = 0.0;
        int columnIndex = 0;

        while (std::getline(ss, value, ',')) {
            if (columnIndex == om_col) temp_om = std::stoi(value);
            if (columnIndex == mem_col) temp_mem = std::stoi(value);
            if (columnIndex == calib_col) temp_calib = std::stod(value);
            if (columnIndex == redchi_col) temp_red_chi2 = std::stod(value);
            columnIndex++;
        }

        if (shift_memory) temp_mem -= 1;  // Shift Memory Index from 1-1024 to 0-1023

        if (temp_om < 0 || temp_mem < 0) {
            std::cout << "Skipping invalid row: OM=" << temp_om << ", Mem=" << temp_mem << std::endl;
            continue;
        }

        calibration_map[{temp_om, temp_mem}] = {temp_calib, temp_red_chi2};
    }

    std::cout << "CSV file loaded successfully: " << filename << std::endl;
    file.close();
}

// Function to save the data from the calibration maps to a ROOT file
void save_to_root(const std::unordered_map<std::pair<int, int>, std::pair<double, double>, boost::hash<std::pair<int, int>>>& map1,
                  const std::unordered_map<std::pair<int, int>, std::pair<double, double>, boost::hash<std::pair<int, int>>>& map2,
                  const std::string& root_filename) {
    // Create a ROOT file
    TFile root_file(root_filename.c_str(), "RECREATE");
    if (!root_file.IsOpen()) {
        std::cout << "Error: Could not open ROOT file!" << std::endl;
        return;
    }

    // Create a ROOT tree
    TTree tree("CalData", "Raw Calibration Data");

    int om1, mem1, om2, mem2;
    double calib1, redchi1, calib2, redchi2, diff_calib;

    // Define branches for both sets of data
    tree.Branch("OM1", &om1, "OM1/I");
    tree.Branch("MemCell1", &mem1, "MemoryIndex1/I");
    tree.Branch("Calib1", &calib1, "CalibrationConstant1/D");
    tree.Branch("RedChi1", &redchi1, "ReducedChi21/D");

    tree.Branch("OM2", &om2, "OM2/I");
    tree.Branch("MemCell2", &mem2, "MemoryIndex2/I");
    tree.Branch("Calib2", &calib2, "CalibrationConstant2/D");
    tree.Branch("RedChi2", &redchi2, "ReducedChi22/D");

    tree.Branch("CalibDiff", &diff_calib, "CalibrationDifference/D");

    // Populate the tree with data from both maps
    for (const auto& [key, value1] : map1) {
        om1 = key.first;
        mem1 = key.second;
        calib1 = value1.first;
        redchi1 = value1.second;

        // Look up the same key in map2
        auto it2 = map2.find(key);
        if (it2 != map2.end()) {
            om2 = it2->first.first;
            mem2 = it2->first.second;
            calib2 = it2->second.first;
            redchi2 = it2->second.second;
        } else {
            om2 = om1;  // Maintain alignment
            mem2 = mem1;
            calib2 = -999.0;  // Indicate missing data
            redchi2 = -999.0;
        }

        // Calculate the difference between the calibration constants
        diff_calib = std::abs(calib1 - calib2); 

        // Fill the tree
        tree.Fill();
    }

    // Save the tree to the ROOT file
    root_file.Write();
    root_file.Close();
    std::cout << "Data saved to ROOT file: " << root_filename << std::endl;
}

int main() {
    std::unordered_map<std::pair<int, int>, std::pair<double, double>, boost::hash<std::pair<int, int>>> calibration_map1;
    std::unordered_map<std::pair<int, int>, std::pair<double, double>, boost::hash<std::pair<int, int>>> calibration_map2;

    // // Load File 1 (Memory Index 1-1024 to 0-1023, shift_memory = true)
    // open_csv_file("calibration_output_get_cell_test.csv", calibration_map1, 2, 3, 4, 6, true);  

    // Load File 1 (Memory Index 0-1023, no shift needed)
    open_csv_file("calibration_output_get_cell_test.csv", calibration_map2, 0, 1, 2, 4, false);
    
    // Load File 2 (Memory Index 0-1023, no shift needed)
    open_csv_file("calibration_output_reduced.csv", calibration_map2, 0, 1, 2, 4, false);

    // Save to ROOT file
    save_to_root(calibration_map1, calibration_map2, "calibration_data_raw.root");

    return 0;
}
