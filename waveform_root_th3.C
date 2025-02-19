#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <memory>
#include <cmath>         // For math functions (with _USE_MATH_DEFINES)
#include <sys/stat.h>    // For file system operations
#include <sys/types.h>   // For file system types
#include <cassert>       // For debugging assertions

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

// ROOT Random Number Generator
#include <TRandom3.h>

// ROOT Fitting Libraries
#include <TF1.h>

// Additional ROOT Utilities
#include <TPaveText.h>
#include <TArrayF.h>

// Macro Definition
#define _USE_MATH_DEFINES

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// FUNCTIONS ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// Get memory cell from bin number and FCR
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024; // Now ranges from 0 to 1023
}

////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// MAIN CODE ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char const *argv[]) {
    // Get run number and OM range from command-line arguments
    std::string run = "";
    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "-r" && i + 1 < argc) {
            run = std::string(argv[i + 1]);
        }
    }

    if (run.empty()) {
        std::cerr << "Error: No run number provided! Use -r <run_number>" << std::endl;
        return 1;
    }

    // Open the ROOT file
    TFile *file = new TFile(Form("snemo_run-%s_udd.root", run.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file snemo_run-" << run << "_udd.root" << std::endl;
        return 1;
    }

    // Get TTree
    TTree *tree = (TTree *)file->Get("SimData");
    if (!tree) {
        std::cerr << "Error: TTree 'SimData' not found in file!" << std::endl;
        file->Close();
        return 1;
    }

    // Variables for branches
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *calo_wall = nullptr;
    std::vector<int> *calo_side = nullptr;
    std::vector<int> *calo_column = nullptr;
    std::vector<int> *calo_row = nullptr;
    std::vector<int> *calo_type = nullptr;
    std::vector<int> *fcr = nullptr;
    int calo_nohits = 0;

    // Check if branches exist before setting addresses
    if (!tree->GetBranch("digicalo.nohits") || 
        !tree->GetBranch("digicalo.waveform") ||
        !tree->GetBranch("digicalo.wall") ||
        !tree->GetBranch("digicalo.side") ||
        !tree->GetBranch("digicalo.column") ||
        !tree->GetBranch("digicalo.row") ||
        !tree->GetBranch("digicalo.type") ||
        !tree->GetBranch("digicalo.fcr")) {
        std::cerr << "Error: One or more required branches are missing!" << std::endl;
        file->Close();
        return 1;
    }

    // Set branch addresses
    tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
    tree->SetBranchAddress("digicalo.waveform", &wave);
    tree->SetBranchAddress("digicalo.wall", &calo_wall);
    tree->SetBranchAddress("digicalo.side", &calo_side);
    tree->SetBranchAddress("digicalo.column", &calo_column);
    tree->SetBranchAddress("digicalo.row", &calo_row);
    tree->SetBranchAddress("digicalo.type", &calo_type);
    tree->SetBranchAddress("digicalo.fcr", &fcr);

    // Define histogram
    TH3D *hist = new TH3D("waveform_hist", "OM vs Memory Cell vs Waveform;OM Number;Memory Cell;Waveform Value", 
                           520, 0, 519,   
                           1024, 0, 1023, 
                           235, 3475, 3710);

    // Create output file
    TFile *outfile = new TFile("waveform_data_TH3v1.root", "RECREATE");
    hist->Reset();

    // Counter for total occurrences
    int total_occurrences = 0;

    int max_entries = tree->GetEntries();
    for (int event = 0; event < max_entries; ++event) {
        tree->GetEntry(event);

        if (calo_nohits <= 0) continue;  // Skip events with no hits

        for (int k = 0; k < calo_nohits; ++k) {
            // Validate vector sizes before accessing
            if (k >= wave->size() || k >= calo_type->size() || k >= fcr->size()) {
                std::cerr << "Warning: Index out of range in event " << event << std::endl;
                continue;
            }

            // Calculate OM number
            int om_num;
            if (calo_type->at(k) == 0) {
                om_num = calo_side->at(k) * 260 + calo_column->at(k) * 13 + calo_row->at(k);
            } else if (calo_type->at(k) == 1) {
                om_num = 520 + calo_side->at(k) * 64 + calo_wall->at(k) * 32 +
                        calo_column->at(k) * 16 + calo_row->at(k);
            } else if (calo_type->at(k) == 2) {
                om_num = 648 + calo_side->at(k) * 32 + calo_wall->at(k) * 16 + calo_column->at(k);
            } else {
                continue;
            }

            int fcr_value = fcr->at(k);
            const std::vector<short> &wave_k = wave->at(k);

            if (wave_k.empty()) continue;

            // Track whether we've seen this OM & memory cell in the current event
            std::set<int> seen_events;

            for (int bin_no = 0; bin_no < wave_k.size(); ++bin_no) {
                short wave_value = wave_k[bin_no];
                if (wave_value > 3000) {
                    int mem_cell_no = get_cell(bin_no, fcr_value);

                    if (om_num == 64 && mem_cell_no == 799) {
                        total_occurrences++;  // Count occurrences across all events

                        // Check if we've seen this event before
                        if (seen_events.count(event)) {
                            std::cout << "WARNING: Duplicate entry in Event " << event 
                                    << " for OM 64, Mem Cell 799!" << std::endl;
                        } else {
                            seen_events.insert(event);  // Mark this event as seen
                        }

                        hist->Fill(om_num, mem_cell_no, wave_value);
                    }
                }
            }
        }

        if (event % 1000 == 0) {
            outfile->Flush();
        }
    }

    hist->Write();
    outfile->Close();
    file->Close();

    // Print the final total occurrences at the end
    std::cout << "Total occurrences of OM 64, Mem Cell 799: " << total_occurrences << std::endl;

    return 0;
}
