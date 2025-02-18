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
    int om_min = -1;  // Default: process all OMs
    int om_max = -1;  // Default: process all OMs

    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "-r") {
            run = std::string(argv[i + 1]);
        }
        if (std::string(argv[i]) == "-om_min") {
            om_min = std::stoi(argv[i + 1]);
        }
        if (std::string(argv[i]) == "-om_max") {
            om_max = std::stoi(argv[i + 1]);
        }
    }

    // Open the ROOT file
    TFile *file = new TFile(Form("snemo_run-%s_udd.root", run.c_str()), "READ");
    TTree *tree = (TTree *)file->Get("SimData");

    // Variables for branches
    std::vector<std::vector<short>> *wave = nullptr;
    std::vector<int> *calo_wall = nullptr;
    std::vector<int> *calo_side = nullptr;
    std::vector<int> *calo_column = nullptr;
    std::vector<int> *calo_row = nullptr;
    std::vector<int> *calo_type = nullptr;
    std::vector<int> *fcr = nullptr;
    int calo_nohits = 0;

    // Set branch addresses
    tree->SetBranchAddress("digicalo.nohits", &calo_nohits);
    tree->SetBranchAddress("digicalo.waveform", &wave);
    tree->SetBranchAddress("digicalo.wall", &calo_wall);
    tree->SetBranchAddress("digicalo.side", &calo_side);
    tree->SetBranchAddress("digicalo.column", &calo_column);
    tree->SetBranchAddress("digicalo.row", &calo_row);
    tree->SetBranchAddress("digicalo.type", &calo_type);
    tree->SetBranchAddress("digicalo.fcr", &fcr);

    // CSV output file
    ofstream csv_file("calibration_output_vTEST1.csv");
    csv_file << "om_num,mem_cell_no,mean,sigma,reduced_chi2\n"; // CSV header

    // Loop through events
    int max_entries = tree->GetEntries();
    map<int, map<int, TH1D*>> histograms; // Nested map: om_num -> mem_cell_no -> histogram

    for (int event = 0; event < max_entries; ++event) {
        tree->GetEntry(event);

        for (int k = 0; k < calo_nohits; ++k) {
            // Calculate the OM number
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

            // Skip this OM if it is outside the specified range
            if ((om_min != -1 && om_num < om_min) || (om_max != -1 && om_num > om_max)) {
                continue;
            }

            // Extract FCR and waveform
            int fcr_value = fcr->at(k);
            std::vector<short> wave_k = wave->at(k);

            if (wave_k.empty()) {
                std::cout << "Skipping empty waveform for OM: " << om_num << std::endl;
                continue;
            }

            // Fill histograms for each memory cell
            for (int bin_no = 0; bin_no < wave_k.size(); ++bin_no) {
                int mem_cell_no = get_cell(bin_no, fcr_value);
                short wave_value = wave_k[bin_no];

                // Outlier suppression
                if (wave_value > 3300 && wave_value < 3750) {
                    // Create histogram for this OM and memory cell if not already done
                    if (histograms[om_num][mem_cell_no] == nullptr) {
                        // Calculate min and max values for the waveform
                        auto min_max = std::minmax_element(wave_k.begin(), wave_k.end());
                        double min_wave = *min_max.first;
                        double max_wave = *min_max.second;
                        int num_bins = static_cast<int>(max_wave - min_wave + 1);

                        histograms[om_num][mem_cell_no] = new TH1D(
                            Form("h_om%d_cell%d_fcr%d", om_num, mem_cell_no, fcr_value),
                            Form("Waveform: OM %d, Cell %d, FCR %d", om_num, mem_cell_no, fcr_value),
                            num_bins, min_wave, max_wave);
                    }

                    // Fill the histogram
                    histograms[om_num][mem_cell_no]->Fill(wave_value);
                }
            }
        }
    }

    // Perform Gaussian fits and write results to CSV
    for (const auto &om_entry : histograms) {
        int om_num = om_entry.first;
        std::cout << "Processing OM: " << om_num << std::endl;

        for (const auto &cell_entry : om_entry.second) {
            int mem_cell_no = cell_entry.first;
            TH1D *hist = cell_entry.second;

            if (hist->GetEntries() <= 0) {
                std::cout << "Skipping empty histogram for OM: " << om_num << ", Cell: " << mem_cell_no << std::endl;
                delete hist;
                continue;
            }

            // Check bins for non-zero content
            bool hasNonZeroBins = false;
            for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
                if (hist->GetBinContent(bin) > 0) {
                    hasNonZeroBins = true;
                    break;
                }
            }

            if (!hasNonZeroBins) {
                std::cout << "Skipping histogram with zero-filled bins for OM: " << om_num << ", Cell: " << mem_cell_no << std::endl;
                delete hist;
                continue;
            }

            // Fit the histogram
            if (hist->GetEntries() > 10) {
            //std::cout << "Fitting histogram for OM: " << om_num << std::endl;

            // // Determine the range of the fit using the minimum and maximum values in the histogram
            // double gaus_min = hist->GetXaxis()->GetXmin(); // Lowest x-axis value
            // double gaus_max = hist->GetXaxis()->GetXmax(); // Highest x-axis value

            TF1 *gaus_fit = new TF1("gaus_fit", "gaus");
            hist->Fit(gaus_fit, "Q"); // Quiet mode fit

            double mean = gaus_fit->GetParameter(1);
            double sigma = gaus_fit->GetParameter(2);
            double chi2 = gaus_fit->GetChisquare();
            double ndof = gaus_fit->GetNDF();
            double red_chi2 = chi2 / ndof;

            // Write to CSV file
            csv_file << om_num << "," << mem_cell_no << "," << mean << "," << sigma << "," << red_chi2 << "\n";
            delete gaus_fit;
            }

            delete hist;
        }
    }

    // Close the CSV file and ROOT file
    csv_file.close();
    file->Close();

    return 0;
}
