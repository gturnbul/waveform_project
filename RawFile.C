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

//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Functions //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

// Calculate the OM number from the calo udd variables
int calculate_om_num(std::vector<int> *calo_type, std::vector<int> *calo_side, 
                     std::vector<int> *calo_wall, std::vector<int> *calo_column, 
                     std::vector<int> *calo_row, int k) {
    int om_num = -1;
    if (calo_type->at(k) == 0) {
        om_num = calo_side->at(k) * 260 + calo_column->at(k) * 13 + calo_row->at(k);
    } else if (calo_type->at(k) == 1) {
        om_num = 520 + calo_side->at(k) * 64 + calo_wall->at(k) * 32 + calo_column->at(k) * 16 + calo_row->at(k);
    } else if (calo_type->at(k) == 2) {
        om_num = 648 + calo_side->at(k) * 32 + calo_wall->at(k) * 16 + calo_column->at(k);
    }
    return om_num;
}

void PlotWaveform(int event, int om_num, std::vector<short>& wave_k) {
    if (wave_k.empty()) {
        std::cout << "Warning: empty waveform for event " << event << " om_num " << om_num << std::endl;
        return; // skip drawing empty waveform
    }

    TCanvas *c1 = new TCanvas("c1", "Waveform", 800, 600);
    TGraph *graph = new TGraph(wave_k.size());

    // Fill the TGraph with waveform data: x = index, y = ADC count
    for (size_t i = 0; i < wave_k.size(); ++i) {
        graph->SetPoint(i, i, wave_k[i]);
    }

    graph->SetTitle(Form("Event %d - OM %d;Clock Ticks;ADC Counts", event, om_num));
    graph->GetXaxis()->SetLimits(0, wave_k.size() - 1);
    graph->SetLineColor(kBlue);
    graph->Draw("AL");

    TString filename = Form("waveform_reordered_run%d_event%d_om%d.png", 0, event, om_num); // run number is 0, change if needed
    c1->SaveAs(filename);

    delete graph;
    delete c1;
}
////////////////////////////// FUNCTIONS from the WAVEFORMS ////////////////////////////////////
// Function to calculate baseline
double calculate_baseline(const std::vector<short>& waveform) {
    double baseline = 0;
    for (int i = 0; i < 976; ++i) {
        baseline += waveform[i];
    }
    return baseline / 976;
}

// Function to calculate chi2df
double calculate_chi2ndf(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 0; i < 976; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 976);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Get memory cell from bin number and FCR
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024; // Now ranges from 0 to 1023
}


//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Main ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
int main(){
    // Open the ROOT file
TFile *file = new TFile("snemo_run-1143_udd.root", "READ");
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

// Get the number of entries in the tree
int max_entries = tree->GetEntries();

// Create output ROOT file and TTree for baseline data
TFile *outfile = new TFile("baseline_output_tree.root", "RECREATE");
TTree *baseline_tree = new TTree("baseline_tree", "OM baseline data");

// initialise output variables
int om_num_out = -1;
double baseline_out = 0;
double chi2ndf_out = 0;
double max_baseline = -1;
int event_num = -1;
int count_chi2ndf_gt3 = 0;
std::set<int> valid_om_nums; // stores unique OMs with valid waveform


// Create branches in the output tree
baseline_tree->Branch("om_num", &om_num_out, "om_num/I");
baseline_tree->Branch("baseline", &baseline_out, "baseline/D");
baseline_tree->Branch("chi2ndf", &chi2ndf_out, "chi2ndf/D");
baseline_tree->Branch("event_num", &event_num, "event_num/I");

// //Define the range of events to process
// std:set<int> selected_events = {0, 1, 2, 3, 4, 5}; // Example: process events 0 to 5
// // Define your allowed OMs:
// std::set<int> allowed_oms = {64,65,66};  // put the OM numbers you want here

// Create a histogram for the baseline distribution
TH1D *baseline_hist = new TH1D("baseline_hist", "Baseline Distribution;Baseline [ADC];Entries", 100, 0, 4000);

// loop over the entries in the tree
for (int event = 0; event < max_entries; ++event) {
    // Skip events not in the selected set
    // if (selected_events.find(event) == selected_events.end()) {
    //     continue;
    // }

    tree->GetEntry(event);

    for (int k = 0; k < calo_nohits; ++k) {
        int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);

        // // Check if om_num is allowed, skip if not
        // if (allowed_oms.find(om_num) == allowed_oms.end()) {
        //     continue;
        // }

        std::vector<short>& wave_k = wave->at(k);
        double baseline = calculate_baseline(wave_k);
        double chi2ndf = calculate_chi2ndf(wave_k, baseline);

        // //print specific waveform
        // if (event == 5 && om_num == 123 && chi2ndf < 3) {
        //     PlotWaveform(event, om_num, wave_k);  // Save the waveform as a .png
        // }
        if (event == 0 && om_num == 123 && chi2ndf < 3) {
            std::vector<short> reordered;
            reordered.reserve(wave_k.size());
            for (size_t i = 0; i < wave_k.size(); ++i) {
                reordered.push_back(wave_k[get_cell(i, fcr->at(k))]);
            }
            PlotWaveform(event, om_num, reordered);  // plots reordered waveform
        }

        // Only consider waveforms with chi2ndf < 3 and baseline > 3000
        if (chi2ndf < 3 && baseline > 3000) {
            // Update max baseline
            if (baseline > max_baseline) {
                max_baseline = baseline;
            }

            valid_om_nums.insert(om_num); // Store unique OM numbers with valid waveforms

            // Fill histogram and output tree
            baseline_hist->Fill(baseline);

            om_num_out = om_num;
            baseline_out = baseline;
            chi2ndf_out = chi2ndf;
            event_num = event;

            baseline_tree->Fill();
        }

    }
}
//output the number of entries to the terminal
cout << "Max entries: " << max_entries << "\n";
cout << "Max baseline: " << max_baseline << "\n";
//cout << "Number of events with chi2ndf > 3: " << count_chi2ndf_gt3 << "\n";
cout << "Number of unique OM numbers with valid waveforms: " << valid_om_nums.size() << "\n";

// // Save histogram to a ROOT file
// TFile *outfile = new TFile("baseline_output_all.root", "RECREATE");
// baseline_hist->Write();

// Write and close output
baseline_tree->Write();
outfile->Close();
delete outfile;
delete baseline_hist;

file->Close();
}
