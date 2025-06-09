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
#include "TMultiGraph.h"

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

// Function to calculate timestamp difference per OM
double calculate_timestamp_diff(int om_num, double current_timestamp) {
    static std::unordered_map<int, double> last_timestamp_per_om;

    double timestamp_diff = -1;
    if (last_timestamp_per_om.find(om_num) != last_timestamp_per_om.end()) {
        timestamp_diff = current_timestamp - last_timestamp_per_om[om_num];
    }
    last_timestamp_per_om[om_num] = current_timestamp;

    return timestamp_diff;
}

// Function to plot waveform
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

    TString filename = Form("eom_elimination_run%d_event%d_om%d.png", 1143, event, om_num); // run number is 0, change if needed
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
double calculate_stddev(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 0; i < 976; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 976);
}

// Function to calculate error on the mean (EOM), with sigma calculated inside
double calculate_eom(const std::vector<short>& waveform, double baseline) {
    int N = waveform.size();
    double sum_sq_diff = 0.0;

    for (int i = 0; i < N; ++i) {
        double diff = waveform[i] - baseline;
        sum_sq_diff += diff * diff;
    }

    double variance = sum_sq_diff / N;      // Population variance
    double sigma = std::sqrt(variance);     // Standard deviation
    double eom = sigma / std::sqrt(N);      // Error on the mean

    return eom;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Get memory cell from bin number and FCR
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024; // Now ranges from 0 to 1023
}

std::vector<short> reorder_waveform(const std::vector<short>& waveform, int fcr) {
    std::vector<short> reordered(1024, -999);  // Use -999 as a null marker

    int max_bins = std::min((int)waveform.size(), 976);  // eliminate end of waveform spike
    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);  // [0, 1023]
        reordered[reordered_index] = waveform[i];
    }

    return reordered;
}
// Clean and reorder the waveform
std::vector<short> clean_and_reorder(const std::vector<short>& waveform, int fcr) {
    std::vector<short> cleaned = waveform;

    for (int i = 0; i < 48; ++i) {
        cleaned[1024 - 48 + i] = cleaned[i];
    }

    std::vector<short> reordered;
    reordered.reserve(1024);
    for (int i = 0; i < 1024; ++i) {
        reordered.push_back(cleaned[get_cell(i, fcr)]);
    }

    return reordered;
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
std::vector<int> *timestamp = nullptr;
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
tree->SetBranchAddress("digicalo.timestamp", &timestamp);

// Get the number of entries in the tree
int max_entries = tree->GetEntries();

// Create output ROOT file and TTree for baseline data
TFile *outfile = new TFile("baseline_output_tree1143.root", "RECREATE");
TTree *baseline_tree = new TTree("baseline_tree", "OM baseline data");


// initialise output variables
int om_num_out = -1;
double baseline_out = 0;
double stddev_out = 0;
int event_num = -1;
double eom = 0;
double timestamp_out = 0;
double timestamp_diff = -1;
std::vector<int> tail_bin_indices;
std::vector<double> tail_bin_distances;


// Create branches in the output tree
baseline_tree->Branch("event_num", &event_num, "event_num/I");
baseline_tree->Branch("om_num", &om_num_out, "om_num/I");
baseline_tree->Branch("baseline", &baseline_out, "baseline/D");
baseline_tree->Branch("stddev", &stddev_out, "stddev/D");
baseline_tree->Branch("eom", &eom, "eom/D");
baseline_tree->Branch("timestamp", &timestamp_out, "timestamp/D");
baseline_tree->Branch("timestamp_diff", &timestamp_diff, "timestamp_diff/D");
baseline_tree->Branch("tail_bin_indices", &tail_bin_indices);
baseline_tree->Branch("tail_bin_distances", &tail_bin_distances);



// Create a histogram for the baseline distribution
TH1D *baseline_hist = new TH1D("baseline_hist", "Baseline Distribution;Baseline [ADC];Entries", 100, 0, 4000);
std::unordered_map<int, double> last_timestamp_per_om;

// loop over the entries in the tree
for (int event = 0; event < max_entries; ++event) {
    tree->GetEntry(event);

    for (int k = 0; k < calo_nohits; ++k) {
        int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);

        std::vector<short>& wave_k = wave->at(k);
        double baseline = calculate_baseline(wave_k);

        ///////////////////////////////////////////////////////////////////
        tail_bin_indices.clear();
        tail_bin_distances.clear();

        int tail_start = 900;
        int tail_end = std::min((int)wave_k.size(), 1024);  // safeguard

        for (int i = tail_start; i < tail_end; ++i) {
            tail_bin_indices.push_back(i);
            tail_bin_distances.push_back(wave_k[i] - baseline);
        }
        ////////////////////////////////////////////////////////////////////


        double stddev = calculate_stddev(wave_k, baseline);
        eom = calculate_eom(wave_k, baseline);

        if (eom >= 0.075 && eom <= 0.076) {
            std::cout << "Event " << event << ", OM " << om_num << ", EOM = " << eom << std::endl;
            PlotWaveform(event, om_num, wave_k);
        }
    

        om_num_out = om_num;
        baseline_out = baseline;
        stddev_out = stddev;
        event_num = event;
        timestamp_out = timestamp->at(k)* 6.25e-9; // Convert to seconds;
        timestamp_diff = calculate_timestamp_diff(om_num, timestamp_out);

        baseline_tree->Fill();  

        }
}
//output the number of entries to the terminal
cout << "Max entries: " << max_entries << "\n";

// Write and close output
baseline_tree->Write();
outfile->Close();
delete outfile;

file->Close();
}
