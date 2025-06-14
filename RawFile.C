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
#include "TSystem.h"
bool overlayWritten = false;
void dbg(const char* msg) { std::cout << msg << std::endl; }



const short EMPTY_BIN = -1;
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
void PlotWaveform(int event, int om_num, const std::vector<short>& wave_k) {
    if (wave_k.empty()) {
        std::cout << "Warning: empty waveform for event " << event << " om_num " << om_num << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "Waveform", 800, 600);
    TGraph *graph = new TGraph();

    int point_index = 0;
    for (size_t i = 0; i < wave_k.size(); ++i) {
        if (wave_k[i] != EMPTY_BIN) {  // skip empty bins
            graph->SetPoint(point_index, i, wave_k[i]);
            point_index++;
        }
    }

    if (point_index == 0) {
        std::cout << "Warning: waveform for event " << event << " om_num " << om_num << " contains only empty bins." << std::endl;
        delete graph;
        delete c1;
        return;
    }

    graph->SetTitle(Form("Event %d - OM %d;Clock Ticks;ADC Counts", event, om_num));
    graph->GetXaxis()->SetLimits(0, wave_k.size() - 1);
    graph->SetLineColor(kBlue);
    graph->Draw("AL");

    TString filename = Form("reordered_run%d_event%d_om%d.png", 1143, event, om_num);
    c1->SaveAs(filename);

    delete graph;
    delete c1;
}
// Function to plot overlayed waveform
void PlotOverlayWaveform(int event, int om_num,
                         const std::vector<short>& reordered)
{
    static const int GROUP_SIZE = 16;      // cells per overlay
    static const int GROUPS     = 4;       // 0‑15, 16‑31, 32‑47, 48‑63
    static const Color_t COLORS[GROUPS] =
        { kRed+1, kBlue+1, kGreen+2, kMagenta+1 };

    // Guard: make sure we *have* those memory cells
    if (reordered.size() < GROUP_SIZE * GROUPS) {
        std::cerr << "Waveform too short for overlay\n";
        return;
    }

    TCanvas    *c  = new TCanvas("c_overlay","Overlayed waveform",800,600);
    TMultiGraph mg;

    for (int g = 0; g < GROUPS; ++g) {
        std::unique_ptr<TGraph> gr(new TGraph());
        gr->SetLineColor(COLORS[g]);
        gr->SetLineWidth(2);

        for (int i = 0; i < GROUP_SIZE; ++i) {
            int cell = g * GROUP_SIZE + i;        // absolute mem‑cell index
            if (reordered[cell] != EMPTY_BIN) {
                gr->SetPoint(gr->GetN(), i, reordered[cell]); // x = 0‑15
            }
        }
        mg.Add(gr.release()); // TMultiGraph now owns the graph
    }

    mg.SetTitle(Form("Event %d - OM %d;Memory-cell index (0-15);ADC counts",
                     event, om_num));
    mg.Draw("AL");

    // Simple legend so you know which colour is which block
    const char *labels[GROUPS] = { "cells 0-15",
                                   "cells 16-31",
                                   "cells 32-47",
                                   "cells 48-63" };
    TLegend leg(0.70,0.70,0.90,0.88);
    for (int g = 0; g < GROUPS; ++g)
        leg.AddEntry(mg.GetListOfGraphs()->At(g), labels[g], "l");
    leg.Draw();

    c->SaveAs(Form("mem_cell_overlay_event%d_om%d.png", event, om_num));
    delete c; // graphs are owned by mg/canvas and get cleaned up
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
    double sum_sq_diff = 0.0;

    for (int i = 0; i < 976; ++i) {
        double diff = waveform[i] - baseline;
        sum_sq_diff += diff * diff;
    }

    double variance = sum_sq_diff / 976;      // Population variance
    double sigma = std::sqrt(variance);     // Standard deviation
    double eom = sigma / std::sqrt(976);      // Error on the mean

    return eom;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Get memory cell from bin number and FCR
int get_cell(int bin_no, int fcr) {
    return (bin_no + fcr) % 1024; // Now ranges from 0 to 1023
}

std::vector<short> reorder_waveform(const std::vector<short>& waveform, int fcr) {
    std::vector<short> reordered(1024, EMPTY_BIN);  // initialize all bins as empty

    int max_bins = std::min((int)waveform.size(), 976);  // only use bins 0 to 975

    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);
        reordered[reordered_index] = waveform[i];
    }

    // bins 976 to 1023 remain EMPTY_BIN, indicating gaps

    return reordered;
}


//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Main ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
int main(){
gROOT->SetBatch(kTRUE);   // write canvases without opening X‑windows
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


// Create branches in the output tree
baseline_tree->Branch("event_num", &event_num, "event_num/I");
baseline_tree->Branch("om_num", &om_num_out, "om_num/I");
baseline_tree->Branch("baseline", &baseline_out, "baseline/D");
baseline_tree->Branch("stddev", &stddev_out, "stddev/D");
baseline_tree->Branch("eom", &eom, "eom/D");
baseline_tree->Branch("timestamp", &timestamp_out, "timestamp/D");
baseline_tree->Branch("timestamp_diff", &timestamp_diff, "timestamp_diff/D");


// Create a histogram for the baseline distribution
TH1D *baseline_hist = new TH1D("baseline_hist", "Baseline Distribution;Baseline [ADC];Entries", 100, 0, 4000);
std::unordered_map<int, double> last_timestamp_per_om;

const int  TARGET_EVENT = 25;  // <‑‑ pick an event ≥ 9
const int  TARGET_OM     = 428; // your OM of interest

// loop over the entries in the tree
for (int event = 0; event < max_entries; ++event) {
    tree->GetEntry(event);

    for (int k = 0; k < calo_nohits; ++k) {
        int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);

        std::vector<short>& wave_k = wave->at(k);
        int fcr_k = fcr->at(k);

        // Reorder waveform (this is just reordering, not preprocessing)
        std::vector<short> r_waveform = reorder_waveform(wave_k, fcr_k); // or clean_and_reorder()

        double baseline = calculate_baseline(wave_k);
        double stddev = calculate_stddev(wave_k, baseline);
        eom = calculate_eom(wave_k, baseline);

        ///////////////////////////////////////////////////////////////////////////////////////
        if (event == TARGET_EVENT && om_num == TARGET_OM) {
            if (!overlayWritten && event == TARGET_EVENT && om_num == TARGET_OM) {
                dbg("→ Found candidate hit");
                PlotOverlayWaveform(event, om_num, r_waveform);

                TString fname = Form("mem_cell_overlay_event%d_om%d.png",
                                    event, om_num);
                if (!gSystem->AccessPathName(fname)) {
                    dbg(("✓ PNG written: " + std::string(fname)).c_str());
                } else {
                    dbg(("⚠️  Tried to write " + std::string(fname) + " but it is NOT on disk").c_str());
                }
                overlayWritten = true;            // don’t try again
            }

            PlotOverlayWaveform(event, om_num, r_waveform);
        }
        if (!overlayWritten) {
            std::cout << "⚠️  No hit on OM " << TARGET_OM 
                    << " found in event " << TARGET_EVENT 
                    << " (or waveform empty)" << std::endl;
        }

        ///////////////////////////////////////////////////////////////////////////////////////
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
