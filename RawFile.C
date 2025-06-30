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
#include <unordered_map>

// Macro Definition
#define _USE_MATH_DEFINES

using namespace std;

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
std::ofstream csv_file("stddev_output.csv");
csv_file << "event,om_num,original_stddev,processed_stddev\n";  // CSV header

    //oms to analyse
    std::vector<int> om_list = {37,39,47,66,107,112,131,168,170,173,181,184,214,217,229,250,258,262,266,269,270,292,311,335,363,386,389,391,402,408,446,455,458,476,482,485,511,523,524,527,528,531,540,542,542,555,559,573,580,596,614,617,625,633,638,640,672,685,688,695};          // OMs to analyse
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

    // *** containers keyed by OM ***
    std::unordered_map<int,TCanvas*>     c_overlay;
    std::unordered_map<int,TMultiGraph*> mg_overlay;
    std::map<std::pair<int,int>,short>   adc_min, adc_max; // (om,bin)->min/max
    std::map<std::pair<int,int>,TH1D*>   h_bin;            // (om,bin)->hist

    // Create output ROOT file and TTree for baseline data
    TFile *outfile = new TFile("baseline_output_tree1143.root", "RECREATE");
    TTree *baseline_tree = new TTree("baseline_tree", "OM baseline data");

    // Initialise output variables
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

    // vectors for printing end of waveform spike information
    std::vector<double> stddevs;
    std::vector<int> bin_numbers;
    std::vector<double> red_chi2;
    std::vector<double> g_mean;


    //////////////////////////////////////////////////////////////////////////
   #include <unordered_map>  // add at the top with other headers if not already included

for (int target_om : om_list)
{
    std::map<int, short> min_adc_per_bin, max_adc_per_bin;
    std::vector<double> means, sigmas;
    std::vector<int> bins;

    for (int bin = 0; bin <= 1023; ++bin) {
        min_adc_per_bin[bin] = std::numeric_limits<short>::max();
        max_adc_per_bin[bin] = std::numeric_limits<short>::min();
    }

    // Pass-1: Get min/max ADCs
    for (int event = 10; event < max_entries; ++event) {
        tree->GetEntry(event);
        for (int k = 0; k < calo_nohits; ++k) {
            int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);
            if (om_num != target_om) continue;

            const std::vector<short>& wave_k = wave->at(k);
            if (wave_k.size() <= 1023) continue;

            for (int bin = 0; bin <= 1023; ++bin) {
                short adc = wave_k[bin];
                if (adc < min_adc_per_bin[bin]) min_adc_per_bin[bin] = adc;
                if (adc > max_adc_per_bin[bin]) max_adc_per_bin[bin] = adc;
            }
        }
    }

    // Create histograms
    std::map<int, TH1D*> bin_histograms;
    for (int bin = 0; bin <= 1023; ++bin) {
        int lo = min_adc_per_bin[bin] - 5;
        int hi = max_adc_per_bin[bin] + 5;
        int nb = hi - lo + 1;
        TString hname = Form("om%d_bin%d", target_om, bin);
        TString htitle = Form("OM %d  Bin %d;ADC;Counts", target_om, bin);
        bin_histograms[bin] = new TH1D(hname, htitle, nb, lo, hi);
    }

    // Pass-2: Fill histograms
    for (int event = 10; event < max_entries; ++event) {
        tree->GetEntry(event);
        for (int k = 0; k < calo_nohits; ++k) {
            int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);
            if (om_num != target_om) continue;

            const std::vector<short>& wave_k = wave->at(k);
            if (wave_k.size() <= 1023) continue;

            for (int bin = 0; bin <= 1023; ++bin)
                bin_histograms[bin]->Fill(wave_k[bin]);
        }
    }

    // Gaussian fits
    std::unordered_map<int, double> bin_means;
    means.clear();
    sigmas.clear();
    bins.clear();

    for (int bin = 0; bin <= 1023; ++bin) {
        TH1D* h = bin_histograms[bin];
        if (h->GetEntries() < 10) continue;
        h->Fit("gaus", "Q0");
        TF1* f = h->GetFunction("gaus");
        if (!f) continue;

        double mean = f->GetParameter(1);
        double sigma = f->GetParameter(2);

        bin_means[bin] = mean;
        means.push_back(mean);
        sigmas.push_back(sigma);
        bins.push_back(bin);
    }

    // OPTIONAL: plot mean per bin
    TGraphErrors* g = new TGraphErrors(bins.size());
    for (size_t i = 0; i < bins.size(); ++i) {
        g->SetPoint(i, bins[i], means[i]);
        g->SetPointError(i, 0, sigmas[i]);
    }
    TCanvas* c_mean = new TCanvas(Form("c_mean_om%d", target_om),
                                  Form("OM%d mean per bin", target_om), 800, 600);
    g->SetTitle(Form("OM %d - Mean from Gaussian Fit;Bin;Mean [ADC]", target_om));
    g->SetLineColor(kRed);
    g->Draw("ALP");
    //c_mean->SaveAs(Form("om%d_mean_per_bin.png", target_om));
    delete c_mean;

    // Event-by-event subtracted plots
    for (int event = 10; event < 20; ++event) {
        tree->GetEntry(event);
        for (int k = 0; k < calo_nohits; ++k) {
            int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);
            if (om_num != target_om) continue;

            const std::vector<short>& wave_k = wave->at(k);
            if (wave_k.size() <= 1023) continue;

            std::vector<double> subtracted;
            for (int bin = 0; bin <= 1023; ++bin) {
                double mean = bin_means.count(bin) ? bin_means[bin] : 0.0;
                subtracted.push_back(wave_k[bin] - mean);
            }
            // Add the new stddev calculations here:
            double baseline_orig = calculate_baseline(wave_k);
            double stddev_orig = calculate_stddev(wave_k, baseline_orig);

            double mean_sub = 0;
            for (auto val : subtracted) mean_sub += val;
            mean_sub /= subtracted.size();

            double sum_sq_diff_sub = 0;
            for (auto val : subtracted) sum_sq_diff_sub += (val - mean_sub) * (val - mean_sub);
            double stddev_sub = std::sqrt(sum_sq_diff_sub / subtracted.size());

            csv_file << event << "," << om_num << "," << stddev_orig << "," << stddev_sub << "\n";

  
            // -----  ORIGINAL waveform plot  -----
            TCanvas* c_orig = new TCanvas(
                Form("c_orig_om%d_evt%d", target_om, event),
                Form("OM %d Event %d - Original", target_om, event), 800, 600);

            TGraph* g_orig = new TGraph(1024);  // all bins 0 to 1023
            for (int i = 0; i < 1024; ++i)
                g_orig->SetPoint(i, i, wave_k[i]);

            g_orig->SetTitle(
                Form("OM %d Event %d - Original;Bin;ADC", target_om, event));
            g_orig->SetLineColor(kBlack);
            g_orig->Draw("AL");
            c_orig->SaveAs(Form("original_om%d_evt%d.png", target_om, event));
            delete c_orig;
            delete g_orig;


            TCanvas* c_sub = new TCanvas(Form("c_sub_om%d_evt%d", target_om, event),
                                         Form("OM %d Event %d - Mean-subtracted", target_om, event), 800, 600);
            TGraph* g_sub = new TGraph(subtracted.size());
            for (size_t i = 0; i < subtracted.size(); ++i)
                g_sub->SetPoint(i, 0 + static_cast<int>(i), subtracted[i]);

            g_sub->SetTitle(Form("OM %d Event %d - Mean-subtracted;Bin;ADC-Mean", target_om, event));
            g_sub->SetLineColor(kBlue);
            g_sub->Draw("AL");
            c_sub->SaveAs(Form("subtracted_om%d_evt%d.png", target_om, event));
            delete c_sub;
            delete g_sub;

        }
    }

    // Clean up histograms
    for (auto& kv : bin_histograms) delete kv.second;
}

    //////////////////////////////////////////////////////////////////////////
    // Write and close output
    baseline_tree->Write();
    outfile->Close();
    delete outfile;
    file->Close();
}
