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
    for (int i = 976; i < 1024; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 48);
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

////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> compute_mod16(const std::vector<short>& waveform) {
    std::vector<double> averages(16, 0.0);
    std::vector<int> counts(16, 0);

    for (int i = 0; i < 976; ++i) {
        int batch_index = i % 16;
        averages[batch_index] += waveform[i];
        counts[batch_index]++;
    }

    for (int i = 0; i < 16; ++i) {
        if (counts[i] > 0)
            averages[i] /= counts[i];
    }

    return averages;
}

//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Main ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
int main(){
std::ofstream csv_file("stddev_output.csv");
csv_file << "om_num,bin_num,original_stddev,processed_stddev\n";  // CSV header

TFile* stddev_outfile = new TFile("baseline_stddev_1143.root", "RECREATE");
TTree* stddev_tree = new TTree("stddev_tree", "OM standard deviation data");

int om_num_stddev;
int bin_number;  // This will hold the bin number, can be adjusted to a range if needed
double original_stddev;
double processed_stddev;

stddev_tree->Branch("om_num", &om_num_stddev, "om_num/I");
stddev_tree->Branch("bin_number", &bin_number, "bin_number/I");  // Branch for bin number
stddev_tree->Branch("original_stddev", &original_stddev, "original_stddev/D");
stddev_tree->Branch("processed_stddev", &processed_stddev, "processed_stddev/D");

TH1D* h_original_stddev = new TH1D("h_original_stddev", "Original StdDev;StdDev;Entries", 100, 0, 100);
TH1D* h_processed_stddev = new TH1D("h_processed_stddev", "Processed StdDev;StdDev;Entries", 100, 0, 100);



    //oms to analyse
    std::vector<int> om_list;          // OMs to analyse
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
    for (int v = 0; v <= 711; ++v){
        om_list.push_back(v);
    }


    for (int target_om : om_list){
        std::map<int, short> min_adc_per_bin, max_adc_per_bin;
        std::vector<double> means, sigmas;
        std::vector<int> bins;

        for (int bin = 976; bin <= 1023; ++bin) {
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

                for (int bin = 976; bin <= 1023; ++bin) {
                    short adc = wave_k[bin];
                    if (adc < min_adc_per_bin[bin]) min_adc_per_bin[bin] = adc;
                    if (adc > max_adc_per_bin[bin]) max_adc_per_bin[bin] = adc;
                }
            }
        }

        // Create histograms
        std::map<int, TH1D*> bin_histograms;
        for (int bin = 976; bin <= 1023; ++bin) {
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

                for (int bin = 976; bin <= 1023; ++bin)
                    bin_histograms[bin]->Fill(wave_k[bin]);
            }
        }

        // Gaussian fits
        std::unordered_map<int, double> bin_means;
        means.clear();
        sigmas.clear();
        bins.clear();

        for (int bin = 976; bin <= 1023; ++bin) {
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

        const int BIN_START = 976;
        const int BIN_END = 1023;
        const int N_BINS = BIN_END - BIN_START + 1;

        for (int event = 10; event < max_entries; ++event) {
            tree->GetEntry(event);

            for (int k = 0; k < calo_nohits; ++k) {
                int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);
                if (om_num != target_om) continue;

                const std::vector<short>& wave_k = wave->at(k);
                if (wave_k.size() <= BIN_END) continue;

                std::vector<short> mod16_wave_k = wave_k; // Copy original waveform
                std::vector<double> avg_mod16 = compute_mod16(wave_k);

                // Replace bins 976–1023 with the computed averages
                for (int i = 976; i <= 1023; ++i) {
                    int avg_index = (i - 976) % 16;
                    mod16_wave_k[i] = static_cast<short>(avg_mod16[avg_index]);
                }


                // Compute baseline and stddev of full waveform
                double baseline_orig = calculate_baseline(wave_k);
                double stddev_orig = calculate_stddev(wave_k, baseline_orig);

                // Subtract bin mean for tail region
                std::vector<double> subtracted(N_BINS);
                double mean_sub = 0.0;

                for (int i = 0; i < N_BINS; ++i) {
                    int bin = BIN_START + i;
                    double mean = bin_means.count(bin) ? bin_means[bin] : 0.0;
                    subtracted[i] = wave_k[bin] - mean;
                    mean_sub += subtracted[i];
                }
                mean_sub /= N_BINS;

                double sum_sq_diff_sub = 0.0;
                for (int i = 0; i < N_BINS; ++i) {
                    sum_sq_diff_sub += std::pow(subtracted[i] - mean_sub, 2);
                }

                double stddev_sub = std::sqrt(sum_sq_diff_sub / N_BINS);

                // Write results
                csv_file << om_num << "," << BIN_START << "-" << BIN_END << "," 
                        << stddev_orig << "," << stddev_sub << "\n";

                h_original_stddev->Fill(stddev_orig);
                h_processed_stddev->Fill(stddev_sub);

                // Fill tree
                om_num_stddev = om_num;
                bin_number = BIN_START; // If you want to use a range, adjust variable type or format
                original_stddev = stddev_orig;
                processed_stddev = stddev_sub;
                stddev_tree->Fill();

        


                // for (event = 10; event <15; ++event){
                //     // -----  ORIGINAL waveform plot  -----
                //     TCanvas* c_orig = new TCanvas(
                //         Form("c_orig_om%d_evt%d", target_om, event),
                //         Form("OM %d Event %d - Original", target_om, event), 800, 600);

                //     TGraph* g_orig = new TGraph(1024);  // all bins 0 to 1023
                //     for (int i = 0; i < 1024; ++i)
                //         g_orig->SetPoint(i, i, wave_k[i]);

                //     g_orig->SetTitle(
                //         Form("OM %d Event %d - Original;Bin;ADC", target_om, event));
                //     g_orig->SetLineColor(kBlack);
                //     g_orig->Draw("AL");
                //     c_orig->SaveAs(Form("original_om%d_evt%d.png", target_om, event));
                //     delete c_orig;
                //     delete g_orig;


                //     TCanvas* c_sub = new TCanvas(Form("c_sub_om%d_evt%d", target_om, event),
                //                                 Form("OM %d Event %d - Mean-subtracted", target_om, event), 800, 600);
                //     TGraph* g_sub = new TGraph(subtracted.size());
                //     for (size_t i = 0; i < subtracted.size(); ++i)
                //         g_sub->SetPoint(i, 0 + static_cast<int>(i), subtracted[i]);

                //     g_sub->SetTitle(Form("OM %d Event %d - Mean-subtracted;Bin;ADC-Mean", target_om, event));
                //     g_sub->SetLineColor(kBlue);
                //     g_sub->Draw("AL");
                //     c_sub->SaveAs(Form("subtracted_om%d_evt%d.png", target_om, event));
                //     delete c_sub;
                //     delete g_sub;
                // }

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
    stddev_outfile->cd();
    stddev_tree->Write();
    h_original_stddev->Write();
    h_processed_stddev->Write();
    stddev_outfile->Close();
    delete stddev_outfile;
    delete h_original_stddev;
    delete h_processed_stddev;
    csv_file.close();
    file->Close();
}
