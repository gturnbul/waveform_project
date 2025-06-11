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
    // Pass 1: Determine min/max ADC values for OM 152 bins
    std::map<int, short> min_adc_per_bin;
    std::map<int, short> max_adc_per_bin;

    for (int bin = 950; bin <= 1023; ++bin) {
        min_adc_per_bin[bin] = std::numeric_limits<short>::max();
        max_adc_per_bin[bin] = std::numeric_limits<short>::min();
    }

    // Overlay plot of all OM152 waveforms (all bins of each event)
    TCanvas* c_overlay = new TCanvas("c_overlay", "OM152 Overlay Waveforms", 1000, 600);

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("OM152 Waveforms Overlay;Bin Number;ADC Value");


    for (int event = 0; event < max_entries; ++event) {
        tree->GetEntry(event);

        for (int k = 0; k < calo_nohits; ++k) {
            int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);
            std::vector<short>& wave_k = wave->at(k);

            if (om_num == 152 && wave_k.size() > 1023) {
                for (int bin = 950; bin <= 1023; ++bin) {
                    short adc = wave_k[bin];
                    if (adc < min_adc_per_bin[bin]) min_adc_per_bin[bin] = adc;
                    if (adc > max_adc_per_bin[bin]) max_adc_per_bin[bin] = adc;
                }
                // Create TGraph for overlay
                int n_bins = wave_k.size();
                TGraph* gr = new TGraph(n_bins-950);
                for (int bin = 950; bin < n_bins; ++bin) {
                    gr->SetPoint(bin-950, bin, wave_k[bin]);
                }
                gr->SetLineColorAlpha(kBlue, 0.1);  // semi-transparent blue
                gr->SetLineWidth(1);
                mg->Add(gr);
            }
        }
    }
    // Now draw and save overlay plot
    mg->Draw("AL");  // Draw axes + all graphs
    c_overlay->Update();
    c_overlay->SaveAs("om152_overlay.png");

    // You can delete mg and c_overlay here if no longer needed
    delete mg;
    delete c_overlay;

    //////////////////////////////////////////////////////////////////////////
    // Allocate histograms using dynamic ranges
    std::map<int, TH1D*> bin_histograms;
    for (int bin = 950; bin <= 1023; ++bin) {
        short min_adc = min_adc_per_bin[bin] - 5;  // padding
        short max_adc = max_adc_per_bin[bin] + 5;

        TString hist_name = Form("om152_bin%d", bin);
        TString hist_title = Form("OM 152 - Bin %d;ADC Value;Counts", bin);
        bin_histograms[bin] = new TH1D(hist_name, hist_title, 100, min_adc, max_adc);
    }

    //////////////////////////////////////////////////////////////////////////
    // Pass 2: Fill baseline tree and histograms
    for (int event = 0; event < max_entries; ++event) {
        tree->GetEntry(event);
        for (int k = 0; k < calo_nohits; ++k) {
            int om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, k);
            std::vector<short>& wave_k = wave->at(k);
            int fcr_k = fcr->at(k);

            std::vector<short> r_waveform = reorder_waveform(wave_k, fcr_k);
            double baseline = calculate_baseline(wave_k);
            double stddev = calculate_stddev(wave_k, baseline);
            eom = calculate_eom(wave_k, baseline);

            if (om_num == 152 && wave_k.size() > 1023) {
                for (int bin = 950; bin <= 1023; ++bin) {
                    bin_histograms[bin]->Fill(wave_k[bin]);
                }
            }

            om_num_out = om_num;
            baseline_out = baseline;
            stddev_out = stddev;
            event_num = event;
            timestamp_out = timestamp->at(k) * 6.25e-9;
            timestamp_diff = calculate_timestamp_diff(om_num, timestamp_out);

            baseline_tree->Fill();  
        }
    }

    // Output the number of entries to the terminal
    cout << "Max entries: " << max_entries << "\n";

    //////////////////////////////////////////////////////////////////////////
    // Fit and draw each histogram
    TCanvas *c_fit = new TCanvas("c_fit", "Fits", 800, 600);
    for (int bin = 950; bin <= 1023; ++bin) {
        TH1D* hist = bin_histograms[bin];
        if (hist->GetEntries() > 10) {
            double xmin = hist->GetXaxis()->GetXmin();
            double xmax = hist->GetXaxis()->GetXmax();

            hist->Fit("gaus", "Q", "", xmin, xmax);
            hist->Draw();
            hist->GetFunction("gaus")->SetLineColor(kRed);
            hist->GetFunction("gaus")->SetLineWidth(2);

            TF1* fit_func = hist->GetFunction("gaus");
            if (fit_func) {
                double chi2 = fit_func->GetChisquare();
                int ndf = fit_func->GetNDF();
                double chi2_ndf = chi2 / ndf;
                double sigma = fit_func->GetParameter(2);
                double mean = fit_func->GetParameter(1);

                stddevs.push_back(sigma);
                bin_numbers.push_back(bin);
                red_chi2.push_back(chi2_ndf);
                g_mean.push_back(mean);

                TPaveText* pave = new TPaveText(0.6, 0.75, 0.88, 0.85, "NDC");
                pave->SetFillColor(0);
                pave->SetBorderSize(1);
                pave->AddText(Form("#chi^{2}/NDF = %.2f", chi2_ndf));
                pave->Draw();

                std::cout << "Bin " << bin
                          << ": chi2 = " << chi2
                          << ", NDF = " << ndf
                          << ", chi2/NDF = " << chi2_ndf << std::endl;
            }

            TString outname = Form("om152_bin%d_fit.png", bin);
            c_fit->SaveAs(outname);
        }

        delete hist;  // clean up memory
    }
    delete c_fit;

    // Plot standard deviation per bin
    TCanvas* c_stddev = new TCanvas("c_stddev", "StdDev per Bin", 800, 600);
    TGraph* gr_stddev = new TGraph(bin_numbers.size());

    for (size_t i = 0; i < bin_numbers.size(); ++i) {
        gr_stddev->SetPoint(i, bin_numbers[i], stddevs[i]);
    }

    gr_stddev->SetTitle("Standard Deviation from Gaussian Fit;Bin Number;#sigma [ADC]");
    gr_stddev->SetMarkerStyle(21);
    gr_stddev->SetLineColor(kRed);
    gr_stddev->SetLineWidth(2);
    gr_stddev->Draw("ALP");

    c_stddev->SaveAs("om152_stddev_per_bin.png");

    delete gr_stddev;
    delete c_stddev;

    // Plot reduced chi-squared per bin
    TCanvas* c_chi2 = new TCanvas("c_chi2", "Reduced Chi2 per Bin", 800, 600);
    TGraph* gr_chi2 = new TGraph(bin_numbers.size());

    for (size_t i = 0; i < bin_numbers.size(); ++i) {
        gr_chi2->SetPoint(i, bin_numbers[i], red_chi2[i]);
    }

    gr_chi2->SetTitle("Reduced #chi^{2} from Gaussian Fit;Bin Number;#chi^{2}/ndf");
    gr_chi2->SetMarkerStyle(22);
    gr_chi2->SetLineColor(kBlue + 2);
    gr_chi2->SetLineWidth(2);
    gr_chi2->Draw("ALP");

    c_chi2->SaveAs("om152_reduced_chi2_per_bin.png");

    // Clean up
    delete gr_chi2;
    delete c_chi2;


    // Plot reduced mean per bin
    TCanvas* c_mean = new TCanvas("c_mean", "Mean per Bin", 800, 600);
    TGraph* gr_mean = new TGraph(bin_numbers.size());

    for (size_t i = 0; i < bin_numbers.size(); ++i) {
        gr_mean->SetPoint(i, bin_numbers[i], g_mean[i]);
    }

    gr_mean->SetTitle("Mean from Gaussian Fit;Bin Number;Mean [ADC]");
    gr_mean->SetMarkerStyle(22);
    gr_mean->SetLineColor(kGreen + 2);
    gr_mean->SetLineWidth(2);
    gr_mean->Draw("ALP");

    c_mean->SaveAs("om152_mean_per_bin.png");

    delete gr_mean;
    delete c_mean;



    //////////////////////////////////////////////////////////////////////////
    // Write and close output
    baseline_tree->Write();
    outfile->Close();
    delete outfile;
    file->Close();
}

