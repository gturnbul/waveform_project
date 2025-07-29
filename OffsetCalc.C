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

////////////////////////////// FUNCTIONS from the WAVEFORMS ////////////////////////////////////
// Function to calculate baseline
double calculate_baseline976(const std::vector<short>& waveform) {
    double baseline = 0;
    for (int i = 0; i < 976; ++i) {
        baseline += waveform[i];
    }
    return baseline / 976;
}

double calculate_baseline48(const std::vector<short>& waveform) {
    double baseline = 0;
    for (int i = 976; i < 1024; ++i) {
        baseline += waveform[i];
    }
    return baseline / 48;
}

// Function to calculate chi2df
double calculate_stddev976(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 0; i < 976; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 976);
}

double calculate_stddev48(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0;
    for (int i = 976; i < 1024; ++i) {
        sum_sq_diff += pow(waveform[i] - baseline, 2);
    }
    return sqrt(sum_sq_diff / 48);
}
// Function to calculate error on the mean (EOM), with sigma calculated inside
double calculate_eom976(const std::vector<short>& waveform, double baseline) {
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
double calculate_eom48(const std::vector<short>& waveform, double baseline) {
    double sum_sq_diff = 0.0;

    for (int i = 976; i < 1024; ++i) {
        double diff = waveform[i] - baseline;
        sum_sq_diff += diff * diff;
    }

    double variance = sum_sq_diff / 48;          // Population variance
    double sigma = std::sqrt(variance);          // Standard deviation
    double eom = sigma / std::sqrt(48);          // Error on the mean

    return eom;
}

struct Mod16Stats {
    std::vector<double> offsets;
    std::vector<double> sigmas;
    std::vector<double> chi2ndfs;
};

Mod16Stats compute_mod16_from_fits(const std::vector<double>& means,
                                    const std::vector<double>& sigmas,
                                    const std::vector<double>& chi2ndfs) {
    std::vector<double> mean_acc(16, 0.0);
    std::vector<double> sigma_acc(16, 0.0);
    std::vector<double> chi2ndf_acc(16, 0.0);
    std::vector<int> counts(16, 0);

    for (int i = 0; i < 976; ++i) {
        int mod16_index = i % 16;
        mean_acc[mod16_index] += means[i];
        sigma_acc[mod16_index] += sigmas[i];
        chi2ndf_acc[mod16_index] += chi2ndfs[i];
        counts[mod16_index]++;
    }

    for (int i = 0; i < 16; ++i) {
        if (counts[i] > 0) {
            mean_acc[i] /= counts[i];
            sigma_acc[i] /= counts[i];
            chi2ndf_acc[i] /= counts[i];
        }
    }

    return { mean_acc, sigma_acc, chi2ndf_acc };
}


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
std::vector<short> inverse_reorder_waveform(const std::vector<short>& reordered, int fcr) {
    int max_bins = 976;  // original waveform length
    
    std::vector<short> original(max_bins, EMPTY_BIN);
    
    for (int i = 0; i < max_bins; ++i) {
        int reordered_index = get_cell(i, fcr);
        // Safety check: reordered_index should be within bounds
        if (reordered_index >= 0 && reordered_index < (int)reordered.size()) {
            original[i] = reordered[reordered_index];
        } else {
            // handle unexpected index (optional)
            original[i] = EMPTY_BIN;
        }
    }
    
    return original;
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
    TFile *outfile = new TFile("baseline_offset_calc_1143.root", "RECREATE");
    TTree *baseline_tree = new TTree("baseline_tree", "OM baseline data");

    // Initialise output variables
    int event_num = -1;
    int om_num = -1;
    double baseline_orig = 0;
    double stddev_orig = 0;
    double eom_orig = 0;

    // Variables for the 1024 FEB offset correction
    double baseline_feb = 0;
    double stddev_feb = 0;
    double eom_feb = 0;
  
    //first 976 variables after wavecatcher (mod16) correction
    double baseline_976 = 0;
    double stddev_976 = 0;
    double eom_976 = 0;

    //end spike variables after end of waveform correction
    double baseline_end = 0;
    double stddev_end = 0;
    double eom_end = 0;

    //Quality assurance variables
    std::vector<double> mean_mem;
    std::vector<double> sigma_mem;
    std::vector<double> chi2ndf_mem; 
    
    std::vector<double> mean_timeo;
    std::vector<double> sigma_timeo;
    std::vector<double> chi2ndf_timeo;

    // Create branches in the output tree
    baseline_tree->Branch("event_num", &event_num, "event_num/I");
    baseline_tree->Branch("om_num", &om_num, "om_num/I");
    baseline_tree->Branch("baseline", &baseline_orig, "baseline/D");
    baseline_tree->Branch("stddev", &stddev_orig, "stddev/D");
    baseline_tree->Branch("eom", &eom_orig, "eom/D");

    baseline_tree->Branch("baseline_feb", &baseline_feb, "baseline_feb/D");
    baseline_tree->Branch("stddev_feb", &stddev_feb, "stddev_feb/D");
    baseline_tree->Branch("eom_feb", &eom_feb, "eom_feb/D");

    baseline_tree->Branch("baseline_976", &baseline_976, "baseline_976/D");
    baseline_tree->Branch("stddev_976", &stddev_976, "stddev_976/D");
    baseline_tree->Branch("eom_976", &eom_976, "eom_976/D");

    baseline_tree->Branch("baseline_end", &baseline_end, "baseline_end/D");
    baseline_tree->Branch("stddev_end", &stddev_end, "stddev_end/D");
    baseline_tree->Branch("eom_end", &eom_end, "eom_end/D");

    baseline_tree->Branch("mean_mem", &mean_mem);
    baseline_tree->Branch("sigma_mem", &sigma_mem);
    baseline_tree->Branch("chi2ndf_mem", &chi2ndf_mem);    
    
    baseline_tree->Branch("mean_timeo", &mean_timeo);
    baseline_tree->Branch("sigma_timeo", &sigma_timeo);
    baseline_tree->Branch("chi2ndf_timeo", &chi2ndf_timeo);
    //////////////////////////// MAPS ////////////////////////////////////////
    std::map<int, std::vector<std::vector<float>>> adc_values;
    std::map<int, std::vector<double>> mod16_offsets_per_om;
    std::map<int, std::vector<double>> mod16_sigmas_per_om;
    std::map<int, std::vector<double>> mod16_chi2ndf_per_om;

    std::map<int, std::vector<std::vector<float>>> feb_adc_values;
    std::map<int, std::vector<double>> feb_offsets_per_om;
    std::map<int, std::vector<double>> feb_sigmas_per_om;
    std::map<int, std::vector<double>> feb_chi2ndf_per_om;

    std::map<int, std::vector<double>> eow_offsets_per_om;
    std::map<int, std::vector<double>> eow_sigmas_per_om;
    std::map<int, std::vector<double>> eow_chi2ndf_per_om;

    //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////// Finding and saving ADC_VAL range limit per OM //////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < max_entries; ++i){
            tree->GetEntry(i);
            event_num = i;
            if (i < 10) continue; // Skip first 10 events due to cleanliness of early data

            for (int j = 0; j <calo_nohits; ++j){
                if (wave->at(j).size() < 1024) continue; // Skip if waveform is too short

                //calculate the om number
                int om_num_out = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, j);

                // Calculate original baseline, stddev, and eom
                baseline_orig = calculate_baseline976(wave->at(j));
                stddev_orig = calculate_stddev976(wave->at(j), baseline_orig);
                eom_orig = calculate_eom976(wave->at(j), baseline_orig);

                // only include waveforms where the eom is less than 0.076
                if (eom_orig > 0.076) continue;
                
                //////////////////////////// Wavecatcher data /////////////////////////////////////
                // Initialise storage for the adc counts of bins per optical module
                if (adc_values.find(om_num_out) == adc_values.end()) {
                    adc_values[om_num_out] = std::vector<std::vector<float>>(1024, std::vector<float>());
                }
                
                // loop each bin in the waveform and collect the ADC values per bin
                for (int bin = 0; bin < 1024; ++bin) {
                    adc_values[om_num_out][bin].push_back(wave->at(j)[bin]);
                }

                ///////////////////////////// FEB data /////////////////////////////////////////////
                // initialise storage for the ADC counts of bins per om for FEB
                if (feb_adc_values.find(om_num_out) == feb_adc_values.end()) {
                    feb_adc_values[om_num_out] = std::vector<std::vector<float>>(1024, std::vector<float>());

                }
                std::vector<short> reordered_waveform = reorder_waveform(wave->at(j), fcr->at(j));
                for (int bin = 0; bin < 1024; ++bin) {
                    feb_adc_values[om_num_out][bin].push_back(reordered_waveform[bin]);
                }

            }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// FEB Calclulations /////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    std::ofstream csv_file_feb("feb_offsets.csv");
    csv_file_feb << "om,bin,mean,sigma,chi2ndf\n";

    // Process FEB_ADC_VALUES for histogramming
    for (auto& om_pair : feb_adc_values) {
        int om_num = om_pair.first;

        // if (om_num != 1) continue; // Skip all OMs except OM 1 for now -----------------
        // TCanvas *c = new TCanvas("c", "OM 1 Histograms", 800, 600); //------------------

        auto& adc_bins = om_pair.second;

        std::vector<double> fit_means_feb_om; // filled with 976 Gaussian means
        std::vector<double> fit_sigmas_feb_om; // filled with 976 Gaussian sigmas
        std::vector<double> fit_chi2_ndf_feb_om; // filled with 976 Gaussian chi2/ndf values

        for (int bin = 0; bin < 1024; ++bin){
            std::vector<float>& adc_values_bin = adc_bins[bin];

            // Filter out EMPTY_BIN values
            std::vector<float> valid_adc_values;
            for (float val : adc_values_bin) {
                if (static_cast<short>(val) != EMPTY_BIN) {
                    valid_adc_values.push_back(val);
                }
            }

            if (valid_adc_values.empty()) {
                // No valid ADC values in this bin; skip histogramming/fitting
                continue;
            }

            // Calculate histogram range from valid ADC values only
            float min_adc = *std::min_element(valid_adc_values.begin(), valid_adc_values.end()) - 5;
            float max_adc = *std::max_element(valid_adc_values.begin(), valid_adc_values.end()) + 5;
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));

            TString hist_name = Form("feb_om%d_bin%d", om_num, bin);
            TString hist_title = Form("FEB OM %d - Bin %d;ADC Value;Counts", om_num, bin);

            TH1D *hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);

            for (float adc_value : valid_adc_values) {
                hist->Fill(adc_value);
            }

            TF1 *fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);

            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = (ndf != 0) ? chi2 / ndf : 0;

            fit_means_feb_om.push_back(mean);
            fit_sigmas_feb_om.push_back(sigma);
            fit_chi2_ndf_feb_om.push_back(chi2_ndf);

            csv_file_feb << om_num << ","
                        << bin << ","
                        << mean << ","
                        << sigma << ","
                        << chi2_ndf << "\n";

            hist->Draw();
            // TString filename = Form("feb_bin%d_histogram.png", bin);
            // c->SaveAs(filename);

            delete hist;
            delete fit_func;
        }
        // Store the results for this OM
        feb_offsets_per_om[om_num] = fit_means_feb_om;
        feb_sigmas_per_om[om_num] = fit_sigmas_feb_om;
        feb_chi2ndf_per_om[om_num] = fit_chi2_ndf_feb_om;

        // delete c; // Delete the canvas after use
 
    }
    csv_file_feb.close();

    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// Wavecater Calclulations ///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    // Process ADC_VALUES for histogramming, fitting and mod16 calculations
    for (auto& om_pair : adc_values) {
        int om_num = om_pair.first;
        auto& adc_bins = om_pair.second;

        std::vector<double> fit_means_om; // filled with 976 Gaussian means
        std::vector<double> fit_sigmas_om; // filled with 976 Gaussian sigmas
        std::vector<double> fit_chi2_ndf_om; // filled with 976 Gaussian chi2/ndf values

        for (int bin = 0; bin < 976; ++bin){
            std::vector<float>& adc_values_bin = adc_bins[bin];

            // Create and fill histogram
            TString hist_name = Form("om%d_bin%d", om_num, bin);
            TString hist_title = Form("OM %d - Bin %d;ADC Value;Counts", om_num, bin);
            // set histogram range based on the ADC values per om
            float min_adc = *std::min_element(adc_values_bin.begin(), adc_values_bin.end())-5;
            float max_adc = *std::max_element(adc_values_bin.begin(), adc_values_bin.end())+5;
            // Ensure bin count is an integer:
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));
            // create histogram with nbins between min and max ADC values
            TH1D *hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);
            for (float adc_value : adc_values_bin) {
                hist->Fill(adc_value);
            }

            // Fit gaussian, extract mean/sigma
            TF1 *fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);
            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = chi2 / ndf;
            fit_means_om.push_back(mean);
            fit_sigmas_om.push_back(sigma);
            fit_chi2_ndf_om.push_back(chi2_ndf);

            delete hist;
            delete fit_func;

        }

        // Calculate mod16 offsets for this OM
        Mod16Stats stats = compute_mod16_from_fits(fit_means_om, fit_sigmas_om, fit_chi2_ndf_om);
        mod16_offsets_per_om[om_num] = stats.offsets;
        // optionally store sigmas and chi2ndf too:
        mod16_sigmas_per_om[om_num] = stats.sigmas;
        mod16_chi2ndf_per_om[om_num] = stats.chi2ndfs;
    }
    // Open CSV file for writing
    std::ofstream csv_file("mod16_offsets.csv");

    // Write header
    csv_file << "om,offset_index,offset_value,sigma_value,chi2ndf_value\n";

    // Loop through the map and write data
    for (const auto& pair : mod16_offsets_per_om) {
        int om = pair.first;
        const std::vector<double>& offsets = pair.second;
        for (int i = 0; i < (int)offsets.size(); ++i) {
            csv_file << om << ","
                 << i << ","
                 << mod16_offsets_per_om[om][i] << ","
                 << mod16_sigmas_per_om[om][i] << ","
                 << mod16_chi2ndf_per_om[om][i] << "\n";
        }
    }

    csv_file.close();
    // fit_file->Close();

    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  End spike offsets  /////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    std::ofstream spike_csv_file("eow_spike_offsets.csv");
    spike_csv_file << "om,bin,offset_value,sigma_value,chi2ndf_value\n";

    for (auto& om_pair : adc_values) {
    int om_num = om_pair.first;
    auto& adc_bins = om_pair.second;

        // vectors to store spike offsets for this OM
        std::vector<double> spike_means_om;
        std::vector<double> spike_sigmas_om;
        std::vector<double> spike_chi2ndf_om;

        auto it_mod16 = mod16_offsets_per_om.find(om_num);
        const std::vector<double>& mod16_offsets = (it_mod16 != mod16_offsets_per_om.end()) ? it_mod16->second : std::vector<double>{};


        // Loop bins 976-1023 
        for (int bin = 976; bin <= 1023; ++bin) {
            std::vector<float>& adc_values_bin = adc_bins[bin];

            TString hist_name = Form("om%d_bin%d", om_num, bin);
            TString hist_title = Form("OM %d - Bin %d;ADC Value;Counts", om_num, bin);
            float min_adc = *std::min_element(adc_values_bin.begin(), adc_values_bin.end()) - 5;
            float max_adc = *std::max_element(adc_values_bin.begin(), adc_values_bin.end()) + 5;
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));
            TH1D* hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);
            for (float adc_value : adc_values_bin) {
                hist->Fill(adc_value);
            }

            TF1* fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);
            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = (ndf > 0) ? chi2 / ndf : -1;

            // Calculate difference: spike mean - mod16 offset for (bin % 16)
            int mod16_index = bin % 16;
            double mod16_value = (mod16_index < (int)mod16_offsets.size()) ? mod16_offsets[mod16_index] : 0.0;
            double corrected_spike_offset = mean - mod16_value;

            // Save results directly instead of computing mod16 offsets
            spike_means_om.push_back(corrected_spike_offset);
            spike_sigmas_om.push_back(sigma);
            spike_chi2ndf_om.push_back(chi2_ndf);

            delete hist;
            delete fit_func;
        }
        // Store the results for this OM
        eow_offsets_per_om[om_num] = spike_means_om;
        eow_sigmas_per_om[om_num] = spike_sigmas_om;
        eow_chi2ndf_per_om[om_num] = spike_chi2ndf_om;


    }

    // Write to CSV

    for (const auto& om_pair : eow_offsets_per_om) {
        int om_num = om_pair.first;
        const auto& offsets = om_pair.second;
        const auto& sigmas = eow_sigmas_per_om[om_num];
        const auto& chi2s = eow_chi2ndf_per_om[om_num];

        for (int i = 0; i < (int)offsets.size(); ++i) {
            int bin = 976 + i;
            spike_csv_file << om_num << ","
                        << bin << ","
                        << offsets[i] << ","
                        << sigmas[i] << ","
                        << chi2s[i] << "\n";
        }
    }
    spike_csv_file.close();



    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Baseline Tree Fill /////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////

    
    for (int i = 0; i < max_entries; ++i){
        tree->GetEntry(i);
        event_num = i;
        if (i < 10) continue; // Skip first 10 events due to cleanliness of early data

        for (int j = 0; j <calo_nohits; ++j){
            if (wave->at(j).size() < 1024) continue; // Skip if waveform is too short

            //calculate the om number
            om_num = calculate_om_num(calo_type, calo_side, calo_wall, calo_column, calo_row, j);

            // Calculate original baseline, stddev, and eom
            baseline_orig = calculate_baseline976(wave->at(j));
            stddev_orig = calculate_stddev976(wave->at(j), baseline_orig);
            eom_orig = calculate_eom976(wave->at(j), baseline_orig);

            // only include waveforms where the eom is less than 0.076
            if (eom_orig > 0.076) continue;

            // feb waveform correction //////////////////////////////////////////////////////////////
            // Step 1: reorder the original waveform
            std::vector<short> reordered_wave = reorder_waveform(wave->at(j), fcr->at(j));

            // Step 2: find the FEB offsets for this OM and subtract them from reordered waveform
            auto it_feb = feb_offsets_per_om.find(om_num);
            if (it_feb != feb_offsets_per_om.end()) {
                const std::vector<double>& feb_offsets = it_feb->second;
                for (int bin = 0; bin < 1024; ++bin) {
                    reordered_wave[bin] -= feb_offsets[bin];
                }
            } else {
                std::cerr << "Warning: No FEB offsets found for OM " << om_num << std::endl;
            }

            // Step 3: invert reorder to get corrected waveform back in original order
            std::vector<short> corrected_wave = inverse_reorder_waveform(reordered_wave, fcr->at(j));

            // Step 4: calculate baseline, stddev, and eom on corrected waveform in original order
            baseline_feb = calculate_baseline976(corrected_wave);
            stddev_feb   = calculate_stddev976(corrected_wave, baseline_feb);
            eom_feb      = calculate_eom976(corrected_wave, baseline_feb);

            // 976 waveform correction ///////////////////////////////////////////////////////////////
            auto it = mod16_offsets_per_om.find(om_num);

            if (it != mod16_offsets_per_om.end()) {
                const std::vector<double>& offsets = it->second;

                for (int bin = 0; bin < 976; ++bin) {
                    wave->at(j)[bin] -= offsets[bin % 16]; // Apply mod16 offset
                }
            } else {
                // OM not found in offsets map, handle error or skip correction
                std::cerr << "Warning: No mod16 offsets found for OM " << om_num << std::endl;
            }

            // Calculate adjusted baseline, stddev, and eom
            baseline_976 = calculate_baseline976(wave->at(j));
            stddev_976 = calculate_stddev976(wave->at(j), baseline_976);
            eom_976 = calculate_eom976(wave->at(j), baseline_976);
         
            // Spike correction for bins 976–1023 /////////////////////////////////////////////////
            auto it_spike = eow_offsets_per_om.find(om_num);

            if (it_spike != eow_offsets_per_om.end()) {
                const std::vector<double>& spike_offsets = it_spike->second;

                for (int bin = 976; bin <= 1023; ++bin) {
                    int spike_bin_index = bin - 976;
                    if (spike_bin_index < (int)spike_offsets.size()) {
                        wave->at(j)[bin] -= spike_offsets[spike_bin_index]; //remove eow spike to mod16 offset

                        // Apply mod16 offset to spike-corrected bin
                        auto it_mod16 = mod16_offsets_per_om.find(om_num);  
                        if (it_mod16 != mod16_offsets_per_om.end()) {
                            const std::vector<double>& offsets = it_mod16->second;
                            wave->at(j)[bin] -= offsets[bin % 16];         // Apply mod16 offset
                        }
                    }
                }
            } else {
                std::cerr << "Warning: No spike offsets found for OM " << om_num << std::endl;
            }

            // Recalculate baseline, stddev, and eom over the full waveform (if needed)
            baseline_end = calculate_baseline48(wave->at(j));
            stddev_end   = calculate_stddev48(wave->at(j), baseline_end);
            eom_end      = calculate_eom48(wave->at(j), baseline_end);

            // storing quality of fit assurance data ////////////////////////////////////////////////////
            // Clear previous event's data
            mean_timeo.clear();
            sigma_timeo.clear();
            chi2ndf_timeo.clear();

            auto it_offset = mod16_offsets_per_om.find(om_num);
            auto it_sigma  = mod16_sigmas_per_om.find(om_num);
            auto it_chi2   = mod16_chi2ndf_per_om.find(om_num);

            if (it_offset != mod16_offsets_per_om.end() &&
                it_sigma != mod16_sigmas_per_om.end() &&
                it_chi2 != mod16_chi2ndf_per_om.end()) {

                mean_timeo = it_offset->second;
                sigma_timeo = it_sigma->second;
                chi2ndf_timeo = it_chi2->second;

            } else {
                std::cerr << "Warning: Missing QA data for OM " << om_num << std::endl;

                // Optionally fill with NaNs or -999 if QA data missing
                mean_timeo.assign(16, -999);
                sigma_timeo.assign(16, -999);
                chi2ndf_timeo.assign(16, -999);
            }

            // feb QA ///////////////////////////////////////////////////////////////
            mean_mem.clear();
            sigma_mem.clear();
            chi2ndf_mem.clear();

            auto it_feb_offset = feb_offsets_per_om.find(om_num);
            auto it_feb_sigma  = feb_sigmas_per_om.find(om_num);
            auto it_feb_chi2   = feb_chi2ndf_per_om.find(om_num);

            if (it_feb_offset != feb_offsets_per_om.end() &&
                it_feb_sigma != feb_sigmas_per_om.end() &&
                it_feb_chi2 != feb_chi2ndf_per_om.end()) {

                mean_mem = it_feb_offset->second;
                sigma_mem = it_feb_sigma->second;
                chi2ndf_mem = it_feb_chi2->second;

            } else {
                std::cerr << "Warning: Missing FEB QA data for OM " << om_num << std::endl;

                // Optionally fill with NaNs or -999 if QA data missing
                mean_mem.assign(1024, -999);
                sigma_mem.assign(1024, -999);
                chi2ndf_mem.assign(1024, -999);

            }
            
            // Fill the tree 
            baseline_tree->Fill();
        } 
    }
    

    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// End of code, writes and closes ////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // write and close
    outfile->cd();
    baseline_tree->Write();
    outfile->Close();
    file->Close();

    return 0;
}
