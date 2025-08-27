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
#include <TVirtualFFT.h>
#include <complex>
#include <TComplex.h>
#include <map>
#include <vector>
#include <cmath>
#include <THStack.h>

#include "TSpectrum.h"

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
        if (chi2ndfs[i] >= 0) {
            chi2ndf_acc[mod16_index] += chi2ndfs[i];
            counts[mod16_index]++;
            } else {
                static std::ofstream badfit_log("bad_fits_mod16.csv", std::ios::app);
                badfit_log << "OM_UNKNOWN," << i << ",ndf=0" << std::endl;  
                // If you want OM numbers, pass it in as an argument
        }
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

    TTree *qaTree = new TTree("qaTree", "FEB/Wavecatcher QA");

    // Initialise output variables
    int event_num = -1;
    int om_num = -1;
    double baseline_orig = 0;
    double stddev_orig = 0;
    double eom_orig = 0;
    double end_baseline = 0;
    double end_stddev = 0;
    double end_eom = 0;

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

    // FFT variables
    double baseline_fft = 0;
    double stddev_fft = 0;
    double eom_fft = 0;

    //Quality assurance variables
    int qa_om_num = -1;
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
    baseline_tree->Branch("end_baseline", &end_baseline, "end_baseline/D");
    baseline_tree->Branch("end_stddev", &end_stddev, "end_stddev/D");
    baseline_tree->Branch("end_eom", &end_eom, "end_eom/D");

    baseline_tree->Branch("baseline_feb", &baseline_feb, "baseline_feb/D");
    baseline_tree->Branch("stddev_feb", &stddev_feb, "stddev_feb/D");
    baseline_tree->Branch("eom_feb", &eom_feb, "eom_feb/D");

    baseline_tree->Branch("baseline_976", &baseline_976, "baseline_976/D");
    baseline_tree->Branch("stddev_976", &stddev_976, "stddev_976/D");
    baseline_tree->Branch("eom_976", &eom_976, "eom_976/D");

    baseline_tree->Branch("baseline_end", &baseline_end, "baseline_end/D");
    baseline_tree->Branch("stddev_end", &stddev_end, "stddev_end/D");
    baseline_tree->Branch("eom_end", &eom_end, "eom_end/D");

    baseline_tree->Branch("baseline_fft", &baseline_fft, "baseline_fft/D");
    baseline_tree->Branch("stddev_fft", &stddev_fft, "stddev_fft/D");
    baseline_tree->Branch("eom_fft", &eom_fft, "eom_fft/D");

    qaTree->Branch("om_num", &qa_om_num, "om_num/I");
    qaTree->Branch("mean_mem", &mean_mem);
    qaTree->Branch("sigma_mem", &sigma_mem);
    qaTree->Branch("chi2ndf_mem", &chi2ndf_mem);    
    
    qaTree->Branch("mean_timeo", &mean_timeo);
    qaTree->Branch("sigma_timeo", &sigma_timeo);
    qaTree->Branch("chi2ndf_timeo", &chi2ndf_timeo);

    //////////////////////////// MAPS ////////////////////////////////////////
    std::map<int, std::vector<std::vector<float>>> adc_values;
    std::map<int, std::vector<double>> mod16_offsets_per_om;
    std::map<int, std::vector<double>> mod16_sigmas_per_om;
    std::map<int, std::vector<double>> mod16_chi2ndf_per_om;
    std::map<int, std::vector<double>> mod16_data_per_om;

    std::map<int, std::vector<std::vector<float>>> feb_adc_values;
    std::map<int, std::vector<double>> feb_offsets_per_om;
    std::map<int, std::vector<double>> feb_sigmas_per_om;
    std::map<int, std::vector<double>> feb_chi2ndf_per_om;

    std::map<int, std::vector<double>> eow_offsets_per_om;
    std::map<int, std::vector<double>> eow_sigmas_per_om;
    std::map<int, std::vector<double>> eow_chi2ndf_per_om;
    std::map<int, std::vector<double>> eow_gausses_per_om;

    static std::map<int, std::vector<std::vector<double>>> om_waveforms;
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

    ////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////// FEB Calclulations /////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////
    std::ofstream csv_file_feb("feb_offsets.csv");
    csv_file_feb << "om,bin,data_mean,fit_mean,sigma,chi2ndf\n";

    // Process FEB_ADC_VALUES for histogramming
    for (auto& om_pair : feb_adc_values) {
        int om_num = om_pair.first;

        // if (om_num != 1) continue; // Skip all OMs except OM 1 for now -----------------
        // TCanvas *c = new TCanvas("c", "OM 1 Histograms", 800, 600); //------------------

        auto& adc_bins = om_pair.second;

        std::vector<double> fit_means_feb_om; // filled with 976 Gaussian means
        std::vector<double> fit_sigmas_feb_om; // filled with 976 Gaussian sigmas
        std::vector<double> fit_chi2_ndf_feb_om; // filled with 976 Gaussian chi2/ndf values
        std::vector<double> data_means_feb_om; // filled with 976 data means

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
            float min_adc = *std::min_element(valid_adc_values.begin(), valid_adc_values.end()) - 2;
            float max_adc = *std::max_element(valid_adc_values.begin(), valid_adc_values.end()) + 2;
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));

            TString hist_name = Form("feb_om%d_bin%d", om_num, bin);
            TString hist_title = Form("FEB OM %d - Bin %d;ADC Value;Counts", om_num, bin);

            TH1D *hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);

            for (float adc_value : valid_adc_values) {
                hist->Fill(adc_value);
            }

            // extract data mean
            double data_mean = hist->GetMean();
            data_means_feb_om.push_back(data_mean);

            // Fit gaussian, extract mean/sigma
            TF1 *fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);

            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = -1;
            if (ndf !=0){
                chi2_ndf = chi2 / ndf;
            } else {
                std::cout << "[BAD FIT] ndf=0 for OM=" << om_num
                        << " Bin=" << bin << std::endl;
            }

            if (chi2_ndf > 9) {
                std::cout << "[WARNING] High chi2/ndf in FEB fit: "
                        << "OM = " << om_num
                        << ", Bin = " << bin
                        << ", chi2/ndf = " << chi2_ndf << std::endl;

                // Create canvas and draw histogram
                TCanvas *c = new TCanvas("c", "High chi2/ndf Fit", 800, 600);
                gStyle->SetOptStat(0);  // Disable stats box
                hist->Draw();
                fit_func->Draw("same");

                // Save to PNG with unique name
                TString filename = Form("feb_fit_om%d_bin%d_chi2%.2f.png", om_num, bin, chi2_ndf);
                c->SaveAs(filename);

                delete c;  // Clean up canvas
            }

            fit_means_feb_om.push_back(mean);
            fit_sigmas_feb_om.push_back(sigma);
            fit_chi2_ndf_feb_om.push_back(chi2_ndf);

            csv_file_feb << om_num << ","
                        << bin << ","
                        << data_mean << ","
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

    // //////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////////// Wavecater Calclulations ///////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////

    // Process ADC_VALUES for histogramming, fitting and mod16 calculations
    for (auto& om_pair : adc_values) {
        int om_num = om_pair.first;
        auto& adc_bins = om_pair.second;

        std::vector<double> fit_means_om; // filled with 976 Gaussian means
        std::vector<double> fit_sigmas_om; // filled with 976 Gaussian sigmas
        std::vector<double> fit_chi2_ndf_om; // filled with 976 Gaussian chi2/ndf values
        std::vector<double> data_means_om; // filled with 976 data means

        for (int bin = 0; bin < 976; ++bin){
            std::vector<float>& adc_values_bin = adc_bins[bin];

            // Create and fill histogram
            TString hist_name = Form("om%d_bin%d", om_num, bin);
            TString hist_title = Form("OM %d - Bin %d;ADC Value;Counts", om_num, bin);
            // set histogram range based on the ADC values per om
            float min_adc = *std::min_element(adc_values_bin.begin(), adc_values_bin.end())-2;
            float max_adc = *std::max_element(adc_values_bin.begin(), adc_values_bin.end())+2;
            // Ensure bin count is an integer:
            int nbins = static_cast<int>(std::ceil(max_adc - min_adc + 1));
            // create histogram with nbins between min and max ADC values
            TH1D *hist = new TH1D(hist_name, hist_title, nbins, min_adc, max_adc);
            for (float adc_value : adc_values_bin) {
                hist->Fill(adc_value);
            }
            // extract data mean
            double data_mean = hist->GetMean();
            data_means_om.push_back(data_mean);

            // Fit gaussian, extract mean/sigma
            TF1 *fit_func = new TF1("fit_func", "gaus", min_adc, max_adc);
            hist->Fit(fit_func, "Q", "", min_adc, max_adc);
            double mean = fit_func->GetParameter(1);
            double sigma = fit_func->GetParameter(2);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = -1;
            if (ndf !=0){
                chi2_ndf = chi2 / ndf;
            } else {
                std::cout << "[BAD FIT] ndf=0 for OM=" << om_num
                        << " Bin=" << bin << std::endl;
            }
            
            ////////////////////////////////////////// Print interesting chi2 ////////////////////////////////////////////
            static int saved_count = 0;
            if (om_num == 652 && chi2_ndf > 9 && saved_count < 10) {
                std::cout << "[WARNING] High chi2/ndf in FEB fit: "
                        << "OM = " << om_num
                        << ", Bin = " << bin
                        << ", chi2/ndf = " << chi2_ndf << std::endl;

                // Create canvas and draw histogram
                TCanvas *c = new TCanvas("c", "High chi2/ndf Fit", 800, 600);
                hist->Draw();
                fit_func->Draw("same");

                // Add chi2/ndf text label in top-right corner
                TLatex latex;
                latex.SetNDC(); // use normalized device coords (0–1)
                latex.SetTextSize(0.04);
                latex.DrawLatex(0.65, 0.85, Form("#chi^{2}/ndf = %.2f", chi2_ndf));

                // Save to PNG with unique name
                TString filename = Form("wcatch_fit_om%d_bin%d_chi2%.2f.png", om_num, bin, chi2_ndf);
                c->SaveAs(filename);

                delete c;  // Clean up canvas
                saved_count++;
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
            
            fit_means_om.push_back(mean);
            fit_sigmas_om.push_back(sigma);
            fit_chi2_ndf_om.push_back(chi2_ndf);

            delete hist;
            delete fit_func;

        }

        // Calculate mod16 offsets for this OM
        Mod16Stats stats = compute_mod16_from_fits(fit_means_om, fit_sigmas_om, fit_chi2_ndf_om);
        mod16_offsets_per_om[om_num] = stats.offsets;

        // calculate mod 16 data means for this OM
        mod16_data_per_om[om_num] = std::vector<double>(16, 0.0);   
        std::vector<int> counts(16, 0);
        for (int i = 0; i < 976; ++i) {
            int mod16_index = i % 16;
            mod16_data_per_om[om_num][mod16_index] += data_means_om[i];
            counts[mod16_index]++;
        }
        for (int i = 0; i < 16; ++i) {
            if (counts[i] > 0) {
                mod16_data_per_om[om_num][i] /= counts[i];
            }
        }

        // optionally store sigmas and chi2ndf too:
        mod16_sigmas_per_om[om_num] = stats.sigmas;
        mod16_chi2ndf_per_om[om_num] = stats.chi2ndfs;
    }
    // Open CSV file for writing
    std::ofstream csv_file("mod16_offsets.csv");

    // Write header
    csv_file << "om,offset_index,data_mean,offset_value,sigma_value,chi2ndf_value\n";

    // Loop through the map and write data
    for (const auto& pair : mod16_offsets_per_om) {
        int om = pair.first;
        const std::vector<double>& offsets = pair.second;
        for (int i = 0; i < (int)offsets.size(); ++i) {
            csv_file << om << ","
                 << i << ","
                 << mod16_data_per_om[om][i] << ","
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
        std::vector<double> spike_gaus_om;

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
            double chi2_ndf = -1;
              if (ndf !=0){
                chi2_ndf = chi2 / ndf;
            } else {
                std::cout << "[BAD FIT] ndf=0 for OM=" << om_num
                        << " Bin=" << bin << std::endl;
            }

            // Calculate difference: spike mean - mod16 offset for (bin % 16)
            int mod16_index = bin % 16;
            double mod16_value = (mod16_index < (int)mod16_offsets.size()) ? mod16_offsets[mod16_index] : 0.0;
            double corrected_spike_offset = mean - mod16_value;

            // Save results directly instead of computing mod16 offsets
            spike_gaus_om.push_back(mean);
            spike_means_om.push_back(corrected_spike_offset);
            spike_sigmas_om.push_back(sigma);
            spike_chi2ndf_om.push_back(chi2_ndf);

            delete hist;
            delete fit_func;
        }
        // Store the results for this OM
        eow_gausses_per_om[om_num] = spike_gaus_om;
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
            end_baseline = calculate_baseline48(wave->at(j));
            end_stddev = calculate_stddev48(wave->at(j), end_baseline);
            end_eom = calculate_eom48(wave->at(j), end_baseline);

            // ////////////////////////////////////// FOURIER NOTCH FILTERING ///////////////////////////////////////
            // // Convert waveform to double precision for FFT input
            // std::vector<double> waveform_double(1024);
            // for (int bin = 0; bin < 1024; ++bin) {
            //     waveform_double[bin] = static_cast<double>(wave->at(j)[bin]);
            // }

            // // Convert to double for FFT
            // int N = 1024;
            // double fs = 3.2e9;  // sampling freq Hz
            // int n_freqs = N / 2 + 1;
            

            // // FFT
            // TVirtualFFT::SetTransform(0);
            // TVirtualFFT* fft = TVirtualFFT::FFT(1, &N, "R2C EX K");
            // fft->SetPoints(&waveform_double[0]);
            // fft->Transform();

            // std::vector<TComplex> freq_domain_before(n_freqs);
            // for (int k = 0; k < n_freqs; ++k) {
            //     double re, im;
            //     fft->GetPointComplex(k, re, im);
            //     freq_domain_before[k] = TComplex(re, im);
            // }

            // // Notch filter with broadened bins (±1 bin)
            // std::vector<TComplex> freq_domain = freq_domain_before;
            // std::vector<double> notch_freqs = { 201, 601, 801, 1001, 1201, 1401};
            // std::vector<int> notch_bins;
            // for (double f : notch_freqs) {
            //     int bin = static_cast<int>(round(f * 1e6 * N / fs));
            //     if (bin >= 0 && bin < n_freqs) notch_bins.push_back(bin);
            // }

            // for (int bin : notch_bins) {
            //     for (int b = bin - 0.5; b <= bin + 0.5; ++b) {
            //         if (b >= 0 && b < n_freqs) {
            //             freq_domain[b] = TComplex(0.0, 0.0);
            //             // Zero symmetric bin for real FFT except DC and Nyquist
            //             if (b != 0 && b != n_freqs - 1) {
            //                 int sym_bin = N - b;
            //                 if (sym_bin < n_freqs) freq_domain[sym_bin] = TComplex(0.0, 0.0);
            //             }
            //         }
            //     }
            // }

            // // Inverse FFT
            // TVirtualFFT* ifft = TVirtualFFT::FFT(1, &N, "C2R EX K");
            // for (int k = 0; k < n_freqs; ++k) {
            //     ifft->SetPointComplex(k, freq_domain[k]);
            // }
            // ifft->Transform();

            // double* filtered = new double[N];
            // ifft->GetPoints(filtered);

            // // Normalize filtered waveform
            // std::vector<double> filtered_wave(N);
            // for (int bin = 0; bin < N; ++bin) {
            //     filtered_wave[bin] = filtered[bin] / N;
            // }

            // delete[] filtered;
            // delete fft;
            // delete ifft;

            // // Convert the filtered waveform back to short integers
            // for (int bin = 0; bin < 1024; ++bin) {
            //     wave->at(j)[bin] = static_cast<short>(std::round(filtered_wave[bin]));
            // }

            //////////////////////////////////////////////////////////////////////////////////////////

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

            // Static counter to track saved FEB-corrected waveforms
            static int saved_feb_waveforms = 0;

            if (saved_feb_waveforms < 20) {
                int n_bins = 976;

                // Copy original waveform bins before correction
                std::vector<double> original_waveform_bins(wave->at(j).begin(), wave->at(j).begin() + n_bins);

                // Create canvas and divide into 2 pads (1 row, 2 columns)
                TCanvas* c_feb = new TCanvas("c_feb", "FEB Correction Region", 1200, 400);
                c_feb->Divide(2, 1);

                // Left pad: original waveform (before correction)
                c_feb->cd(1);
                TGraph* g_before = new TGraph(n_bins);
                for (int bin = 0; bin < n_bins; ++bin) {
                    g_before->SetPoint(bin, bin, original_waveform_bins[bin]);
                }
                g_before->SetLineColor(kBlue);
                g_before->SetMarkerSize(0.5);
                g_before->SetTitle("Bins 0-975 Before FEB Correction;Bin;ADC Count");
                g_before->Draw("APL");

                // Right pad: corrected waveform
                c_feb->cd(2);
                TGraph* g_after = new TGraph(n_bins);
                for (int bin = 0; bin < n_bins; ++bin) {
                    g_after->SetPoint(bin, bin, corrected_wave[bin]);
                }
                g_after->SetLineColor(kRed);
                g_after->SetMarkerSize(0.5);
                g_after->SetTitle("Bins 0-975 After FEB Correction;Bin;ADC Count");
                g_after->Draw("APL");

                // Save plot as PNG, including event and waveform IDs
                std::string filename = "waveform_feb_comparison_event" + std::to_string(i) + "_waveform" + std::to_string(j) + ".png";
                c_feb->SaveAs(filename.c_str());

                // Clean up
                delete g_before;
                delete g_after;
                delete c_feb;

                saved_feb_waveforms++;
            }


            // Step 4: calculate baseline, stddev, and eom on corrected waveform in original order
            baseline_feb = calculate_baseline976(corrected_wave);
            stddev_feb   = calculate_stddev976(corrected_wave, baseline_feb);
            eom_feb      = calculate_eom976(corrected_wave, baseline_feb);

            // 976 waveform correction ///////////////////////////////////////////////////////////////

            // --- Save original waveform BEFORE applying correction ---
            std::vector<double> original_waveform_bins(wave->at(j).begin(), wave->at(j).begin() + 976);

            // --- Apply Mod16 correction ---
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

            // --- Plot before/after comparison ---
            static int saved_mod16_waveforms = 0;

            if (saved_mod16_waveforms < 20) {
                int n_bins = 976;

                // Create canvas and divide into 2 pads (1 row, 2 columns)
                TCanvas* c_mod16 = new TCanvas("c_mod16", "Mod16 Correction Region", 1200, 400);
                c_mod16->Divide(2, 1);

                // Draw original (before correction) graph on left pad
                c_mod16->cd(1);
                TGraph* g_before = new TGraph(n_bins);
                for (int bin = 0; bin < n_bins; ++bin) {
                    g_before->SetPoint(bin, bin, original_waveform_bins[bin]);
                }
                g_before->SetLineColor(kBlue);
                g_before->SetMarkerSize(0.5);
                g_before->SetTitle("Bins 0-975 Before Mod16 Correction;Bin;ADC Count");
                g_before->Draw("APL");

                // Draw corrected waveform on right pad
                c_mod16->cd(2);
                TGraph* g_after = new TGraph(n_bins);
                for (int bin = 0; bin < n_bins; ++bin) {
                    g_after->SetPoint(bin, bin, wave->at(j)[bin]);
                }
                g_after->SetLineColor(kRed);
                g_after->SetMarkerSize(0.5);
                g_after->SetTitle("Bins 0-975 After Mod16 Correction;Bin;ADC Count");
                g_after->Draw("APL");

                // Save plot as PNG
                std::string filename = "waveform_mod16_comparison_event" + std::to_string(i) + "_waveform" + std::to_string(j) + ".png";
                c_mod16->SaveAs(filename.c_str());

                // Clean up
                delete g_before;
                delete g_after;
                delete c_mod16;

                saved_mod16_waveforms++;
            }

            // --- Calculate adjusted baseline, stddev, and eom ---
            baseline_976 = calculate_baseline976(wave->at(j));
            stddev_976 = calculate_stddev976(wave->at(j), baseline_976);
            eom_976 = calculate_eom976(wave->at(j), baseline_976);

         
            // Spike correction for bins 976–1023 /////////////////////////////////////////////////
            static int saved_wave_snippets = 0;

            int start_bin = 976;
            int end_bin = 1023;
            int n_bins = end_bin - start_bin + 1;

            // Copy original bins BEFORE correction
            std::vector<double> original_spike_bins(wave->at(j).begin() + start_bin, wave->at(j).begin() + end_bin + 1);

            // Apply spike + mod16 correction
            auto it_spike = eow_offsets_per_om.find(om_num);

            if (it_spike != eow_offsets_per_om.end()) {
                const std::vector<double>& spike_offsets = it_spike->second;

                for (int bin = start_bin; bin <= end_bin; ++bin) {
                    int spike_bin_index = bin - start_bin;
                    if (spike_bin_index < (int)spike_offsets.size()) {
                        wave->at(j)[bin] -= spike_offsets[spike_bin_index]; // remove spike offset

                        // Apply mod16 offset to spike-corrected bin
                        auto it_mod16 = mod16_offsets_per_om.find(om_num);  
                        if (it_mod16 != mod16_offsets_per_om.end()) {
                            const std::vector<double>& offsets = it_mod16->second;
                            wave->at(j)[bin] -= offsets[bin % 16]; // apply mod16 offset
                        }
                    }
                }
            } else {
                std::cerr << "Warning: No spike offsets found for OM " << om_num << std::endl;
            }

            // Plot comparison
            if (saved_wave_snippets < 20) {
                // Create canvas with 2 pads
                TCanvas* c_spike = new TCanvas("c_spike", "Spike Region Correction", 1200, 400);
                c_spike->Divide(2, 1);

                // Draw original (before correction)
                c_spike->cd(1);
                TGraph* g_before = new TGraph(n_bins);
                for (int i = 0; i < n_bins; ++i) {
                    int bin = start_bin + i;
                    g_before->SetPoint(i, bin, original_spike_bins[i]);
                }
                g_before->SetLineColor(kBlue);
                g_before->SetMarkerSize(0.5);
                g_before->SetTitle(Form("Bins %d-%d Before Spike Correction;Bin;ADC Count", start_bin, end_bin));
                g_before->Draw("APL");

                // Draw corrected (after correction)
                c_spike->cd(2);
                TGraph* g_after = new TGraph(n_bins);
                for (int i = 0; i < n_bins; ++i) {
                    int bin = start_bin + i;
                    g_after->SetPoint(i, bin, wave->at(j)[bin]);
                }
                g_after->SetLineColor(kRed);
                g_after->SetMarkerSize(0.5);
                g_after->SetTitle(Form("Bins %d-%d After Spike Correction;Bin;ADC Count", start_bin, end_bin));
                g_after->Draw("APL");

                // Save as PNG
                std::string filename = "waveform_spike_comparison_event" + std::to_string(i) + "_waveform" + std::to_string(j) + ".png";
                c_spike->SaveAs(filename.c_str());

                // Clean up
                delete g_before;
                delete g_after;
                delete c_spike;

                saved_wave_snippets++;
            }

            // Recalculate baseline, stddev, and eom
            baseline_end = calculate_baseline48(wave->at(j));
            stddev_end   = calculate_stddev48(wave->at(j), baseline_end);
            eom_end      = calculate_eom48(wave->at(j), baseline_end);
            
            ////////////////////////// Fill the Baseline Tree ////////////////////////////
            baseline_tree->Fill();
        } 
    }
    
    // storing quality of fit assurance data ////////////////////////////////////////////////////
    for (int omnum = 0; omnum < 712; ++omnum) {
        qa_om_num = omnum;
        // Clear previous event's data
        mean_timeo.clear();
        sigma_timeo.clear();
        chi2ndf_timeo.clear();

        auto it_offset = mod16_offsets_per_om.find(omnum);
        auto it_sigma  = mod16_sigmas_per_om.find(omnum);
        auto it_chi2   = mod16_chi2ndf_per_om.find(omnum);

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

        auto it_feb_offset = feb_offsets_per_om.find(omnum);
        auto it_feb_sigma  = feb_sigmas_per_om.find(omnum);
        auto it_feb_chi2   = feb_chi2ndf_per_om.find(omnum);

        if (it_feb_offset != feb_offsets_per_om.end() &&
            it_feb_sigma != feb_sigmas_per_om.end() &&
            it_feb_chi2 != feb_chi2ndf_per_om.end()) {

            mean_mem = it_feb_offset->second;
            sigma_mem = it_feb_sigma->second;
            chi2ndf_mem = it_feb_chi2->second;

        } else {
            std::cerr << "Warning: Missing FEB QA data for OM " << omnum << std::endl;

            // Optionally fill with NaNs or -999 if QA data missing
            mean_mem.assign(1024, -999);
            sigma_mem.assign(1024, -999);
            chi2ndf_mem.assign(1024, -999);

        }
        qaTree->Fill();
    }   
    //////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// End of code, writes and closes ////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // write and close
    outfile->cd();
    baseline_tree->Write();
    qaTree->Write();
    outfile->Close();
    file->Close();

    return 0;
}
