#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

void treecollect_waveform(const vector<int>& e_events_to_plot) {

    // Open the existing ROOT file with the original tree
    TFile *file = TFile::Open("filtered_gamma_e.root", "READ");  // input file path
    if (!file || file->IsZombie()) {
        cout << "Error: Could not open input file!" << endl;
        return;
    }

    // Access the original tree
    TTree *tree = (TTree*)file->Get("Result_tree");
    if (!tree) {
        cout << "Error: Could not find tree in file!" << endl;
        file->Close();
        return;
    }

    // Set up variables to read the branches of the original tree
    std::vector<int> *om_number = nullptr;
    std::vector<int> *electron = nullptr;
    std::vector<double> *energy = nullptr;
    std::vector<int> *charge = nullptr;
    std::vector<std::vector<short>> *wave = nullptr;

    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchAddress("electron", &electron);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("wave", &wave);

    TFile *newFile = new TFile("treebuild.root", "RECREATE");
    TTree *newTree = new TTree("treebuild", "treebuild");

    // Declare electron parameters
    int e_event = -1;
    int e_om_number = 0;
    double e_energy = 0;
    int e_charge = 0;
    double e_baseline = 0;
    double e_wave_peak_value = 0;
    int e_wave_peak_time = 0;
    double e_calc_amplitude = 0;

    // Fit parameters
    double e_chi2 = 0;
    double e_ndf = 0;
    double e_chi2ndf = 0;
    double e_fit_a_factor = 0;
    double e_fit_pulse_onset = 0;
    double e_fit_r_factor = 0;
    double e_fit_d_factor = 0;

    // Branches for electron data
    newTree->Branch("e_event", &e_event);
    newTree->Branch("e_om_number", &e_om_number);
    newTree->Branch("e_energy", &e_energy);
    newTree->Branch("e_charge", &e_charge);
    newTree->Branch("e_baseline", &e_baseline);
    newTree->Branch("e_wave_peak_value", &e_wave_peak_value);
    newTree->Branch("e_wave_peak_time", &e_wave_peak_time);
    newTree->Branch("e_calc_amplitude", &e_calc_amplitude);
    newTree->Branch("e_chi2", &e_chi2);
    newTree->Branch("e_ndf", &e_ndf);
    newTree->Branch("e_chi2ndf", &e_chi2ndf);
    newTree->Branch("e_fit_a_factor", &e_fit_a_factor);
    newTree->Branch("e_fit_pulse_onset", &e_fit_pulse_onset);
    newTree->Branch("e_fit_r_factor", &e_fit_r_factor);
    newTree->Branch("e_fit_d_factor", &e_fit_d_factor);

    Long64_t nEntries = tree->GetEntries();
    int ticks = 1024;

    for (Long64_t j = 0; j < nEntries; j++) {
        tree->GetEntry(j);
        
        for (size_t i = 0; i < electron->size(); i++) {
            if (electron->at(i) == 1) {  // Electron only
                
                // Baseline calculation
                e_baseline = 0;
                for (int k = 1; k < 97; k++) {
                    e_baseline += wave->at(i).at(k);
                }
                e_baseline /= 96;

                // Standard deviation calculation
                double e_sum_sq_diff = 0;
                for (int k = 1; k < 97; k++) {
                    e_sum_sq_diff += pow(wave->at(i).at(k) - e_baseline, 2);
                }
                double e_std_dev = sqrt(e_sum_sq_diff / 95);

                // Find waveform peak
                e_wave_peak_value = wave->at(i).at(0);
                e_wave_peak_time = 0;
                for (int k = 1; k < wave->at(i).size(); k++) {
                    if (wave->at(i).at(k) < e_wave_peak_value) {
                        e_wave_peak_value = wave->at(i).at(k);
                        e_wave_peak_time = k;
                    }
                }
                e_calc_amplitude = e_wave_peak_value - e_baseline;

                // Histogram for electron waveform
                TH1D *e_h_wave = new TH1D(Form("e_h_wave_%lld", j), Form("Waveform Histogram for Electron Event %lld", j), ticks, 0, ticks);
                for (int k = 0; k < wave->at(i).size(); k++) {
                    e_h_wave->SetBinContent(k+1, wave->at(i).at(k));
                    e_h_wave->SetBinError(k+1, 1.6);
                }

                // Define fit function for waveform
                TF1 *e_fitFunc = new TF1("e_fitFunc", Form("(1/(1-exp(-[4]*(x-[1]))))*[0]*(exp(-(x-[1])/[2]) - exp(-(x-[1])/[3])) + %f", e_baseline), 200, 900);
                e_fitFunc->SetParameters(e_calc_amplitude, 266, 18, 10, 5);

                // Fit and retrieve fit parameters
                e_h_wave->Fit("e_fitFunc", "RQ");
                e_chi2 = e_fitFunc->GetChisquare();
                e_ndf = e_fitFunc->GetNDF();
                e_chi2ndf = e_chi2 / e_ndf;
                e_fit_a_factor = e_fitFunc->GetParameter(0);
                e_fit_pulse_onset = e_fitFunc->GetParameter(1);
                e_fit_r_factor = e_fitFunc->GetParameter(2);
                e_fit_d_factor = e_fitFunc->GetParameter(3);

                // Fill new tree
                e_event = j;
                e_om_number = om_number->at(i);
                e_energy = energy->at(i);
                e_charge = charge->at(i);
                newTree->Fill();

                // Plot if event is in e_events_to_plot
                if (std::find(e_events_to_plot.begin(), e_events_to_plot.end(), j) != e_events_to_plot.end()) {
                    TCanvas *c = new TCanvas(Form("c_event_%lld", j), Form("Electron Event %lld", j), 800, 600);
                    e_h_wave->SetTitle(Form("Waveform and Fit for Electron Event %lld", j));
                    e_h_wave->Draw();
                    e_fitFunc->Draw("SAME");

                    // Add legend in bottom right corner
                    TLegend *leg = new TLegend(0.7, 0.1, 0.9, 0.3);
                    leg->AddEntry(e_h_wave, "Waveform", "l");
                    leg->AddEntry(e_fitFunc, "Fit", "l");
                    leg->AddEntry((TObject*)0, Form("Amplitude: %.2f", e_calc_amplitude), "");
                    leg->Draw();  // Avoid "same" in Draw() for the legend
                    
                    // Update the canvas to ensure the legend is displayed
                    c->Update();
                    
                    // Save and clean up
                    c->SaveAs(Form("treecollect_waveform_Event_%lld.png", j));
                    delete c;
                }


                delete e_h_wave;
            }
        }
    }

    newTree->Write();
    file->Close();
}

// Main function to call treebuilding with selected events to plot
int main() {
    vector<int> e_events_to_plot = {7, 9, 20};  // Add specific electron events you want to plot
    treecollect_waveform(e_events_to_plot);
    return 0;
}
