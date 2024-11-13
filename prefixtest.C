#include <vector>
#include <string>
#include <iostream>

void prefixtest() {

    // Initialize vectors to hold the data with 8 entries each (4 for e_ and 4 for g_)
    std::vector<double> energy = {10.5, 20.2, 30.5, 40.3, 15.7, 25.4, 35.1, 45.6};  // e_ and g_ respectively
    std::vector<double> charge = {1.1, 2.2, 3.3, 4.4, 1.5, 2.6, 3.7, 4.8};    // e_ and g_ respectively
    std::vector<int> om_number = {101, 102, 103, 104, 105, 106, 107, 108};    // e_ and g_ respectively
    std::vector<int> event = {201, 201, 202, 202, 203, 203, 204, 204};        // e_ and g_ respectively

    // Vector to store the electron indicator (true for e_, false for g_)
    std::vector<bool> electron_indicator;

    // Open the input file to read the boolean indicators
    std::ifstream input_file("prefix_e_vector.txt");
    
    // Read the 1D vector from the file (1 or 0)
    int value;
    while (input_file >> value) {
        if (value == 1 || value == 0) {
            electron_indicator.push_back(value);
        } else {
            std::cerr << "Invalid value in input file. Only 1 and 0 are allowed." << std::endl;
            return;
        }
    }

    // Ensure the size of prefix_indicator matches the size of other data vectors
    if (electron_indicator.size() != energy.size()) {
        std::cerr << "Mismatch in data size between prefix_indicator and other data vectors." << std::endl;
        return;
    }
    // Vectors to store the formatted output strings
    std::vector<std::string> electron_vars;
    std::vector<std::string> gamma_vars;

    // Loop through the data
    for (int i = 0; i < electron_indicator.size(); i++) {

        // Determine the prefix based on the value in electron_indicator
        std::string prefix = electron_indicator[i] ? "e_" : "g_";

        // Format each variable with its corresponding value
        std::string energy_str = prefix + "energy: " + std::to_string(energy[i]);
        std::string charge_str = prefix + "charge: " + std::to_string(charge[i]);
        std::string om_number_str = prefix + "om_number: " + std::to_string(om_number[i]);
        std::string event_str = prefix + "event: " + std::to_string(event[i]);
        
        // Store each formatted string in the appropriate vector
        if (electron_indicator[i]) {
            electron_vars.push_back(event_str);
            electron_vars.push_back(energy_str);
            electron_vars.push_back(charge_str);
            electron_vars.push_back(om_number_str);
            
        } else {
            gamma_vars.push_back(event_str);
            gamma_vars.push_back(energy_str);
            gamma_vars.push_back(charge_str);
            gamma_vars.push_back(om_number_str);
            
        }
    }

    // Print the electron-prefixed data
    std::cout << "Electron Variables (e_):" << std::endl;
    for (const auto& var : electron_vars) {
        std::cout << var << std::endl;
    }

    // Print the gamma-state-prefixed data
    std::cout << "\ngamma State Variables (g_):" << std::endl;
    for (const auto& var : gamma_vars) {
        std::cout << var << std::endl;
    }
}
============================ included for later, but not needed now =============================
// Main function to run the test
int main() {
    prefixtest();
    return 0;
}
