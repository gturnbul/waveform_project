#include <vector>
#include <string>
#include <iostream>

void prefixtest() {
    // Vector where true represents "e_" prefix and false represents "g_" prefix
    std::vector<bool> electron_indicator = {true, false};

    // Corresponding values for each attribute
    std::vector<double> energy = {10.5, 20.2};  // e_ and g_ respectively
    std::vector<double> charge = {1.1, 2.2};    // e_ and g_ respectively
    std::vector<int> om_number = {101, 102};    // e_ and g_ respectively
    std::vector<int> event = {201, 202};        // e_ and g_ respectively

    // Vectors to store the formatted output strings
    std::vector<std::string> electron_vars;
    std::vector<std::string> gamma_vars;

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
            electron_vars.push_back(energy_str);
            electron_vars.push_back(charge_str);
            electron_vars.push_back(om_number_str);
            electron_vars.push_back(event_str);
        } else {
            gamma_vars.push_back(energy_str);
            gamma_vars.push_back(charge_str);
            gamma_vars.push_back(om_number_str);
            gamma_vars.push_back(event_str);
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

// Main function to run the test
int main() {
    prefixtest();
    return 0;
}
