
// Copyright 2020, Diogo Costa, diogo.costa@uevora.pt
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "models_SI/headerfile_SI.hpp"

// Function to calculate sorption based on the Langmuir isotherm

void OpenWQ_SI_model::freundlich(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){

    // TO REMOVE -> only testing the cpassing of relevant variable lass-nested/nested/nested variable 
    OpenWQ_wqconfig.SI_model->FREUNDLICH->test_parameter_freundlich_2_delete = 8;

    /*

    #include <iostream>
    #include <cmath>
    #include <vector>
    #include <iomanip>

    using namespace std;

    // Function to calculate sorption based on the Freundlich model
    double calculateSorption(double Kf, double n, double concentration) {
        return Kf * pow(concentration, n);
    }

    // Function to calculate desorption
    // (Inverse calculation if required; for illustration, assumes desorption is the same process)
    double calculateDesorption(double Kf, double n, double sorbedAmount) {
        return pow(sorbedAmount / Kf, 1.0 / n);
    }

    int main() {
        // Freundlich parameters
        double Kf, n;
        cout << "Enter Freundlich constant Kf: ";
        cin >> Kf;
        cout << "Enter Freundlich constant n: ";
        cin >> n;

        // Input equilibrium concentrations
        vector<double> concentrations;
        cout << "Enter equilibrium concentrations (space-separated, end with -1): ";
        double concentration;
        while (cin >> concentration && concentration != -1) {
            concentrations.push_back(concentration);
        }

        // Calculate sorption for each concentration
        cout << "\nSorption Results:\n";
        cout << setw(15) << "Concentration" << setw(15) << "Sorption" << endl;
        for (double C : concentrations) {
            double sorbedAmount = calculateSorption(Kf, n, C);
            cout << setw(15) << C << setw(15) << sorbedAmount << endl;
        }

        // Example for desorption (optional)
        cout << "\nEnter sorbed amount to calculate desorption: ";
        double sorbedAmount;
        cin >> sorbedAmount;
        double desorbedConcentration = calculateDesorption(Kf, n, sorbedAmount);
        cout << "Desorbed concentration: " << desorbedConcentration << endl;

        return 0;
    }

    */

}
