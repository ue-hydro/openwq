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

void OpenWQ_SI_model::langmuir(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){



}



/*

double calculateSorption(double qmax, double KL, double concentration) {
    return (qmax * KL * concentration) / (1 + KL * concentration);
}

// Function to calculate equilibrium concentration for a given sorbed amount
double calculateDesorption(double qmax, double KL, double sorbedAmount) {
    return (sorbedAmount / (qmax - sorbedAmount)) / KL;
}

int main() {
    // Langmuir parameters
    double qmax, KL;
    cout << "Enter Langmuir parameter qmax (maximum adsorption capacity): ";
    cin >> qmax;
    cout << "Enter Langmuir constant KL: ";
    cin >> KL;

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
        double sorbedAmount = calculateSorption(qmax, KL, C);
        cout << setw(15) << C << setw(15) << sorbedAmount << endl;
    }

    // Example for desorption (optional)
    cout << "\nEnter sorbed amount to calculate desorption: ";
    double sorbedAmount;
    cin >> sorbedAmount;
    double desorbedConcentration = calculateDesorption(qmax, KL, sorbedAmount);
    cout << "Desorbed concentration: " << desorbedConcentration << endl;

    return 0;
}
*/