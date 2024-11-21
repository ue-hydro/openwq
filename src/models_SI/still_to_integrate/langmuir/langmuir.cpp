#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

// Function to calculate sorption based on the Langmuir isotherm
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
