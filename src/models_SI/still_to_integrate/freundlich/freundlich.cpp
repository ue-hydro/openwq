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
