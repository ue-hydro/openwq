
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

#include <cmath>
#include <algorithm>


 // Calculate eroded particles from soil with a method based on the Morgan-Morgan-Finney erosion model. (from HYPE)
//
// Reference: R. P. C. Morgan, D. D. V. Morgan, and H. J. Finney, 1984. 
// A predictive model for the assessment of soil erosion risk, Journal 
// of Agricultural Engineering Research, vol. 30, pp. 245�253
// 


// Function to calculate current crop and ground cover (stub, needs actual implementation)
void get_current_crop_and_ground_cover(int i, int j, double& common_cropcover, double& common_groundcover) {
    // Placeholder: Set default values for testing
    common_cropcover = 0.5;
    common_groundcover = 0.3;
}

void hype_mmf(int i, int j, int pdayno, double prec, double surfacerunoff, double cohesion, 
                           double erodibility, double snow, double sreroexp, double flow, double& erodedsed) {

    // Local variables
    double Rainfall_energy = 0.0;
    double intensity = 1.0;
    double common_cropcover = 0.0;
    double common_groundcover = 0.0;
    double transportfactor = 1.0;
    double mobilisedsed = 0.0;

    // Local parameters
    const double trans1 = 4.0;
    const double trans2 = 1.3;

    // Initialize output variable
    erodedsed = 0.0;

    // No erosion if cohesion or erodibility is zero
    if (cohesion == 0.0 || erodibility == 0.0) return;

    // Get current crop and ground cover
    get_current_crop_and_ground_cover(i, j, common_cropcover, common_groundcover);

    // Check for snow limiting erosion
    if (snow > 0.0) intensity = 0.0;

    // Calculate mobilised sediment
    if (prec > 0.0) {
        if (intensity > 0.0) {
            if (prec > 5.0) {
                // Calculate rainfall energy
                Rainfall_energy = 8.95 + 8.44 * log10(prec * (0.257 + sin(2 * M_PI * ((pdayno - 70.0) / 365.0)) * 0.09) * 2.0);
            } else {
                Rainfall_energy = 0.0;
            }
            Rainfall_energy *= prec; // J/m2
            mobilisedsed = Rainfall_energy * (1.0 - common_cropcover) * erodibility; // g/m2
        }
    }

    if (surfacerunoff > 0.0) {
        mobilisedsed += pow((surfacerunoff * 365.0), sreroexp) * (1.0 - common_groundcover) * 
                        (1.0 / (0.5 * cohesion)) * sin(basin[i].slope / 100.0) / 365.0; // g/m2
    }

    // Calculate transport capacity
    if (flow > 0.0) {
        transportfactor = std::min(1.0, pow(flow / trans1, trans2));
    } else {
        transportfactor = 1.0;
    }

    // Calculate eroded sediment
    erodedsed = 1000.0 * mobilisedsed * transportfactor; // kg/km2
}
