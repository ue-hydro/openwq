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

// Calculate eroded particles from soil by HBV-sed based model (based on HYPE)

// Constants and parameters for the model (these need to be defined in your program)
struct Basin {
    double eroindex; // Erosion index for the basin
    double slope;    // Slope percentage for the basin
};

struct Time {
    struct {
        int month; // Current month (1-12)
    } date;
};

// External arrays or parameters (placeholders, implement as needed)
double monthpar(int param, int month) {
    // Stub for monthly erosion factor parameter lookup
    // Replace with your data structure implementation
    return 0.0; // Example: default value
}

// Function to calculate eroded particles using the HBV-sed model
void hype_hbvsed(int i, double prec, double lusepar, 
                              double soilpar, double slopepar, 
                              double precexppar, double eroindexpar, double snow, double frostdepth, 
                              double& erodedsed, const Basin* basin, const Time& current_time) {

    // USE HYPEVARIABLES, ONLY : m_erodmon
    // USE MODVAR, ONLY : basin,monthpar,current_time

    //  !Argument declarations
    //  INTEGER, INTENT(IN) :: i             !<index of current subbasin
    //  REAL, INTENT(IN)    :: prec          !<precipitation (rainfall only)    
    //  REAL, INTENT(IN)    :: lusepar       !<soil erosion factor (land use dependence)
    //  REAL, INTENT(IN)    :: soilpar       !<soil erosion factor (soil dependence)
    //  REAL, INTENT(IN)    :: slopepar      !<slope erosion factor (exponent) 
    //  REAL, INTENT(IN)    :: precexppar    !<erosion precipitation dependence factor (exponent)
    //  REAL, INTENT(IN)    :: eroindexpar   !<model parameter for scaling of erosion index
    //  REAL, INTENT(IN)    :: snow          !<snow water (mm)
    // REAL, INTENT(IN)    :: frostdepth    !<frost depth (cm)
    // REAL, INTENT(OUT)   :: erodedsed     !<eroded (transported) sediment (kg/km2)

    // Local variables
    double erosionindex = 0.0;
    double mobilisedsed = 0.0;
    double erodmonth = 0.0;

    // Snow and soil frost turn off particle mobilization
    if (snow > 0.0 || (frostdepth < 0.0 && frostdepth > -9999.0)) {
        erodedsed = 0.0;
        return;
    }

    // Calculate mobilised sediments from erosion index, slope, and rainfall
    erosionindex = basin[i].eroindex / eroindexpar;
    mobilisedsed = std::pow((basin[i].slope / 5.0), slopepar) * lusepar 
                   * soilpar * erosionindex * 
                   std::pow(prec, precexppar); // tonnes/km2 = g/m2

    // Eroded sediment calculated from mobilised sediment with a monthly factor
    erodmonth = 1.0 + monthpar(0, current_time.date.month); // Replace `0` with the appropriate `m_erodmon` parameter
    erodedsed = 1000.0 * mobilisedsed * erodmonth; // kg/km2
}
