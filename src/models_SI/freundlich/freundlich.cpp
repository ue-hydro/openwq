
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
      const double Limit = 0.00001; // Threshold for breaking iterations

    // Local variables
    double totalP, PPequi_conc;
    double conc_sol, adsdes;
    double x0, xn, xn_1, fxn, fprimxn, dx;
    double coeff;
    double nfrloc;
    double help;

    if (Vol == 0) return;

    if (Kfr == 0 || Nfr == 0 || Kadsdes == 0) {
        std::cerr << "ERROR: Values for freundlich parameters missing" << std::endl;
        exit(1);
    }

    nfrloc = Nfr;
    totalP = poolPP + SRP_Conc * Vol; // Total amount of P (kg/km2)
    if (totalP == 0.0) return;

    conc_sol = (poolPP / bulkdensity) / LayerThick; // mg P / kg soil
    if (conc_sol <= 0.0) {
        nfrloc = 1.0;
        if (conductwarning) {
            std::cerr << "Warning: soil partP <= 0. Freundlich will give error, take shortcut." << std::endl;
        }
    }

    coeff = Kfr * bulkdensity * LayerThick;

    if (nfrloc == 1.0) {
        xn_1 = totalP / (Vol + coeff);
        PPequi_conc = Kfr * xn_1;
    } else {
        // Newton-Raphson method to calculate equilibrium concentration
        x0 = exp((log(conc_sol) - log(Kfr)) / Nfr); // Initial guess of equilibrium liquid concentration
        fxn = x0 * Vol + coeff * pow(x0, Nfr) - totalP;
        xn = x0;
        xn_1 = xn;
        int i = 0;

        while (fabs(fxn) > Limit && i < 20) {
            fxn = xn * Vol + coeff * pow(xn, Nfr) - totalP;
            fprimxn = Vol + Nfr * coeff * pow(xn, Nfr - 1);
            dx = fxn / fprimxn;
            if (fabs(dx) < 0.000001 * xn) break;
            xn_1 = xn - dx;
            if (xn_1 <= 0.0) {
                xn_1 = 1.0E-10;
            }
            xn = xn_1;
            ++i;
        }
        PPequi_conc = Kfr * pow(xn_1, Nfr);
    }

    // Calculate new pool and concentration, depends on the equilibrium concentration
    if (fabs(PPequi_conc - conc_sol) > 1.0E-6) {
        adsdes = (PPequi_conc - conc_sol) * (1 - exp(-Kadsdes)); // Kinetic adsorption/desorption
        help = adsdes * bulkdensity * LayerThick;
        if (-help > poolPP || SRP_Conc < (help / Vol)) {
            if (-help > poolPP) help = -poolPP;
            if (SRP_Conc < (help / Vol)) help = SRP_Conc * Vol;
            if (conductwarning) {
                std::cerr << "Warning: freundlich flow adjusted, was larger than pool" << std::endl;
            }
        }
        poolPP += help;       // New Pool PP
        SRP_Conc -= (help / Vol); // New liquid concentration
    }

    // Safety check for negative SRP_Conc
    if (SRP_Conc < 0.0) {
        std::cerr << "ERROR: SRP_Conc in freundlich negative! " << SRP_Conc << std::endl;
        std::cerr << iin << " " << jin << " " << lin << std::endl;
        exit(1);
    }

    */

}
