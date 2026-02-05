
// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
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

/* #################################################
// Freundlich Sorption Isotherm
// q = Kfr * C^(1/Nfr)
//
// Uses Newton-Raphson to find equilibrium dissolved concentration C_eq
// given total mass (dissolved + sorbed) conservation:
//     C_eq * V + Kfr * C_eq^Nfr * rho * L = total_mass
//
// Then applies kinetic adsorption/desorption:
//     adsdes_flux = (q_eq - q_current) * (1 - exp(-Kadsdes * dt)) * rho * L
//
// The flux is added to d_chemass_dt_chem (removes from dissolved, adds to sorbed phase).
// In the current simplified implementation, we track only the dissolved phase
// (chemass) and the sorption acts as a sink/source on the dissolved concentration.
################################################# */

void OpenWQ_SI_model::freundlich(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){

    // Newton-Raphson convergence parameters
    const double NR_Limit = 1.0e-5;     // Convergence threshold
    const int NR_MaxIter = 20;           // Maximum iterations
    const double NR_dx_rel_tol = 1.0e-6; // Relative step tolerance

    // Get Freundlich parameters
    auto* FR = OpenWQ_wqconfig.SI_model->FREUNDLICH;
    const unsigned int num_si_species = FR->num_species;

    if (num_si_species == 0) return;  // Nothing to do

    const double bulk_density = FR->bulk_density;    // kg/m3
    const double layer_thick = FR->layer_thickness;  // m
    const double rhoL = bulk_density * layer_thick;  // kg/m2

    // Get timestep for kinetic rate
    const double dt = OpenWQ_hostModelconfig.get_time_step();  // seconds
    if (dt <= 0.0) return;

    // Get number of compartments
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Water volume minimum limit (to avoid division by zero)
    const double watervol_minlim = OpenWQ_hostModelconfig.get_watervol_minlim();

    // Loop over each SI species
    for (unsigned int si = 0; si < num_si_species; si++) {

        const unsigned int chemi = FR->species_index[si];
        const double Kfr = FR->Kfr[si];
        const double Nfr = FR->Nfr[si];
        const double Kadsdes = FR->Kadsdes[si];

        // Skip if parameters are invalid
        if (Kfr <= 0.0 || Nfr <= 0.0 || Kadsdes <= 0.0) continue;

        // Pre-compute kinetic factor: (1 - exp(-Kadsdes * dt))
        const double kinetic_factor = 1.0 - std::exp(-Kadsdes * dt);

        // Coefficient for mass balance: Kfr * rho * L
        const double coeff = Kfr * rhoL;

        // Loop over all compartments and cells
        for (unsigned int icmp = 0; icmp < num_comps; icmp++) {

            const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
            const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
            const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

            // References to state and derivative cubes
            auto& chemass = (*OpenWQ_vars.chemass)(icmp)(chemi);
            auto& d_chemass = (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(chemi);

            for (unsigned int ix = 0; ix < nx; ix++) {
                for (unsigned int iy = 0; iy < ny; iy++) {
                    for (unsigned int iz = 0; iz < nz; iz++) {

                        // Get water volume [m3] in this cell
                        const double Vol = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
                            icmp, ix, iy, iz);

                        // Skip dry cells
                        if (Vol <= watervol_minlim) continue;

                        // Current dissolved mass [g] in this cell
                        const double mass_dissolved = chemass(ix, iy, iz);
                        if (mass_dissolved <= 0.0) continue;

                        // Current dissolved concentration [g/m3 = mg/L]
                        const double C_current = mass_dissolved / Vol;

                        // Current sorbed concentration [mg/kg_soil]
                        // q_current = Kfr * C_current^Nfr  (assuming at current quasi-equilibrium)
                        // Note: We don't track sorbed mass separately; we compute
                        // the equilibrium target and apply kinetic correction

                        // Total mass in this cell (dissolved + sorbed)
                        // total_mass = C * Vol + q * rho * L
                        // where q = Kfr * C^Nfr (current sorbed conc estimate)
                        const double q_current = Kfr * std::pow(C_current, Nfr);
                        const double total_mass = C_current * Vol + q_current * rhoL;

                        if (total_mass <= 0.0) continue;

                        // ############################
                        // Newton-Raphson to find equilibrium dissolved concentration C_eq
                        // f(C) = C * Vol + Kfr * C^Nfr * rho * L - total_mass = 0
                        // f'(C) = Vol + Nfr * Kfr * C^(Nfr-1) * rho * L
                        // ############################

                        double C_eq;
                        double nfrloc = Nfr;

                        // Handle edge case where sorbed conc is zero or negative
                        if (q_current <= 0.0) {
                            nfrloc = 1.0;  // Linear shortcut
                        }

                        if (std::fabs(nfrloc - 1.0) < 1.0e-10) {
                            // Linear case: C_eq * Vol + Kfr * C_eq * rho * L = total_mass
                            C_eq = total_mass / (Vol + coeff);
                        } else {
                            // Newton-Raphson for nonlinear Freundlich
                            // Initial guess from inverse Freundlich
                            double x0;
                            if (C_current > 1.0e-10) {
                                x0 = C_current;  // Use current concentration as initial guess
                            } else {
                                x0 = total_mass / Vol;  // Fallback: assume all dissolved
                            }

                            double xn = x0;
                            double fxn, fprimxn, dx;

                            for (int iter = 0; iter < NR_MaxIter; iter++) {
                                fxn = xn * Vol + coeff * std::pow(xn, Nfr) - total_mass;

                                if (std::fabs(fxn) < NR_Limit) break;

                                fprimxn = Vol + Nfr * coeff * std::pow(xn, Nfr - 1.0);
                                dx = fxn / fprimxn;

                                if (std::fabs(dx) < NR_dx_rel_tol * std::fabs(xn)) break;

                                xn -= dx;

                                // Ensure positive concentration
                                if (xn <= 0.0) xn = 1.0e-10;
                            }

                            C_eq = xn;
                        }

                        // Equilibrium sorbed concentration
                        const double q_eq = Kfr * std::pow(C_eq, Nfr);

                        // ############################
                        // Kinetic adsorption/desorption
                        // Mass flux from dissolved to sorbed phase
                        // positive = adsorption (remove from dissolved)
                        // negative = desorption (add to dissolved)
                        // ############################

                        const double delta_q = (q_eq - q_current) * kinetic_factor;
                        double flux_mass = delta_q * rhoL;  // [g]

                        // Safety: flux cannot remove more mass than available
                        if (flux_mass > 0.0 && flux_mass > mass_dissolved) {
                            flux_mass = mass_dissolved;
                        }
                        // Safety: desorption cannot create negative sorbed pool
                        // (not tracked explicitly, but limit to reasonable bounds)

                        // Apply as a sink on dissolved concentration
                        // (positive flux_mass = adsorption = decrease dissolved)
                        d_chemass(ix, iy, iz) -= flux_mass;

                    } // iz
                } // iy
            } // ix
        } // icmp
    } // si species
}
