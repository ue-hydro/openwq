
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
// Langmuir Sorption Isotherm
// q = (qmax * KL * C) / (1 + KL * C)
//
// The Langmuir isotherm has an analytical solution for equilibrium:
// Given total mass conservation:
//     total_mass = C * Vol + q * rho * L
//     total_mass = C * Vol + (qmax * KL * C) / (1 + KL * C) * rho * L
//
// This leads to a quadratic in C:
//     KL * Vol * C^2 + (Vol + qmax * KL * rhoL - KL * total_mass) * C - total_mass = 0
//
// Then applies kinetic adsorption/desorption:
//     flux = (q_eq - q_current) * (1 - exp(-Kadsdes * dt)) * rho * L
//
// The flux is applied to d_chemass_dt_chem as a sink on dissolved mass.
################################################# */

void OpenWQ_SI_model::langmuir(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){

    // Get Langmuir parameters
    auto* LG = OpenWQ_wqconfig.SI_model->LANGMUIR;
    const unsigned int num_si_species = LG->num_species;

    if (num_si_species == 0) return;  // Nothing to do

    const double bulk_density = LG->bulk_density;    // kg/m3
    const double layer_thick = LG->layer_thickness;  // m
    const double rhoL = bulk_density * layer_thick;  // kg/m2

    // Get timestep for kinetic rate
    const double dt = OpenWQ_hostModelconfig.get_time_step();  // seconds
    if (dt <= 0.0) return;

    // Get number of compartments
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Water volume minimum limit
    const double watervol_minlim = OpenWQ_hostModelconfig.get_watervol_minlim();

    // Loop over each SI species
    for (unsigned int si = 0; si < num_si_species; si++) {

        const unsigned int chemi = LG->species_index[si];
        const double qmax = LG->qmax[si];
        const double KL = LG->KL[si];
        const double Kadsdes = LG->Kadsdes[si];

        // Skip if parameters are invalid
        if (qmax <= 0.0 || KL <= 0.0 || Kadsdes <= 0.0) continue;

        // Pre-compute kinetic factor: (1 - exp(-Kadsdes * dt))
        const double kinetic_factor = 1.0 - std::exp(-Kadsdes * dt);

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

                        // Current sorbed concentration (from Langmuir at current C)
                        // q_current = (qmax * KL * C) / (1 + KL * C)
                        const double q_current = (qmax * KL * C_current) / (1.0 + KL * C_current);

                        // Total mass in this cell (dissolved + sorbed)
                        const double total_mass = C_current * Vol + q_current * rhoL;

                        if (total_mass <= 0.0) continue;

                        // ############################
                        // Analytical solution for equilibrium C_eq via quadratic formula
                        //
                        // From mass conservation with Langmuir:
                        //   C * Vol + (qmax * KL * C / (1 + KL * C)) * rhoL = total_mass
                        //
                        // Multiply through by (1 + KL * C):
                        //   C * Vol * (1 + KL * C) + qmax * KL * C * rhoL = total_mass * (1 + KL * C)
                        //   C * Vol + KL * Vol * C^2 + qmax * KL * rhoL * C = total_mass + KL * total_mass * C
                        //
                        // Rearranging:
                        //   KL * Vol * C^2 + (Vol + qmax * KL * rhoL - KL * total_mass) * C - total_mass = 0
                        //
                        // a = KL * Vol
                        // b = Vol + qmax * KL * rhoL - KL * total_mass
                        // c = -total_mass
                        // ############################

                        const double a = KL * Vol;
                        const double b = Vol + qmax * KL * rhoL - KL * total_mass;
                        const double c = -total_mass;

                        const double discriminant = b * b - 4.0 * a * c;

                        double C_eq;
                        if (discriminant < 0.0 || a <= 0.0) {
                            // Fallback: assume all dissolved
                            C_eq = total_mass / Vol;
                        } else {
                            // Take positive root
                            C_eq = (-b + std::sqrt(discriminant)) / (2.0 * a);
                            if (C_eq < 0.0) {
                                C_eq = (-b - std::sqrt(discriminant)) / (2.0 * a);
                            }
                            if (C_eq < 0.0) {
                                C_eq = total_mass / Vol;  // Fallback
                            }
                        }

                        // Equilibrium sorbed concentration
                        const double q_eq = (qmax * KL * C_eq) / (1.0 + KL * C_eq);

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

                        // Apply as a sink on dissolved concentration
                        // (positive flux_mass = adsorption = decrease dissolved)
                        d_chemass(ix, iy, iz) -= flux_mass;

                    } // iz
                } // iy
            } // ix
        } // icmp
    } // si species
}
