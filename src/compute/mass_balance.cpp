// Copyright 2026, Diogo Costa (diogo.costa@uevora.pt)
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "compute/headerfile_compute.hpp"
#include "output/OpenWQ_output.hpp"
#include <cmath>
#include <iomanip>
#include <sstream>

/* #################################################
// Mass Balance Tracking Implementation
//
// This module tracks mass conservation throughout the simulation:
// - Total dissolved mass (chemass)
// - Total sorbed mass (sorbed_mass)
// - Cumulative sources (IC, SS, EWF)
// - Cumulative sinks (transport out, negative corrections)
// - Mass balance error computation
//
// The mass balance equation:
//   initial_mass + cumulative_sources - cumulative_sinks = current_mass
//
// Any deviation indicates mass leakage.
################################################# */

/* #################################################
// Initialize mass balance tracking
// Call once at simulation start, after IC is applied
################################################# */
void OpenWQ_compute::MassBalance_Initialize(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_output& OpenWQ_output) {

    // Determine number of chemical species
    unsigned int num_chem = 0;
    if (OpenWQ_wqconfig.is_native_bgc_flex) {
        num_chem = OpenWQ_wqconfig.CH_model->NativeFlex->num_chem;
    } else {
        num_chem = OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;
    }

    if (num_chem == 0) return;

    // Initialize mass balance data structure
    OpenWQ_vars.mass_balance.initialize(num_chem);

    // Compute initial total mass
    MassBalance_ComputeTotals(OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_vars);

    // Store as initial mass
    for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
        OpenWQ_vars.mass_balance.initial_mass[chemi] =
            OpenWQ_vars.mass_balance.total_mass[chemi];
    }

    // Log initialization
    std::string msg_string = "<OpenWQ> Mass Balance Tracking initialized for "
        + std::to_string(num_chem) + " species";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Log initial mass for each species
    for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
        double init_mass = OpenWQ_vars.mass_balance.initial_mass[chemi];
        if (init_mass > 0) {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(4);
            oss << "<OpenWQ> Species " << chemi << " initial mass: " << init_mass << " g";
            msg_string = oss.str();
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }
}

/* #################################################
// Compute current total mass across all compartments
################################################# */
void OpenWQ_compute::MassBalance_ComputeTotals(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars) {

    if (!OpenWQ_vars.mass_balance.initialized) return;

    const unsigned int num_chem = OpenWQ_vars.mass_balance.num_species;
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Reset totals
    std::fill(OpenWQ_vars.mass_balance.total_dissolved_mass.begin(),
              OpenWQ_vars.mass_balance.total_dissolved_mass.end(), 0.0);
    std::fill(OpenWQ_vars.mass_balance.total_sorbed_mass.begin(),
              OpenWQ_vars.mass_balance.total_sorbed_mass.end(), 0.0);

    // Sum dissolved mass (chemass)
    for (unsigned int icmp = 0; icmp < num_comps; icmp++) {
        const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

        for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
            double species_sum = 0.0;

            // Sum all cells for this species in this compartment
            const auto& chemass_cube = (*OpenWQ_vars.chemass)(icmp)(chemi);
            for (unsigned int ix = 0; ix < nx; ix++) {
                for (unsigned int iy = 0; iy < ny; iy++) {
                    for (unsigned int iz = 0; iz < nz; iz++) {
                        species_sum += chemass_cube(ix, iy, iz);
                    }
                }
            }

            OpenWQ_vars.mass_balance.total_dissolved_mass[chemi] += species_sum;
        }
    }

    // Sum sorbed mass (if sorbed_mass is allocated and populated)
    if (OpenWQ_vars.sorbed_mass) {
        for (unsigned int icmp = 0; icmp < num_comps; icmp++) {
            // Check if this compartment's sorbed_mass is allocated
            if ((*OpenWQ_vars.sorbed_mass)(icmp).n_elem == 0) continue;

            const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
            const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
            const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

            for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
                // Check if this species has sorbed_mass tracking
                if ((*OpenWQ_vars.sorbed_mass)(icmp).n_elem <= chemi) continue;
                if ((*OpenWQ_vars.sorbed_mass)(icmp)(chemi).n_elem == 0) continue;

                double species_sum = 0.0;
                const auto& sorbed_cube = (*OpenWQ_vars.sorbed_mass)(icmp)(chemi);

                for (unsigned int ix = 0; ix < nx; ix++) {
                    for (unsigned int iy = 0; iy < ny; iy++) {
                        for (unsigned int iz = 0; iz < nz; iz++) {
                            species_sum += sorbed_cube(ix, iy, iz);
                        }
                    }
                }

                OpenWQ_vars.mass_balance.total_sorbed_mass[chemi] += species_sum;
            }
        }
    }

    // Compute total mass = dissolved + sorbed
    for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
        OpenWQ_vars.mass_balance.total_mass[chemi] =
            OpenWQ_vars.mass_balance.total_dissolved_mass[chemi] +
            OpenWQ_vars.mass_balance.total_sorbed_mass[chemi];
    }
}

/* #################################################
// Update cumulative fluxes from derivative arrays
################################################# */
void OpenWQ_compute::MassBalance_UpdateFluxes(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars) {

    if (!OpenWQ_vars.mass_balance.initialized) return;

    const unsigned int num_chem = OpenWQ_vars.mass_balance.num_species;
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Sum sources (SS + EWF with positive values)
    // Sum sinks (SS + EWF with negative values, transport out)
    for (unsigned int icmp = 0; icmp < num_comps; icmp++) {
        const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

        for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
            const auto& ss_cube = (*OpenWQ_vars.d_chemass_ss)(icmp)(chemi);
            const auto& ewf_cube = (*OpenWQ_vars.d_chemass_ewf)(icmp)(chemi);

            for (unsigned int ix = 0; ix < nx; ix++) {
                for (unsigned int iy = 0; iy < ny; iy++) {
                    for (unsigned int iz = 0; iz < nz; iz++) {
                        // SS flux
                        double ss_flux = ss_cube(ix, iy, iz);
                        if (ss_flux > 0) {
                            OpenWQ_vars.mass_balance.cumulative_sources[chemi] += ss_flux;
                        } else {
                            OpenWQ_vars.mass_balance.cumulative_sinks[chemi] -= ss_flux;
                        }

                        // EWF flux (external water flux - always a source)
                        double ewf_flux = ewf_cube(ix, iy, iz);
                        if (ewf_flux > 0) {
                            OpenWQ_vars.mass_balance.cumulative_sources[chemi] += ewf_flux;
                        }
                    }
                }
            }
        }
    }
}

/* #################################################
// Compute mass balance error
################################################# */
void OpenWQ_compute::MassBalance_ComputeError(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars) {

    if (!OpenWQ_vars.mass_balance.initialized) return;

    const unsigned int num_chem = OpenWQ_vars.mass_balance.num_species;

    // First, recompute current totals
    MassBalance_ComputeTotals(OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_vars);

    for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
        // Expected mass = initial + sources - sinks - out_flux - negative_correction
        double expected_mass =
            OpenWQ_vars.mass_balance.initial_mass[chemi]
            + OpenWQ_vars.mass_balance.cumulative_sources[chemi]
            - OpenWQ_vars.mass_balance.cumulative_sinks[chemi]
            - OpenWQ_vars.mass_balance.cumulative_out_flux[chemi]
            - OpenWQ_vars.mass_balance.cumulative_negative_correction[chemi];

        // Actual mass
        double actual_mass = OpenWQ_vars.mass_balance.total_mass[chemi];

        // Error = expected - actual
        double error = expected_mass - actual_mass;
        OpenWQ_vars.mass_balance.mass_balance_error[chemi] = error;

        // Error as percentage of initial + sources
        double reference = OpenWQ_vars.mass_balance.initial_mass[chemi]
                          + OpenWQ_vars.mass_balance.cumulative_sources[chemi];
        if (reference > 1e-10) {
            OpenWQ_vars.mass_balance.mass_balance_error_pct[chemi] =
                (error / reference) * 100.0;
        } else if (actual_mass > 1e-10) {
            OpenWQ_vars.mass_balance.mass_balance_error_pct[chemi] =
                (error / actual_mass) * 100.0;
        } else {
            OpenWQ_vars.mass_balance.mass_balance_error_pct[chemi] = 0.0;
        }
    }
}

/* #################################################
// Print mass balance report to console/log
################################################# */
void OpenWQ_compute::MassBalance_Report(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_output& OpenWQ_output) {

    if (!OpenWQ_vars.mass_balance.initialized) return;

    // First compute error
    MassBalance_ComputeError(OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_vars);

    const unsigned int num_chem = OpenWQ_vars.mass_balance.num_species;

    std::string msg_string;

    msg_string = "<OpenWQ> ============== MASS BALANCE REPORT ==============";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    for (unsigned int chemi = 0; chemi < num_chem; chemi++) {
        // Skip species with no mass at all
        double total_involved = OpenWQ_vars.mass_balance.initial_mass[chemi]
                               + OpenWQ_vars.mass_balance.cumulative_sources[chemi];
        if (total_involved < 1e-10 && OpenWQ_vars.mass_balance.total_mass[chemi] < 1e-10) {
            continue;
        }

        std::ostringstream oss;
        oss << std::scientific << std::setprecision(4);

        msg_string = "<OpenWQ> --- Species " + std::to_string(chemi) + " ---";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Initial mass:     " << OpenWQ_vars.mass_balance.initial_mass[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Cum. sources:     +" << OpenWQ_vars.mass_balance.cumulative_sources[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Cum. sinks:       -" << OpenWQ_vars.mass_balance.cumulative_sinks[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Cum. out-flux:    -" << OpenWQ_vars.mass_balance.cumulative_out_flux[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Cum. neg. corr:   -" << OpenWQ_vars.mass_balance.cumulative_negative_correction[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Current dissolved: " << OpenWQ_vars.mass_balance.total_dissolved_mass[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Current sorbed:    " << OpenWQ_vars.mass_balance.total_sorbed_mass[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << "<OpenWQ>   Current total:     " << OpenWQ_vars.mass_balance.total_mass[chemi] << " g";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        oss.str(""); oss.clear();
        oss << std::fixed << std::setprecision(6);
        oss << "<OpenWQ>   MASS BALANCE ERROR: " << OpenWQ_vars.mass_balance.mass_balance_error[chemi] << " g";
        oss << " (" << std::setprecision(4) << OpenWQ_vars.mass_balance.mass_balance_error_pct[chemi] << "%)";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, oss.str(), true, true);

        // Warning if error is significant
        if (std::fabs(OpenWQ_vars.mass_balance.mass_balance_error_pct[chemi]) > 1.0) {
            msg_string = "<OpenWQ>   WARNING: Mass balance error > 1%!";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }

    msg_string = "<OpenWQ> =====================================================";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
}

/* #################################################
// Track mass lost to negative clamping
// Called from solver when mass goes negative
################################################# */
void OpenWQ_compute::MassBalance_TrackNegativeCorrection(
    OpenWQ_vars& OpenWQ_vars,
    unsigned int chemi,
    double mass_lost) {

    if (!OpenWQ_vars.mass_balance.initialized) return;
    if (chemi >= OpenWQ_vars.mass_balance.num_species) return;

    // Add to cumulative negative correction (thread-safe would need atomic)
    OpenWQ_vars.mass_balance.cumulative_negative_correction[chemi] += mass_lost;
}
