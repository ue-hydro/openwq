

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


#include "models_TD/headerfile_TD.hpp"


/* #################################################
// Mass transport
// Only Advection
// General case (flux exchanges within the model domain)
// OPTIMIZED: cached flags, pre-fetched refs, pre-computed ratio
################################################# */
void OpenWQ_TD_model::Adv(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r,
    double wmass_source){

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    // OPTIMIZED: use cached values instead of string comparisons
    const unsigned int numspec = OpenWQ_wqconfig.cached_num_mobile_species;
    const std::vector<unsigned int>& mobile_species = *OpenWQ_wqconfig.cached_mobile_species_ptr;

    // Advected fraction of the source mass: the share of the source water that
    // leaves with this flux in the step, conc_factor = min(Q*dt/V, 1), with V
    // the water available in the source cell during the step (as passed by the
    // host). This is the explicit upwind (donor-cell) scheme.
    //
    // It is the fraction that is consistent with the rest of the scheme: the
    // flux is taken from the START-OF-STEP mass only (see below), so the mass
    // that enters the cell in a step stays there until the next one. At steady
    // state the stored mass is then (inflow mass)/conc_factor, and the cell
    // concentration equals the inflow concentration only if conc_factor is
    // Q*dt/V. The exponential form 1 - exp(-Q*dt/V) used before removes less
    // than that, and the cell concentration settled at
    // (Q*dt/V)/(1 - exp(-Q*dt/V)) times the inflow concentration: +5% for
    // Q*dt/V = 0.1 and +58% for Q*dt/V = 1 (river reaches and other
    // through-flow cells). The sediment transport (models_TS) already uses the
    // same linear fraction.
    //
    // No oscillations: the fraction is capped at 1 (a cell cannot export more
    // than it holds) and only start-of-step mass is exported, so the result
    // does not depend on the order in which the host processes the cells.
    const double conc_factor = std::fmin(wflux_s2r / wmass_source, 1.0);

    // OPTIMIZED: pre-fetch field references to avoid repeated pointer dereferencing
    auto& chemass_source = (*OpenWQ_vars.chemass)(source);
    auto& d_transp_source = (*OpenWQ_vars.d_chemass_dt_transp)(source);

    // Loop for mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        const unsigned int ichem_mob = mobile_species[chemi];

        // Use only the start-of-step source mass (chemass) for the advective
        // flux. The previous formulation included d_transp_source as well, so
        // that mass arriving from upstream in the same timestep could be
        // routed further downstream in one step -- but that depends on the
        // host calling reaches in strict topological order. Under MPI/OpenMP
        // (or any routing scheme that processes reaches out of order), the
        // dependence flips between timesteps and produces large alternating
        // oscillations. Using chemass only gives the standard explicit upwind
        // scheme: order-independent, mass-conservative, stable for any CFL.
        // Mass advances at most one reach per timestep.
        double current_mass = std::fmax(
            chemass_source(ichem_mob)(ix_s,iy_s,iz_s), 0.0);

        // Chemical mass flux between source and recipient (Advection).
        // Flux magnitude is set by the start-of-step source mass (order-
        // independent). The safety cap below uses LIVE balance (chemass +
        // d_transp) so that a reach with multiple downstream connections
        // cannot have N successive outflow calls each extract current_mass
        // -- without this, the end-of-step non-negativity clamp would
        // silently create mass at the recipients.
        double chemass_flux_adv = conc_factor * current_mass;
        const double src_live = std::fmax(
            chemass_source(ichem_mob)(ix_s,iy_s,iz_s)
            + d_transp_source(ichem_mob)(ix_s,iy_s,iz_s), 0.0);
        chemass_flux_adv = std::fmin(chemass_flux_adv, src_live);

        // Remove Chemical mass flux from SOURCE
        d_transp_source(ichem_mob)(ix_s,iy_s,iz_s) -= chemass_flux_adv;

        // Add Chemical mass flux to RECIPIENT
        // if recipient == -1, then it's an OUT-flux (loss from system)
        if (recipient == -1) {
            // Track mass leaving the system for mass balance
            if (OpenWQ_vars.mass_balance.initialized &&
                ichem_mob < OpenWQ_vars.mass_balance.num_species) {
                OpenWQ_vars.mass_balance.cumulative_out_flux[ichem_mob] += chemass_flux_adv;
            }
            continue;
        }
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
            += chemass_flux_adv;
    }

}
