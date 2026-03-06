

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
// Advection & Dispersion (Fickian dispersive flux approximation)
//
// The advective flux is: F_adv = (wflux_s2r / wmass_source) * chemass_source
//   which moves mass proportionally to the water flux (concentration-based).
//
// The dispersive flux uses a Fickian approximation between cell pairs:
//   F_disp = D_eff * (C_source - C_recipient) * wmass_source
// where:
//   D_eff = D_avg / L^2  [1/s]
//   D_avg = (Dx + Dy + Dz) / 3  [m2/s]  (averaged since we don't know connection orientation)
//   L = characteristic_length_m  [m]  (user-provided distance between cell centers)
//   C_source = chemass_source / wmass_source
//   C_recipient = chemass_recipient / wmass_recipient
//
// The dispersive flux is BIDIRECTIONAL: it can move mass in either direction
// depending on the concentration gradient, smoothing concentration differences
// regardless of flow direction. This is physically correct for dispersion.
//
// OPTIMIZED: cached flags, pre-fetched refs, pre-computed D_eff
################################################# */
void OpenWQ_TD_model::AdvDisp(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
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

    // OPTIMIZED: pre-compute advective concentration factor once
    const double conc_factor = wflux_s2r / wmass_source;

    // Get pre-computed effective dispersion rate [1/s]
    const double D_eff = OpenWQ_wqconfig.TD_model->NativeAdvDisp->D_eff;

    // Get recipient water mass from host model for concentration calculation
    // wmass_recipient is needed to compute the recipient concentration: C_r = chemass_r / wmass_r
    double wmass_recipient = 0.0;
    bool has_recipient = (recipient != -1);
    if (has_recipient){
        wmass_recipient = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
            recipient, ix_r, iy_r, iz_r);
    }

    // OPTIMIZED: pre-fetch field references to avoid repeated pointer dereferencing
    auto& chemass_source = (*OpenWQ_vars.chemass)(source);
    auto& d_transp_source = (*OpenWQ_vars.d_chemass_dt_transp)(source);

    // Loop for mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        const unsigned int ichem_mob = mobile_species[chemi];

        // ===========================
        // 1) ADVECTIVE FLUX
        // F_adv = (wflux_s2r / wmass_source) * chemass_source
        // ===========================
        double chemass_flux_adv = conc_factor * chemass_source(ichem_mob)(ix_s,iy_s,iz_s);

        // Cap at available mass (accounting for prior transport in this time step)
        // to prevent over-extraction when CFL > 1 or multiple connections export from same cell
        double available = chemass_source(ichem_mob)(ix_s,iy_s,iz_s)
                         + d_transp_source(ichem_mob)(ix_s,iy_s,iz_s);
        chemass_flux_adv = std::fmin(chemass_flux_adv, std::fmax(available, 0.0));

        // Remove advective flux from SOURCE
        d_transp_source(ichem_mob)(ix_s,iy_s,iz_s) -= chemass_flux_adv;

        // Add advective flux to RECIPIENT (if not an OUT-flux)
        if (!has_recipient){
            // Track mass leaving the system for mass balance
            if (OpenWQ_vars.mass_balance.initialized &&
                ichem_mob < OpenWQ_vars.mass_balance.num_species) {
                OpenWQ_vars.mass_balance.cumulative_out_flux[ichem_mob] += chemass_flux_adv;
            }
            continue;
        }
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
            += chemass_flux_adv;

        // ===========================
        // 2) DISPERSIVE FLUX (Fickian approximation between cell pairs)
        // F_disp = D_eff * (C_source - C_recipient) * wmass_source
        // Positive F_disp means mass moves from source to recipient (source has higher concentration)
        // Negative F_disp means mass moves from recipient to source (recipient has higher concentration)
        // ===========================
        if (D_eff > 0.0 && wmass_recipient > 0.0){

            // Compute concentrations [mass/mass] in source and recipient
            const double C_source = chemass_source(ichem_mob)(ix_s,iy_s,iz_s) / wmass_source;
            const double C_recipient = (*OpenWQ_vars.chemass)(recipient)(ichem_mob)(ix_r,iy_r,iz_r) / wmass_recipient;

            // Dispersive mass flux [mass/time]
            double chemass_flux_disp = D_eff * (C_source - C_recipient) * wmass_source;

            // Cap dispersive flux at available mass to prevent over-extraction
            if (chemass_flux_disp > 0.0) {
                // Source loses mass: cap at remaining available (after advective removal)
                double avail_after_adv = chemass_source(ichem_mob)(ix_s,iy_s,iz_s)
                                       + d_transp_source(ichem_mob)(ix_s,iy_s,iz_s);
                chemass_flux_disp = std::fmin(chemass_flux_disp, std::fmax(avail_after_adv, 0.0));
            } else if (chemass_flux_disp < 0.0) {
                // Recipient loses mass: cap at recipient's remaining available
                double recip_avail = (*OpenWQ_vars.chemass)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
                                   + (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r);
                chemass_flux_disp = std::fmax(chemass_flux_disp, -std::fmax(recip_avail, 0.0));
            }

            // Apply dispersive flux: remove from source, add to recipient
            d_transp_source(ichem_mob)(ix_s,iy_s,iz_s) -= chemass_flux_disp;
            (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
                += chemass_flux_disp;
        }
    }

}
