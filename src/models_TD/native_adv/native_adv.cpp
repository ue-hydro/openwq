

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

    // OPTIMIZED: pre-compute concentration factor once
    const double conc_factor = wflux_s2r / wmass_source;

    // OPTIMIZED: pre-fetch field references to avoid repeated pointer dereferencing
    auto& chemass_source = (*OpenWQ_vars.chemass)(source);
    auto& d_transp_source = (*OpenWQ_vars.d_chemass_dt_transp)(source);

    // Loop for mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        const unsigned int ichem_mob = mobile_species[chemi];

        // Chemical mass flux between source and recipient (Advection)
        const double chemass_flux_adv = conc_factor * chemass_source(ichem_mob)(ix_s,iy_s,iz_s);

        // Remove Chemical mass flux from SOURCE
        d_transp_source(ichem_mob)(ix_s,iy_s,iz_s) -= chemass_flux_adv;

        // Add Chemical mass flux to RECIPIENT
        // if recipient == -1, then it's an OUT-flux (loss from system)
        if (recipient == -1) continue;
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
            += chemass_flux_adv;
    }

}
