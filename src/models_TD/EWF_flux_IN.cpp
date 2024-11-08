

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


#include "models_TD/headerfile_TD.hpp"

/* #################################################
// Mass transport
// Only Advection
// Special case (in-fluxes from external water flux sources)
################################################# */
void OpenWQ_TD_model::EWF_flux_IN(
    OpenWQ_vars& OpenWQ_vars, 
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r, // water flux (m3/s)
    int chemi,              // chemical id
    double ewf_conc){       // concentration of chemical

    // Local variables
    double chemass_flux_adv;

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f || ewf_conc == 0.0f){return;}

    // Chemical mass flux between source and recipient (Advection)
    chemass_flux_adv = 
        wflux_s2r
        * ewf_conc;
    
    //##########################################
    // Set derivative for source and recipient 
    
    // Add Chemical mass flux to RECIPIENT 
    (*OpenWQ_vars.d_chemass_ewf)(recipient)(chemi)(ix_r,iy_r,iz_r) 
        += chemass_flux_adv;

}