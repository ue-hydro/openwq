

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

void OpenWQ_TD_model::EWF_flux_driver_run(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig, 
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output,
    std::string source_EWF_name,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r){ // water flux (m3/s)
    
    // Dummy variable
    double ewf_conc_chemi;  // concentraton of chemical at EWF at current time step
    unsigned int iewf;      // index of EWF
    unsigned int ichem_mob; // iterative dummy var for mobile chemical species

    // Find specific EWF index in EWF names
    iewf = OpenWQ_utils.FindStrIndexInVectStr(
        OpenWQ_hostModelconfig.get_HydroExtFlux_names(),
        source_EWF_name);

    unsigned int numspec = OpenWQ_wqconfig.CH->NativeFlex->mobile_species.size();

    // Loop over mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        // mobile chemical species index
        ichem_mob = OpenWQ_wqconfig.CH->NativeFlex->mobile_species[chemi];

        // Get appropriate chemi concentraton for this ewf and domain ix, iy, iz
        ewf_conc_chemi = (*OpenWQ_vars.ewf_conc)
            (iewf)
            (ichem_mob)
            (ix_r,iy_r,iz_r);

        // Advection and dispersion
        EWF_flux_IN(
            OpenWQ_vars, OpenWQ_wqconfig,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r,                      // external water flux
            ichem_mob,                          // chemical id (from BGQ file and used in state-varibble data structures)
            ewf_conc_chemi);                // concentration taken from EWF json file

    }

}


// Add actual load from precipitation

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