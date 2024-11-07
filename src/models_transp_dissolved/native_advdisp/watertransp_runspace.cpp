

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


#include "models_transp_dissolved/headerfile_td.hpp"


/* #################################################
// Mass transport
// Only Advection
// General case (flux exchanges within the model domain)
################################################# */
void OpenWQ_watertransp::Adv(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    double wmass_source){

    double chemass_flux_adv;
    unsigned int ichem_mob;

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    // CHANGE THE LINE BELOW: see my notes -> there should be no icmp because all compartments should have the same number of mobile species
    unsigned int numspec = OpenWQ_wqconfig.BGC_general_mobile_species.size();

    // Loop for mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        // mobile chemical species index
        ichem_mob = OpenWQ_wqconfig.BGC_general_mobile_species[chemi];

        // Chemical mass flux between source and recipient (Advection)
        chemass_flux_adv = 
            wflux_s2r
            * (*OpenWQ_vars.chemass)(source)(ichem_mob)(ix_s,iy_s,iz_s) // concentration calculation
             / wmass_source;
        
        //##########################################
        // Set derivative for source and recipient 
        
        // Remove Chemical mass flux from SOURCE 
        (*OpenWQ_vars.d_chemass_dt_transp)(source)(ichem_mob)(ix_s,iy_s,iz_s) 
            -= chemass_flux_adv;

        // Add Chemical mass flux to RECIPIENT 
        // if recipient == -1, then  it's an OUT-flux (loss from system)
        if (recipient == -1) continue;
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r) 
            += chemass_flux_adv;
    }
                
}

/* #################################################
// Mass transport
// Only Advection
// Special case (in-fluxes from external water flux sources)
################################################# */
void OpenWQ_watertransp::Adv_IN(
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

/* #################################################
// Mass transport
// Advection & Dispersion
################################################# */
void OpenWQ_watertransp::AdvDisp(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    double wmass_source){

    double chemass_flux_adv;
    unsigned int ichem_mob;

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    // CHANGE THE LINE BELOW: see my notes -> there should be no icmp because all compartments should have the same number of mobile species
    unsigned int numspec = OpenWQ_wqconfig.BGC_general_mobile_species.size();

    // Loop for mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        // mobile chemical species index
        ichem_mob = OpenWQ_wqconfig.BGC_general_mobile_species[chemi];

        // Chemical mass flux between source and recipient (Advection)
        chemass_flux_adv = 
            wflux_s2r
            * (*OpenWQ_vars.chemass)(source)(ichem_mob)(ix_s,iy_s,iz_s) // concentration calculation
             / wmass_source;
        
        // dcdy(xi,yi,t)=(c(xi,yi,t-1)-c(xi,yi-1,t-1))/delty;
        // dc2dy2(xi,yi,t)=((c(xi,yi+1,t-1)-c(xi,yi,t-1))/delty-(c(xi,yi,t-1)-c(xi,yi-1,t-1))/delty)/(delty);         
        // dc2dx2(xi,yi,t)=((c(xi+1,yi,t-1)-c(xi,yi,t-1))/deltx-(c(xi,yi,t-1)-c(xi-1,yi,t-1))/deltx)/(deltx);
        // dcdt(xi,yi,t)=-U*dcdy(xi,yi,t)+Di(xi,yi,t)*(dc2dx2(xi,yi,t)+dc2dy2(xi,yi,t));

        //##########################################
        // Set derivative for source and recipient 
        
        // Remove Chemical mass flux from SOURCE 
        (*OpenWQ_vars.d_chemass_dt_transp)(source)(ichem_mob)(ix_s,iy_s,iz_s) 
            -= chemass_flux_adv;

        // Add Chemical mass flux to RECIPIENT
        // if recipient == -1, then  it's an OUT-flux (loss from system)
        if (recipient == -1) continue;
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r) 
            += chemass_flux_adv;
    }
                
}