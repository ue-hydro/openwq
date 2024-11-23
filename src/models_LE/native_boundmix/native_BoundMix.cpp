

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


#include "models_LE/headerfile_LE.hpp"

/* #################################################
// Boundary Mixing due to velocity gradients
// due to turbulence and cross-boarder eddies
################################################# */
void OpenWQ_LE_model::native_BoundMix(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    double wmass_source){
    
    // Local Variables
    double chemass_exchange_upper_comprt;
    double chemass_exchange_lower_comprt;
    unsigned int chemi_mob;             // interactive index of mobile species
    std::vector<unsigned int> xyz_source;
    xyz_source.push_back(ix_s);
    xyz_source.push_back(iy_s);
    xyz_source.push_back(iz_s);         // get x,y,z from source comparment in vector

    // Dummy/iterative variables
    unsigned int input_direction_index;
    int input_upper_compartment_index;
    int input_lower_compartment_index;
    double input_k_val;
    unsigned int index_lower_cell;     // index of lower cell for mass exchange
    std::vector<int> xyz_upper_compartment;

    // Return if flux across compartments
    // This lateral mixing only occurs when fluxes occur when fluxes along the interface of compartments or if leaving the system
    if ((source != recipient) && (recipient != -1))
        return;

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    for (unsigned int entry_i = 0; entry_i < OpenWQ_wqconfig.LE->BoundMix->info_vector.size(); entry_i++){

        xyz_upper_compartment.clear();

        // Get inputs of entry
        input_direction_index = std::get<0>(OpenWQ_wqconfig.LE->BoundMix->info_vector[entry_i]);
        input_upper_compartment_index = std::get<1>(OpenWQ_wqconfig.LE->BoundMix->info_vector[entry_i]);
        input_lower_compartment_index = std::get<2>(OpenWQ_wqconfig.LE->BoundMix->info_vector[entry_i]);
        input_k_val = std::get<3>(OpenWQ_wqconfig.LE->BoundMix->info_vector[entry_i]);

        // Num of cells in upper_compartment
        xyz_upper_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(input_upper_compartment_index));
        xyz_upper_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(input_upper_compartment_index));
        xyz_upper_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(input_upper_compartment_index));

        // Ignore if entry is not applicable to the current source compartment
        if (input_upper_compartment_index != source)
            continue;

        // Get number of cells in input_direction_index
        index_lower_cell = xyz_upper_compartment[input_direction_index] - 1;

        // Ignore if current source cell is not the lower adjacent cell where mixing may occur
        if (xyz_source[input_direction_index] != index_lower_cell)
            continue;


        // Loop over mobile species
        for (unsigned int chemi=0;chemi<OpenWQ_wqconfig.CH->NativeFlex->mobile_species.size();chemi++){
   
            // Get index of mobile species
            chemi_mob = OpenWQ_wqconfig.CH->NativeFlex->mobile_species[chemi];
            
            //##########################################
            // Chemical exxhange between upper and lower compartments
            chemass_exchange_upper_comprt = 
                fmin(
                    input_k_val
                    * wflux_s2r
                    * (*OpenWQ_vars.chemass)(input_upper_compartment_index)(chemi_mob)(ix_s,iy_s,iz_s),
                    (*OpenWQ_vars.chemass)(input_upper_compartment_index)(chemi_mob)(ix_s,iy_s,iz_s));
            
            chemass_exchange_lower_comprt = 
                fmin(
                    input_k_val
                    * wflux_s2r
                    * (*OpenWQ_vars.chemass)(input_lower_compartment_index)(chemi_mob)(ix_s,iy_s,iz_s),
                    (*OpenWQ_vars.chemass)(input_lower_compartment_index)(chemi_mob)(ix_s,iy_s,iz_s));
                
            //##########################################
            // Set derivative for source and recipient 
            
            // Remove Chemical mass flux from SOURCE 
            (*OpenWQ_vars.d_chemass_dt_transp)(input_upper_compartment_index)(chemi_mob)(ix_s,iy_s,iz_s) 
                += (chemass_exchange_lower_comprt - chemass_exchange_upper_comprt);

            // Add Chemical mass flux to RECIPIENT
            // if recipient == -1, then  it's an OUT-flux (loss from system)
            if (recipient == -1) continue;
            (*OpenWQ_vars.d_chemass_dt_transp)(input_lower_compartment_index)(chemi_mob)(ix_r,iy_r,iz_r) 
                += (chemass_exchange_upper_comprt - chemass_exchange_lower_comprt);

        }

    }
}

/* #################################################
// Internal mobilization of immobile pools
// Erosion and weathering
################################################# */
/*
void OpenWQ_TD_model::IntMob(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    double wmass_source){

    // Internal Variables
    double chemass_flux;

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    // Loop for immobile chemical species
    for (unsigned int chemi=0;chemi<OpenWQ_wqconfig.CH->NativeFlex->num_chem;chemi++){

        //##########################################
        // Check if the chem species is mobile
        // Skip if mobile (mobile species are moved by OpenWQ_TD_model::AdvDisp)
        std::vector<unsigned int>::iterator itr 
            = std::find(
                OpenWQ_wqconfig.CH->NativeFlex->mobile_species.begin(), 
                OpenWQ_wqconfig.CH->NativeFlex->mobile_species.end(), 
                chemi);

        if (itr != end(OpenWQ_wqconfig.CH->NativeFlex->mobile_species))
            continue;

        //##########################################
        // Chemical mass flux between source and recipient 
        chemass_flux = 
            OpenWQ_wqconfig.OpenWQ_TE_native_IntMob_Erodib_K[source]
            * wflux_s2r
            * (*OpenWQ_vars.chemass)(source)(chemi)(ix_s,iy_s,iz_s);
        
        
        //##########################################
        // Set derivative for source and recipient 
        
        // Remove Chemical mass flux from SOURCE 
        (*OpenWQ_vars.d_chemass_dt_transp)(source)(chemi)(ix_s,iy_s,iz_s) 
            -= chemass_flux;

        // Add Chemical mass flux to RECIPIENT
        // if recipient == -1, then  it's an OUT-flux (loss from system)
        if (recipient == -1) continue;
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(chemi)(ix_r,iy_r,iz_r) 
            += chemass_flux;
    }

}
*/


