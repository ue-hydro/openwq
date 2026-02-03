// Copyright 2026, Diogo Costa (diogo.costa@uevora.pt)
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


#include "compute/headerfile_compute.hpp"

/* #################################################
// Euler solver
################################################# */

void OpenWQ_compute::Solve_with_BE(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_CH_model& OpenWQ_CH_model){

    const unsigned int num_chem = 
        ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0)
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();
    const bool is_first_step = OpenWQ_hostModelconfig.is_first_interaction_step();

    /* #####################################################
    // Compartment loop - parallelized over compartments and chemicals
    ##################################################### */

    #pragma omp parallel num_threads(num_threads)
    {
        // Local variables - thread private
        unsigned int nx, ny, nz;
        unsigned int ix, iy, iz;
        double dm_dt_chem, dm_dt_trans, dm_ic, dm_ss, dm_ewf;

        #pragma omp for schedule(dynamic) collapse(2)
        for (unsigned int icmp = 0; icmp < num_comps; icmp++){
            for (unsigned int chemi = 0; chemi < num_chem; chemi++){

                // Get dimensions for this compartment (only once per icmp-chemi pair)
                nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
                ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
                nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

                // Pre-fetch pointers to reduce dereferencing overhead
                auto& chemass = (*OpenWQ_vars.chemass)(icmp)(chemi);
                auto& d_chemass_ic = (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi);
                auto& d_chemass_ss = (*OpenWQ_vars.d_chemass_ss)(icmp)(chemi);
                auto& d_chemass_ss_out = (*OpenWQ_vars.d_chemass_ss_out)(icmp)(chemi);
                auto& d_chemass_ewf = (*OpenWQ_vars.d_chemass_ewf)(icmp)(chemi);
                auto& d_chemass_ewf_out = (*OpenWQ_vars.d_chemass_ewf_out)(icmp)(chemi);
                auto& d_chemass_dt_chem = (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(chemi);
                auto& d_chemass_dt_chem_out = (*OpenWQ_vars.d_chemass_dt_chem_out)(icmp)(chemi);
                auto& d_chemass_dt_transp = (*OpenWQ_vars.d_chemass_dt_transp)(icmp)(chemi);
                auto& d_chemass_dt_transp_out = (*OpenWQ_vars.d_chemass_dt_transp_out)(icmp)(chemi);

                // X, Y, Z loops - process entire 3D grid for this compartment-chemical pair
                for (ix = 0; ix < nx; ix++){
                    for (iy = 0; iy < ny; iy++){
                        // Vectorization hint for innermost loop
                        #pragma omp simd
                        for (iz = 0; iz < nz; iz++){

                            // ####################################
                            // 1. Single-time change at simulation start
                            dm_ic = is_first_step ? d_chemass_ic(ix, iy, iz) : 0.0;

                            // ####################################
                            // 2. SS (Sink & Sources)
                            dm_ss = d_chemass_ss(ix, iy, iz);
                            d_chemass_ss_out(ix, iy, iz) += dm_ss;

                            // ####################################
                            // 3. EWF (External Water Fluxes)
                            dm_ewf = d_chemass_ewf(ix, iy, iz);
                            d_chemass_ewf_out(ix, iy, iz) += dm_ewf;

                            // ####################################
                            // 4. Dynamic change (derivatives): chemistry and transport
                            dm_dt_chem = d_chemass_dt_chem(ix, iy, iz);
                            d_chemass_dt_chem_out(ix, iy, iz) += dm_dt_chem;

                            dm_dt_trans = d_chemass_dt_transp(ix, iy, iz);
                            d_chemass_dt_transp_out(ix, iy, iz) += dm_dt_trans;

                            // ####################################
                            // 5. Apply all changes to state variable
                            double new_mass = chemass(ix, iy, iz) + dm_ic + dm_ss + dm_ewf + dm_dt_chem + dm_dt_trans;
                            
                            // Ensure non-negativity
                            chemass(ix, iy, iz) = (new_mass > 0.0) ? new_mass : 0.0;
                        }
                    }
                }
            }
        }
    }
}

void OpenWQ_compute::Solve_with_BE_Sediment(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output){

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Find sediment compartment index
    unsigned int sed_icmp = 0;
    const std::string& erod_transp_cmpt = OpenWQ_wqconfig.TS_model->ErodTranspCmpt;
    
    for (unsigned int icmp = 0; icmp < num_comps; icmp++){
        if (erod_transp_cmpt.compare(OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)) == 0) {
            sed_icmp = icmp;
            break;
        }
    }

    // Get dimensions for sediment compartment
    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(sed_icmp);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(sed_icmp);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(sed_icmp);

    // Pre-fetch array references
    auto& sedmass = *OpenWQ_vars.sedmass;
    auto& d_sedmass_dt = *OpenWQ_vars.d_sedmass_dt;
    auto& d_sedmass_mobilized_dt = *OpenWQ_vars.d_sedmass_mobilized_dt;

    /* #####################################################
    // Parallel 3D loop with better granularity
    ##################################################### */

    #pragma omp parallel num_threads(num_threads)
    {
        unsigned int ix, iy, iz;
        double dm_dt;

        // Parallelize over flattened 3D index for better load balancing
        #pragma omp for schedule(static)
        for (unsigned int idx = 0; idx < nx * ny * nz; idx++){
            ix = idx / (ny * nz);
            iy = (idx / nz) % ny;
            iz = idx % nz;

            dm_dt = d_sedmass_dt(ix, iy, iz) + d_sedmass_mobilized_dt(ix, iy, iz);

            // Apply change and ensure non-negativity
            double new_sedmass = sedmass(ix, iy, iz) + dm_dt;
            sedmass(ix, iy, iz) = (new_sedmass > 0.0) ? new_sedmass : 0.0;
        }
    }
}