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
#include "models_CH/headerfile_CH.hpp"
#include "models_TS/headerfile_TS.hpp"
#include "utils/headerfile_UTILS.hpp"

#include <cstring> // for std::memcpy (OPTIMIZED #10)
#include <sundials/sundials_core.hpp>
#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>

/* #################################################
// Utilities for CVode solver - structure to hold OpenWQ structures
################################################# */
struct UserData {
    OpenWQ_hostModelconfig& hostModelconfig;
    OpenWQ_wqconfig& wqconfig;
    OpenWQ_vars& vars;
    OpenWQ_json& json;
    OpenWQ_output& output;
    OpenWQ_CH_model& chem;
    std::vector<std::tuple<int,int,int,int,int,int,int,int,double,double,std::string>>& sediment_calls;
    OpenWQ_TS_model& TS_model;
    OpenWQ_utils& utils;

    // Cache frequently accessed values
    unsigned int num_chem;
    unsigned int num_comps;
    unsigned int num_threads;
    bool is_first_step;
    double time_step;

    // Pre-computed compartment dimensions
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> comp_dims;
    std::vector<unsigned int> comp_offsets; // Flattened index offsets

    // OPTIMIZED #7: cached sediment compartment index
    unsigned int sed_icmp;
};

/* #################################################
// Utilities for CVode solver - function to compute total flux
################################################# */
int totalFlux(sunrealtype t, N_Vector u, N_Vector f, void* udata) {

    auto& user_data = *(static_cast<UserData*>(udata));
    
    sunrealtype* uvec_data = N_VGetArrayPointer(u);
    sunrealtype* fvec_data = N_VGetArrayPointer(f);

    // Reset chemistry derivatives in parallel
    #pragma omp parallel num_threads(user_data.num_threads)
    {
        #pragma omp for schedule(static)
        for (unsigned int idx = 0; idx < user_data.num_comps * user_data.num_chem; idx++){
            const unsigned int icmp = idx / user_data.num_chem;
            const unsigned int chemi = idx % user_data.num_chem;
            (*user_data.vars.d_chemass_dt_chem)(icmp)(chemi).zeros();
        }
    }

    // OPTIMIZED #10: Copy data from SUNDIALS vector to chemass using memcpy
    for (unsigned int icmp = 0; icmp < user_data.num_comps; icmp++){
        const auto& [nx, ny, nz] = user_data.comp_dims[icmp];
        const unsigned int offset = user_data.comp_offsets[icmp];
        const unsigned int cube_size = nx * ny * nz;

        for (unsigned int chemi = 0; chemi < user_data.num_chem; chemi++){
            const unsigned int chem_offset = offset + chemi * cube_size;
            std::memcpy((*user_data.vars.chemass)(icmp)(chemi).memptr(),
                       uvec_data + chem_offset,
                       cube_size * sizeof(double));
        }
    }

    // Run chemistry model
    user_data.chem.CH_driver_run(user_data.json, user_data.vars, user_data.wqconfig, 
                                  user_data.hostModelconfig, user_data.output);

    // Compute flux derivatives in parallel
    #pragma omp parallel num_threads(user_data.num_threads)
    {
        unsigned int nx, ny, nz, ix, iy, iz;
        double dm_dt_chem, dm_dt_trans, dm_ic, dm_ss, dm_ewf;

        #pragma omp for schedule(dynamic)
        for (unsigned int icmp = 0; icmp < user_data.num_comps; icmp++){
            
            std::tie(nx, ny, nz) = user_data.comp_dims[icmp];

            for (unsigned int chemi = 0; chemi < user_data.num_chem; chemi++){
                
                // Pre-fetch array references
                auto& d_chemass = (*user_data.vars.d_chemass)(icmp)(chemi);
                auto& d_chemass_ic = (*user_data.vars.d_chemass_ic)(icmp)(chemi);
                auto& d_chemass_ss = (*user_data.vars.d_chemass_ss)(icmp)(chemi);
                auto& d_chemass_ss_out = (*user_data.vars.d_chemass_ss_out)(icmp)(chemi);
                auto& d_chemass_ewf = (*user_data.vars.d_chemass_ewf)(icmp)(chemi);
                auto& d_chemass_ewf_out = (*user_data.vars.d_chemass_ewf_out)(icmp)(chemi);
                auto& d_chemass_dt_chem = (*user_data.vars.d_chemass_dt_chem)(icmp)(chemi);
                auto& d_chemass_dt_chem_out = (*user_data.vars.d_chemass_dt_chem_out)(icmp)(chemi);
                auto& d_chemass_dt_transp = (*user_data.vars.d_chemass_dt_transp)(icmp)(chemi);
                auto& d_chemass_dt_transp_out = (*user_data.vars.d_chemass_dt_transp_out)(icmp)(chemi);

                d_chemass.zeros();

                for (ix = 0; ix < nx; ix++){
                    for (iy = 0; iy < ny; iy++){
                        #pragma omp simd
                        for (iz = 0; iz < nz; iz++){

                            // 1. Initial conditions (only first step)
                            dm_ic = user_data.is_first_step ? d_chemass_ic(ix, iy, iz) : 0.0;

                            // 2. SS (Sink & Sources)
                            dm_ss = d_chemass_ss(ix, iy, iz);
                            d_chemass_ss_out(ix, iy, iz) += dm_ss;

                            // 3. EWF (External Water Fluxes)
                            dm_ewf = d_chemass_ewf(ix, iy, iz);
                            d_chemass_ewf_out(ix, iy, iz) += dm_ewf;

                            // 4. Chemistry derivatives
                            dm_dt_chem = d_chemass_dt_chem(ix, iy, iz);
                            d_chemass_dt_chem_out(ix, iy, iz) += dm_dt_chem;

                            // Transport derivatives
                            dm_dt_trans = d_chemass_dt_transp(ix, iy, iz);
                            d_chemass_dt_transp_out(ix, iy, iz) += dm_dt_trans;

                            // 5. Total derivative
                            d_chemass(ix, iy, iz) = dm_ic + dm_ss + dm_ewf + dm_dt_chem + dm_dt_trans;
                        }
                    }
                }
            }
        }
    }

    // OPTIMIZED #10: Copy derivatives to SUNDIALS f vector using memcpy + scaling
    const double inv_dt = 1.0 / user_data.time_step;

    for (unsigned int icmp = 0; icmp < user_data.num_comps; icmp++){
        const auto& [nx, ny, nz] = user_data.comp_dims[icmp];
        const unsigned int offset = user_data.comp_offsets[icmp];
        const unsigned int cube_size = nx * ny * nz;

        for (unsigned int chemi = 0; chemi < user_data.num_chem; chemi++){
            const unsigned int chem_offset = offset + chemi * cube_size;
            const double* src = (*user_data.vars.d_chemass)(icmp)(chemi).memptr();

            // memcpy then scale (avoids triple-nested loop overhead)
            std::memcpy(fvec_data + chem_offset, src, cube_size * sizeof(double));
            for (unsigned int i = 0; i < cube_size; i++){
                fvec_data[chem_offset + i] *= inv_dt;
            }
        }
    }

    return 0;
}

int sedimentFlux(sunrealtype t, N_Vector u, N_Vector f, void* udata) {

    auto& user_data = *(static_cast<UserData*>(udata));
    
    sunrealtype* uvec_data = N_VGetArrayPointer(u);
    sunrealtype* fvec_data = N_VGetArrayPointer(f);

    // Reset derivatives
    (*user_data.vars.d_sedmass_mobilized_dt).zeros();
    (*user_data.vars.d_sedmass_dt).zeros();

    // OPTIMIZED #7: use cached sediment compartment index
    const unsigned int sed_icmp = user_data.sed_icmp;

    const auto& [nx, ny, nz] = user_data.comp_dims[sed_icmp];
    auto& sedmass = *user_data.vars.sedmass;
    const unsigned int sed_size = nx * ny * nz;

    // OPTIMIZED #10: Copy from SUNDIALS vector to sedmass using memcpy
    std::memcpy(sedmass.memptr(), uvec_data, sed_size * sizeof(double));

    // Process sediment transport calls
    for (const auto& [source, ix_s, iy_s, iz_s, recipient, ix_r, iy_r, iz_r, wflux_s2r, wmass_source, TS_type] 
         : user_data.sediment_calls) {
        user_data.TS_model.TS_driver_run(user_data.hostModelconfig,
                                          user_data.wqconfig,
                                          user_data.vars,
                                          user_data.utils,
                                          user_data.output,
                                          source, ix_s, iy_s, iz_s,
                                          recipient, ix_r, iy_r, iz_r,
                                          wflux_s2r, wmass_source, TS_type);
    }

    // OPTIMIZED #10: Copy derivatives to SUNDIALS f vector using raw pointers
    const double* dt_ptr = user_data.vars.d_sedmass_dt->memptr();
    const double* mob_ptr = user_data.vars.d_sedmass_mobilized_dt->memptr();

    for (unsigned int i = 0; i < sed_size; i++){
        fvec_data[i] = dt_ptr[i] + mob_ptr[i];
    }

    return 0;
}

/* #################################################
// Solving with CVode
################################################# */
void OpenWQ_compute::Solve_with_CVode(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_CH_model& OpenWQ_CH_model,
    OpenWQ_TS_model& OpenWQ_TS_model,
    OpenWQ_utils& OpenWQ_utils){
    
    // OPTIMIZED: use cached num_chem instead of string comparison
    const unsigned int num_chem = OpenWQ_wqconfig.is_native_bgc_flex
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Pre-compute compartment dimensions and offsets
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> comp_dims;
    std::vector<unsigned int> comp_offsets;
    unsigned int system_size = 0;

    comp_dims.reserve(num_comps);
    comp_offsets.reserve(num_comps);

    for (unsigned int icmp = 0; icmp < num_comps; icmp++) {
        const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

        comp_dims.emplace_back(nx, ny, nz);
        comp_offsets.push_back(system_size);
        system_size += nx * ny * nz * num_chem;
    }

    // Initialize UserData with cached values
    UserData user_data = {
        OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_vars, OpenWQ_json,
        OpenWQ_output, OpenWQ_CH_model, sediment_calls, OpenWQ_TS_model, OpenWQ_utils,
        num_chem, num_comps,
        static_cast<unsigned int>(OpenWQ_wqconfig.get_num_threads_requested()),
        OpenWQ_hostModelconfig.is_first_interaction_step(),
        OpenWQ_hostModelconfig.get_time_step(),
        comp_dims, comp_offsets,
        0 // sed_icmp (not used for chemistry solve)
    };

    // Create SUNDIALS context and vector
    sundials::Context sunctx;
    N_Vector u = N_VNew_Serial(system_size, sunctx);
    sunrealtype* udata_ptr = N_VGetArrayPointer(u);

    // OPTIMIZED #10: Initialize state vector using memcpy where possible
    for (unsigned int icmp = 0; icmp < num_comps; icmp++){
        const auto& [nx, ny, nz] = comp_dims[icmp];
        const unsigned int offset = comp_offsets[icmp];
        const unsigned int cube_size = nx * ny * nz;

        for (unsigned int chemi = 0; chemi < num_chem; chemi++){
            const unsigned int chem_offset = offset + chemi * cube_size;

            if (user_data.is_first_step) {
                // Armadillo stores data contiguously in column-major order
                std::memcpy(udata_ptr + chem_offset,
                           (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi).memptr(),
                           cube_size * sizeof(double));
            } else {
                std::memcpy(udata_ptr + chem_offset,
                           (*OpenWQ_vars.chemass)(icmp)(chemi).memptr(),
                           cube_size * sizeof(double));
            }
        }
    }

    // SUNDIALS/CVode setup
    void* cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
    CVodeInit(cvode_mem, totalFlux, 0, u);
    CVodeSetUserData(cvode_mem, (void*)&user_data);

    SUNLinearSolver LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 0, sunctx);
    CVodeSetLinearSolver(cvode_mem, LS, NULL);
    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    // Solve
    sunrealtype t;
    if (OpenWQ_hostModelconfig.get_time_step() > 0) {
        CVode(cvode_mem, OpenWQ_hostModelconfig.get_time_step(), u, &t, CV_NORMAL);
    }

    // OPTIMIZED #10: Extract solution using memcpy
    for (unsigned int icmp = 0; icmp < num_comps; icmp++){
        const auto& [nx, ny, nz] = comp_dims[icmp];
        const unsigned int offset = comp_offsets[icmp];
        const unsigned int cube_size = nx * ny * nz;

        for (unsigned int chemi = 0; chemi < num_chem; chemi++){
            const unsigned int chem_offset = offset + chemi * cube_size;
            std::memcpy((*OpenWQ_vars.chemass)(icmp)(chemi).memptr(),
                       udata_ptr + chem_offset,
                       cube_size * sizeof(double));
        }
    }

    // Cleanup
    N_VDestroy(u);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
}

/* #################################################
// Solving Sediment with CVode
################################################# */
void OpenWQ_compute::Solve_with_CVode_Sediment(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_CH_model& OpenWQ_CH_model,
    OpenWQ_TS_model& OpenWQ_TS_model,
    OpenWQ_utils& OpenWQ_utils){
    
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // OPTIMIZED #7: Find sediment compartment index once
    unsigned int sed_icmp = 0;
    const std::string& erod_transp_cmpt = OpenWQ_wqconfig.TS_model->ErodTranspCmpt;

    for (unsigned int icmp = 0; icmp < num_comps; icmp++){
        if (erod_transp_cmpt.compare(OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)) == 0) {
            sed_icmp = icmp;
            break;
        }
    }

    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(sed_icmp);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(sed_icmp);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(sed_icmp);
    const unsigned int system_size = nx * ny * nz;

    // Pre-compute compartment dimensions for UserData
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> comp_dims;
    comp_dims.reserve(num_comps);
    for (unsigned int icmp = 0; icmp < num_comps; icmp++) {
        comp_dims.emplace_back(
            OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp),
            OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp),
            OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp)
        );
    }

    std::vector<unsigned int> comp_offsets; // Not used for sediment but needed for struct

    // OPTIMIZED: use cached flag instead of string comparison
    const unsigned int num_chem = OpenWQ_wqconfig.is_native_bgc_flex
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    UserData user_data = {
        OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_vars, OpenWQ_json,
        OpenWQ_output, OpenWQ_CH_model, sediment_calls, OpenWQ_TS_model, OpenWQ_utils,
        num_chem, num_comps,
        static_cast<unsigned int>(OpenWQ_wqconfig.get_num_threads_requested()),
        OpenWQ_hostModelconfig.is_first_interaction_step(),
        OpenWQ_hostModelconfig.get_time_step(),
        comp_dims, comp_offsets,
        sed_icmp  // OPTIMIZED #7: pass cached sediment compartment index
    };

    // Create SUNDIALS context and vector
    sundials::Context sunctx;
    N_Vector u = N_VNew_Serial(system_size, sunctx);
    sunrealtype* udata_ptr = N_VGetArrayPointer(u);

    // OPTIMIZED #10: Initialize state vector using memcpy
    std::memcpy(udata_ptr, OpenWQ_vars.sedmass->memptr(), system_size * sizeof(double));

    // SUNDIALS/CVode setup
    void* cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
    CVodeInit(cvode_mem, sedimentFlux, 0, u);
    CVodeSetUserData(cvode_mem, (void*)&user_data);

    SUNLinearSolver LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 0, sunctx);
    CVodeSetLinearSolver(cvode_mem, LS, NULL);
    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    // Solve
    sunrealtype t;
    if (OpenWQ_hostModelconfig.get_time_step() > 0) {
        CVode(cvode_mem, OpenWQ_hostModelconfig.get_time_step(), u, &t, CV_NORMAL);
    }

    // OPTIMIZED #10: Extract solution using memcpy
    std::memcpy(OpenWQ_vars.sedmass->memptr(), udata_ptr, system_size * sizeof(double));

    // Cleanup
    N_VDestroy(u);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
}