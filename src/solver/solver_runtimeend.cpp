
// Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
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


#include "OpenWQ_solver.hpp"
#include "chem/OpenWQ_chem.hpp"

#include <sundials/sundials_core.hpp>
#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>

/* #################################################
// General numerical solver
################################################# */
struct UserData {
    OpenWQ_hostModelconfig& hostModelconfig;
    OpenWQ_wqconfig& wqconfig;
    OpenWQ_vars& vars;
    OpenWQ_json& json;
    OpenWQ_output& output;
    OpenWQ_chem& chem;
};

int totalFlux(sunrealtype t, N_Vector u, N_Vector f, void* udata) {

    auto user_data = *(static_cast<UserData*>(udata));

    // Local variables
    unsigned int nx, ny, nz;    // interactive compartment domain dimensions
    unsigned int ix, iy, iz;    // interactive compartment domain cell indexes
    double dm_dt_chem;          // ineractive final derivative (mass) - chemistry transformations
    double dm_dt_trans;         // ineractive final derivative (mass) - transport
    double dm_ic;               // interactive dynamic change to state-variable at start of simulations
    double dm_ss;               // interactive ss load/ink (mass)
    double dm_ewf;              // interactive ewf load/ink (mass)
    int idx=0;
    sunrealtype* uvec_data;
    sunrealtype* fvec_data;


    /* #####################################################
    // Compartment loop
    ##################################################### */

    uvec_data = N_VGetArrayPointer(u);
    fvec_data = N_VGetArrayPointer(f);


    #pragma omp parallel for collapse(2) num_threads(user_data.wqconfig.get_num_threads_requested())
        for (unsigned int icmp=0;icmp<user_data.hostModelconfig.get_num_HydroComp();icmp++){

            // Chemical loop
            for (unsigned int chemi=0;chemi<(user_data.wqconfig.BGC_general_num_chem);chemi++){

                // Reset derivatives of state-variables to zero after each time interaction
                (*user_data.vars.d_chemass_dt_chem)(icmp)(chemi).zeros();

            }

        }
    for (unsigned int icmp=0;icmp<user_data.hostModelconfig.get_num_HydroComp();icmp++){

        // Dimensions for compartment icmp
        nx = user_data.hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
        ny = user_data.hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
        nz = user_data.hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements

        // Chemical loop
        for (unsigned int chemi=0;chemi<(user_data.wqconfig.BGC_general_num_chem);chemi++){

            // X, Y, Z loops
            for (ix=0;ix<nx;ix++){
                for (iy=0;iy<ny;iy++){
                    for (iz=0;iz<nz;iz++){                        
                        (*user_data.vars.chemass)(icmp)(chemi)(ix,iy,iz) = uvec_data[idx];
                        idx++;
                    }
                }
            }
        }
    }

    user_data.chem.Run(user_data.json, user_data.vars, user_data.wqconfig, user_data.hostModelconfig, user_data.output);

    #pragma omp parallel for private (nx, ny, nz, ix, iy, iz, dm_ic, dm_ss, dm_ewf, dm_dt_chem, dm_dt_trans) num_threads(user_data.wqconfig.get_num_threads_requested())
    for (unsigned int icmp=0;icmp<user_data.hostModelconfig.get_num_HydroComp();icmp++){

        // Dimensions for compartment icmp
        nx = user_data.hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
        ny = user_data.hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
        nz = user_data.hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements

        // Chemical loop
        for (unsigned int chemi=0;chemi<(user_data.wqconfig.BGC_general_num_chem);chemi++){
            (*user_data.vars.d_chemass)(icmp)(chemi).zeros();
            // X, Y, Z loops
            for (ix=0;ix<nx;ix++){
                for (iy=0;iy<ny;iy++){
                    for (iz=0;iz<nz;iz++){

                        // ####################################
                        // 1
                        // Single-time change at simulation start
                        //  (IC and input unit conversions)
                        if (user_data.hostModelconfig.is_first_interaction_step()){

                            dm_ic = 
                                (*user_data.vars.d_chemass_ic)(icmp)(chemi)(ix,iy,iz);

                        }else{
                            dm_ic = 0.0f;
                        };

                        // ####################################
                        // 2
                        // SS (Sink & Sources)
                        // Discrete or continuous chemical load

                        dm_ss = (*user_data.vars.d_chemass_ss)(icmp)(chemi)(ix,iy,iz);
                        // updating cumulative calc for output in debug mode
                        (*user_data.vars.d_chemass_ss_out)(icmp)(chemi)(ix,iy,iz) += dm_ss;

                        // ####################################
                        // 3
                        // EWF (External Water Fluxes)
                        // Continuous load caused by external fluxes

                        dm_ewf = (*user_data.vars.d_chemass_ewf)(icmp)(chemi)(ix,iy,iz);
                        // updating cumulative calc for output in debug mode
                        (*user_data.vars.d_chemass_ewf_out)(icmp)(chemi)(ix,iy,iz) += dm_ewf;

                        // ####################################
                        // 4
                        // Dynamic change (derivatives): chemistry and tranport
                        // It has already been multiplied by hostmodel dt (OpenWQ_hostModelconfig.time_step)

                        // Chemistry
                        dm_dt_chem = (*user_data.vars.d_chemass_dt_chem)(icmp)(chemi)(ix,iy,iz);
                        // updating cumulative calc for output in debug mode
                        // No need to multiply by timestep because that has been done in chem module
                        (*user_data.vars.d_chemass_dt_chem_out)(icmp)(chemi)(ix,iy,iz) += dm_dt_chem;
                        // Transport (no need to multiply by timestep because the water flux is usually
                        // is provided as a water volume, and not as a water volume per unit of time)
                        dm_dt_trans = (*user_data.vars.d_chemass_dt_transp)(icmp)(chemi)(ix,iy,iz); 
                        // updating cumulative calc for output in debug mode
                        (*user_data.vars.d_chemass_dt_transp_out)(icmp)(chemi)(ix,iy,iz) += dm_dt_trans;  
                        
                        // ####################################
                        // 5
                        // Apply all changes to state variable (all at once)

                        // Change state-variable
                        //(*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz) += 
                        //    (   dm_ic               // Mass change due to IC (= to zero if not 1st timestep)
                        //        + dm_ss             // Mass change due to SS input
                        //        + dm_ewf            // Mass change due to EWF input
                        //        + dm_dt_chem        // Mass change due to chemistry                            
                        //        + dm_dt_trans);     // Mass change due to transport
                        (*user_data.vars.d_chemass)(icmp)(chemi)(ix,iy,iz) = dm_ic+dm_ss+dm_ewf+dm_dt_chem + dm_dt_trans;
                        /*if((*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz) < 0){
                            (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz) = 0;
                        }*/
                    }
                }
            }
        }
    }

    idx = 0;
    for (unsigned int icmp=0;icmp<user_data.hostModelconfig.get_num_HydroComp();icmp++){

        // Dimensions for compartment icmp
        nx = user_data.hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
        ny = user_data.hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
        nz = user_data.hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements

        // Chemical loop
        for (unsigned int chemi=0;chemi<(user_data.wqconfig.BGC_general_num_chem);chemi++){

            // X, Y, Z loops
            for (ix=0;ix<nx;ix++){
                for (iy=0;iy<ny;iy++){
                    for (iz=0;iz<nz;iz++){                        
                        fvec_data[idx] = (*user_data.vars.d_chemass)(icmp)(chemi)(ix,iy,iz)/user_data.hostModelconfig.get_time_step();
                        idx++;
                    }
                }
            }
        }
    }
    return 0;

}

void OpenWQ_solver::Numerical_Solver(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_chem& OpenWQ_chem){
    
    // Local variables
    unsigned int nx, ny, nz;    // interactive compartment domain dimensions
    unsigned int ix, iy, iz;    // interactive compartment domain cell indexes
    double dm_dt_chem;          // ineractive final derivative (mass) - chemistry transformations
    double dm_dt_trans;         // ineractive final derivative (mass) - transport
    double ic;               // interactive dynamic change to state-variable at start of simulations
    double dm_ss;               // interactive ss load/ink (mass)
    double dm_ewf;              // interactive ewf load/ink (mass)
    unsigned int system_size = 0;
    sundials::Context sunctx;
    N_Vector u;
    SUNLinearSolver LS;
    int retval;
    int idx=0;
    sunrealtype* udata;
    void* cvode_mem;
    sunrealtype t;

    UserData user_data = {OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_vars, OpenWQ_json, OpenWQ_output, OpenWQ_chem};

    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++) {
        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
        system_size += nx*ny*nz*OpenWQ_wqconfig.BGC_general_num_chem;
    }

    //retval = SUNContext_Create(NULL, &sunctx);
    u = N_VNew_Serial(system_size, sunctx);
    udata = N_VGetArrayPointer(u);

    // Loop  over compartments
    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){

        // Dimensions for compartment icmp
        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements

        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements

        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements

        // Chemical loop
        for (unsigned int chemi=0;chemi<(OpenWQ_wqconfig.BGC_general_num_chem);chemi++){

            // X, Y, Z loops
            for (ix=0;ix<nx;ix++){
                for (iy=0;iy<ny;iy++){
                    for (iz=0;iz<nz;iz++){  
                        if (user_data.hostModelconfig.is_first_interaction_step()){

                            ic = 
                                (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi)(ix,iy,iz);

                        } else {
                            ic = (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz);
                        } 
                        udata[idx] = ic;
                        idx++;

                    }
                }
            }
        }
    }

    // SUNDIALS/CVode implementation (Victoria Guenter)
    cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
    CVodeInit(cvode_mem, totalFlux, 0, u);

    CVodeSetUserData(cvode_mem, (void*) &user_data);

    LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 0, sunctx);

    CVodeSetLinearSolver(cvode_mem, LS, NULL); 
    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    if (OpenWQ_hostModelconfig.get_time_step() > 0) {
        retval = CVode(cvode_mem, OpenWQ_hostModelconfig.get_time_step(), u, &t, CV_NORMAL);
    }

    idx=0;

    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){

        // Dimensions for compartment icmp
        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements

        // Chemical loop
        for (unsigned int chemi=0;chemi<(OpenWQ_wqconfig.BGC_general_num_chem);chemi++){

            // X, Y, Z loops
            for (ix=0;ix<nx;ix++){
                for (iy=0;iy<ny;iy++){
                    for (iz=0;iz<nz;iz++){     

                        (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz) = udata[idx];
                        idx++;

                    }
                }
            }
        }
    }

}
