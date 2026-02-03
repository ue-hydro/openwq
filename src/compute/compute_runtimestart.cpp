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
// Reset derivatives 
// (needed a the start of each temporal iteraction)
################################################# */

void OpenWQ_compute::Reset_Deriv(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    bool inst_deriv_flag,
    bool cum_deriv_flag){

    const unsigned int num_chem = 
        ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0) 
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Reset Instant Derivative
    if (inst_deriv_flag)
    {   
        // Chemistry variables - Single parallelized loop with better work distribution
        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp for schedule(static) nowait
            for (unsigned int idx = 0; idx < num_comps * num_chem; idx++){
                const unsigned int icmp = idx / num_chem;
                const unsigned int chemi = idx % num_chem;

                // Reset derivatives of state-variables to zero after each time interaction
                (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi).zeros();
                (*OpenWQ_vars.d_chemass_ss)(icmp)(chemi).zeros();
                (*OpenWQ_vars.d_chemass_ewf)(icmp)(chemi).zeros();
                (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(chemi).zeros();
                (*OpenWQ_vars.d_chemass_dt_transp)(icmp)(chemi).zeros();
            }

            // Sediment transport variables - executed by single thread
            #pragma omp single
            {
                (*OpenWQ_vars.d_sedmass_mobilized_dt).zeros();
                (*OpenWQ_vars.d_sedmass_dt).zeros();
            }
        }
    }

    // Reset Cumulative Derivative
    // only used for exporting results in debug mode
    if (cum_deriv_flag)
    {
        // Single flattened loop for better load balancing
        #pragma omp parallel for schedule(static) num_threads(num_threads)
        for (unsigned int idx = 0; idx < num_comps * num_chem; idx++){
            const unsigned int icmp = idx / num_chem;
            const unsigned int chemi = idx % num_chem;

            // Reset derivatives of state-variables to zero after each time interaction
            (*OpenWQ_vars.d_chemass_ss_out)(icmp)(chemi).zeros();
            (*OpenWQ_vars.d_chemass_ewf_out)(icmp)(chemi).zeros();
            (*OpenWQ_vars.d_chemass_dt_chem_out)(icmp)(chemi).zeros();
            (*OpenWQ_vars.d_chemass_dt_transp_out)(icmp)(chemi).zeros();
        }
    }
}

// ########################################
// Reset EWF conc 
// Specially needed for discrete conc requests
// ########################################
void OpenWQ_compute::Reset_EWFconc(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars){

    const unsigned int num_chem = 
        ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0)
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    const unsigned int num_ewf = OpenWQ_hostModelconfig.get_num_HydroExtFlux();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Single flattened loop for better load balancing
    #pragma omp parallel for schedule(static) num_threads(num_threads)
    for (unsigned int idx = 0; idx < num_ewf * num_chem; idx++){
        const unsigned int ewfi = idx / num_chem;
        const unsigned int chemi = idx % num_chem;

        // Reset ewf_conc after each iteraction
        // Specially needed for discrete conc requests
        (*OpenWQ_vars.ewf_conc)(ewfi)(chemi).zeros();
    }
}