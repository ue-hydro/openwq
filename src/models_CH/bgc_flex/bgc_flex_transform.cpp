
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

#include "models_CH/headerfile_CH.hpp"

/* #################################################
// Compute each chemical transformation
// OPTIMIZED: uses pre-built BGC lookup map,
// pre-allocated chemass vector, and OpenMP parallelism
################################################# */
void OpenWQ_CH_model::bgc_flex_transform(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    unsigned int icmp){

    // Local variables
    unsigned int index_cons,index_prod;
    unsigned int num_transf;
    unsigned int num_BGCcycles;
    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
    std::string BGCcycles_name_icmp;
    std::string msg_string;

    // Find compartment icmp name from code (host hydrological model)
    const std::string CompName_icmp =
        OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);

    // Cache frequently accessed pointers
    auto* nativeFlex = OpenWQ_wqconfig.CH_model->NativeFlex;
    const bool is_time_step_0 = OpenWQ_hostModelconfig.is_time_step_0();
    const double time_step = OpenWQ_hostModelconfig.get_time_step();
    const unsigned int num_depend = OpenWQ_hostModelconfig.get_num_HydroDepend();

    // Get cycling_frameworks for compartment icomp
    try{
        std::vector<std::string> BGCcycles_icmp =
            OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"]
                [CompName_icmp]["CYCLING_FRAMEWORK"];

        num_BGCcycles = BGCcycles_icmp.size();

        /* ########################################
        // Loop over biogeochemical cycling frameworks
        ######################################## */

        for (unsigned int bgci=0;bgci<num_BGCcycles;bgci++){

            BGCcycles_name_icmp = BGCcycles_icmp.at(bgci);

            // OPTIMIZED: Use pre-built lookup map instead of O(N) scan
            const auto it = OpenWQ_wqconfig.bgc_cycle_to_transf_indices.find(BGCcycles_name_icmp);

            if (it == OpenWQ_wqconfig.bgc_cycle_to_transf_indices.end() || it->second.empty()){
                if (OpenWQ_wqconfig.invalid_bgc_entry_errmsg){
                    msg_string = "<OpenWQ> Unkown CYCLING_FRAMEWORK with name '"
                        + BGCcycles_name_icmp
                        + "' defined for compartment: "
                        + CompName_icmp;
                    OpenWQ_output.ConsoleLog(
                        OpenWQ_wqconfig, msg_string, true, true);
                }
                continue;
            }

            const std::vector<unsigned int>& transf_index = it->second;
            num_transf = transf_index.size();

            /* ########################################
            // Loop over transformations in biogeochemical cycle bgci
            ######################################## */

            for (unsigned int transi=0;transi<num_transf;transi++){

                const auto& bgc_info = nativeFlex->BGCexpressions_info.at(transf_index[transi]);
                index_cons = std::get<3>(bgc_info);
                index_prod = std::get<4>(bgc_info);
                const std::vector<unsigned int>& index_chemtransf = std::get<5>(bgc_info);
                const unsigned int num_chem_in_transf = index_chemtransf.size();

                // Pre-fetch cube references for consumed and produced species
                auto& d_chem_cons = (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(index_cons);
                auto& d_chem_prod = (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(index_prod);
                auto& chemass_cons = (*OpenWQ_vars.chemass)(icmp)(index_cons);

                /* ########################################
                // Loop over space: nx, ny, nz
                // Note: exprtk expressions are not thread-safe (they use
                // shared symbol tables), so the spatial loop remains sequential.
                // The pre-allocation optimization below still provides a
                // significant speedup by avoiding heap alloc per cell.
                ######################################## */

                // OPTIMIZED: Pre-allocate to correct size once, reuse via indexing
                nativeFlex->chemass_InTransfEq.resize(num_chem_in_transf);

                for (unsigned int ix=0;ix<nx;ix++){
                    for (unsigned int iy=0;iy<ny;iy++){
                        for (unsigned int iz=0;iz<nz;iz++){

                            // OPTIMIZED: overwrite by index instead of clear/push_back
                            for (unsigned int chem=0;chem<num_chem_in_transf;chem++){
                                nativeFlex->chemass_InTransfEq[chem] =
                                    (*OpenWQ_vars.chemass)(icmp)(index_chemtransf[chem])(ix,iy,iz);
                            }

                            // Update current dependencies for current x, y and z
                            for (unsigned int depi=0;depi<num_depend;depi++){
                                OpenWQ_hostModelconfig.set_dependVar_scalar_at(
                                    depi, OpenWQ_hostModelconfig.get_dependVar_at(depi,ix,iy,iz));
                            }

                            // Mass transfered: Consumed -> Produced (using exprtk)
                            double transf_mass = std::fmax(
                                nativeFlex->BGCexpressions_eq[transf_index[transi]].value(),
                                0.0);

                            // Guarantee that removed mass is not larger than existing mass
                            if (!is_time_step_0){
                                transf_mass = std::fmin(
                                        chemass_cons(ix,iy,iz),
                                        transf_mass * time_step);
                            }else{
                                transf_mass = 0.0;
                            }

                            // Update derivatives
                            d_chem_cons(ix,iy,iz) -= transf_mass;
                            d_chem_prod(ix,iy,iz) += transf_mass;

                        }
                    }
                }
            }
        }
    }catch(json::exception& e){

        if (OpenWQ_wqconfig.BGC_Transform_print_errmsg){
            msg_string =
                "<OpenWQ> No CYCLING_FRAMEWORK defined for compartment: "
                + CompName_icmp;
            OpenWQ_output.ConsoleLog(
                OpenWQ_wqconfig, msg_string, true, true);
        }

        return;
    }
}