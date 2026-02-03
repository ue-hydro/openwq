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

#include "headerfile_INIT.hpp"

/* #################################################
// Get Time variables
################################################# */
void OpenWQ_initiate::setTimeVars(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    time_t simtime){

    // Update interaction number
    OpenWQ_hostModelconfig.increment_interaction_step();

    // Update other time variables
    if (OpenWQ_hostModelconfig.is_first_interaction_step()){
        // Initiate nexttime_out with simtime at start of simulation
        OpenWQ_wqconfig.set_nexttime_out(simtime);
    } else {
        // Update time step
        OpenWQ_hostModelconfig.set_time_step(simtime - OpenWQ_wqconfig.get_time_previous());
    }

    // Update time previous
    OpenWQ_wqconfig.set_time_previous(simtime);
}

/* #################################################
// Read and Set simulations conditions and options
################################################# */
void OpenWQ_initiate::readSet(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
        
    /* #################################################
    // Read and set IC conditions
    ################################################# */
   
    if ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0) {
        const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
        
        for (unsigned int icmp = 0; icmp < num_comps; icmp++){
            setIC_driver(
                OpenWQ_json,
                OpenWQ_vars,
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_utils,
                OpenWQ_units,
                OpenWQ_output,
                icmp);
        }
    } else {
        // ICs for PHREEQC come from PHREEQC
        setIC_phreeqc(
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_utils,
            OpenWQ_units,
            OpenWQ_output);
    }
}

void OpenWQ_initiate::setIC_phreeqc(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output){

    std::string errorMsgIdentifier = " Config file";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config, "BIOGEOCHEMISTRY_CONFIGURATION",
        errorMsgIdentifier,
        true);

    const int nxyz = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetGridCellCount();
    const int nc = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponentCount();
    const unsigned int num_chem = OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    std::vector<double> c(nxyz * num_chem);
    std::vector<int> ic1(nxyz * 7, -1);

    int indx = 0;
    for (unsigned int icmp = 0; icmp < num_comps; icmp++){
        
        const std::string CompName_icmp = 
            OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
        
        // Initialize all IDs to -1
        int SOLUTIONS_id = -1;
        int EQUILIBRIUM_PHASES_id = -1;
        int EXCHANGE_id = -1;
        int SURFACE_id = -1;
        int GAS_PHASE_id = -1;
        int SOLID_SOLUTIONS_id = -1;
        int KINETICS_id = -1;
        
        if (OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"].contains(CompName_icmp)){
            const auto& bgc_config = OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"][CompName_icmp];
            
            if (bgc_config.contains("SOLUTIONS"))
                SOLUTIONS_id = bgc_config["SOLUTIONS"];
            if (bgc_config.contains("EQUILIBRIUM_PHASES"))
                EQUILIBRIUM_PHASES_id = bgc_config["EQUILIBRIUM_PHASES"];
            if (bgc_config.contains("EXCHANGE"))
                EXCHANGE_id = bgc_config["EXCHANGE"];
            if (bgc_config.contains("SURFACE"))
                SURFACE_id = bgc_config["SURFACE"];
            if (bgc_config.contains("GAS_PHASE"))
                GAS_PHASE_id = bgc_config["GAS_PHASE"];
            if (bgc_config.contains("SOLID_SOLUTIONS"))
                SOLID_SOLUTIONS_id = bgc_config["SOLID_SOLUTIONS"];
            if (bgc_config.contains("KINETICS"))
                KINETICS_id = bgc_config["KINETICS"];
        } else {
            std::string msg_string = 
                "<OpenWQ> PHREEQC conditions not defined for compartment:" 
                + CompName_icmp + " (set to undefined)";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }

        const int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        const int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        const int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

        for (int ix = 0; ix < nx; ix++){
            for (int iy = 0; iy < ny; iy++){
                for (int iz = 0; iz < nz; iz++){
                    ic1[indx] = SOLUTIONS_id;
                    ic1[nxyz + indx] = EQUILIBRIUM_PHASES_id;
                    ic1[2 * nxyz + indx] = EXCHANGE_id;
                    ic1[3 * nxyz + indx] = SURFACE_id;
                    ic1[4 * nxyz + indx] = GAS_PHASE_id;
                    ic1[5 * nxyz + indx] = SOLID_SOLUTIONS_id;
                    ic1[6 * nxyz + indx] = KINETICS_id;
                    indx++;
                }
            }
        }
    }

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->InitialPhreeqc2Module(ic1);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTime(0.0);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTimeStep(0);
    IRM_RESULT result = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunCells();
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetConcentrations(c);

    // Parallelized concentration extraction
    #pragma omp parallel num_threads(OpenWQ_wqconfig.get_num_threads_requested())
    {
        #pragma omp for schedule(static)
        for (unsigned int icmp = 0; icmp < num_comps; icmp++){
            
            const int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
            const int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
            const int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
            
            // Calculate starting index for this compartment
            int comp_start_idx = 0;
            for (unsigned int prev_icmp = 0; prev_icmp < icmp; prev_icmp++){
                comp_start_idx += OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(prev_icmp) *
                                  OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(prev_icmp) *
                                  OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(prev_icmp);
            }

            for (unsigned int chemi = 0; chemi < num_chem; chemi++){
                int local_indx = 0;
                
                for (int ix = 0; ix < nx; ix++){
                    for (int iy = 0; iy < ny; iy++){
                        for (int iz = 0; iz < nz; iz++){
                            const double den = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ix, iy, iz);
                            const double volume = (den == 0) ? 1.0 : den;
                            
                            (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi)(ix, iy, iz) = 
                                c[chemi * nxyz + comp_start_idx + local_indx] * volume;
                            
                            local_indx++;
                        }
                    }
                }
            }
        }
    }
}

/* #################################################
// Read and set IC conditions
################################################# */
void OpenWQ_initiate::setIC_driver(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    const int icmp){

    // Find compartment name
    const std::string CompName_icmp = 
        OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);

    // Check if BIOGEOCHEMISTRY_CONFIGURATION key exists
    std::string errorMsgIdentifier = " Config file";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config, "BIOGEOCHEMISTRY_CONFIGURATION",
        errorMsgIdentifier,
        true);

    const unsigned int num_chem = OpenWQ_wqconfig.CH_model->NativeFlex->num_chem;

    /* ########################################
    // Loop over chemical species
    ######################################## */

    for (unsigned int chemi = 0; chemi < num_chem; chemi++){

        // Get chem name
        const std::string chemname = 
            (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi];

        // Default all IC chemical mass to zero
        (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi).zeros();

        /* ########################################
        // IC conditions provided
        ######################################## */

        if (OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"].contains(CompName_icmp)){
            
            const std::string ic_data_format = 
                OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"][CompName_icmp]
                ["INITIAL_CONDITIONS"]["DATA_FORMAT"];

            if (ic_data_format.compare("JSON") == 0){
                setIC_json(OpenWQ_json, OpenWQ_vars, OpenWQ_hostModelconfig,
                          OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_units,
                          OpenWQ_output, icmp, chemi);
            }
            else if (ic_data_format.compare("HDF5") == 0){
                setIC_h5(OpenWQ_json, OpenWQ_vars, OpenWQ_hostModelconfig,
                        OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_units,
                        OpenWQ_output, icmp, chemi);
            }
            else {
                std::string msg_string = 
                    "<OpenWQ> WARNING: Unknown data format (" 
                    + ic_data_format
                    + ") for BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp 
                    + " > INITIAL_CONDITIONS (entry skipped)";

                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                continue;
            }
        }
        /* ########################################
        // IC conditions NOT provided - set to ZERO
        ######################################## */
        else {
            std::string msg_string = 
                "<OpenWQ> IC conditions not defined for compartment:" 
                + CompName_icmp + " and chemical:" 
                + chemname + " (set to zero)";

            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }
}

/* #################################################
// Read and set IC conditions: JSON format
################################################# */
void OpenWQ_initiate::setIC_json(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    const int icmp,
    const int chemi){

    // Find compartment name
    const std::string CompName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);

    // Get chem name
    const std::string chemname = 
        (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi];

    // Get compartment dimensions
    const int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
    const int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
    const int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Validate JSON structure
    std::string errorMsgIdentifier = "Config file";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config, "BIOGEOCHEMISTRY_CONFIGURATION",
        errorMsgIdentifier, true);

    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"], CompName_icmp,
        errorMsgIdentifier, true);

    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp;
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"][CompName_icmp], "INITIAL_CONDITIONS",
        errorMsgIdentifier, true);

    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp + "INITIAL_CONDITIONS";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"][CompName_icmp]["INITIAL_CONDITIONS"], "DATA",
        errorMsgIdentifier, true);

    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp + "INITIAL_CONDITIONS > DATA";
    json ic_info_i = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"][CompName_icmp]["INITIAL_CONDITIONS"]["DATA"], chemname,
        errorMsgIdentifier, true);

    const int num_icData = ic_info_i.size();

    // Loop over all IC entries
    for (int nICchem = 0; nICchem < num_icData; nICchem++){

        // Parse spatial indices (with "ALL" support)
        auto parse_index = [&](int idx, const std::string& coord_name) -> int {
            try {
                std::string str_val = ic_info_i[std::to_string(nICchem + 1)].at(idx);
                if (str_val.compare("ALL") == 0){
                    return -1;
                } else {
                    std::string msg_string = 
                        "<OpenWQ> WARNING: Unknown entry (" + str_val
                        + ") for BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp 
                        + " > INITIAL_CONDITIONS > " + chemname 
                        + " > " + std::to_string(nICchem + 1) + " (entry skipped)";
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    return -999; // Error flag
                }
            } catch (...) {
                int val = ic_info_i[std::to_string(nICchem + 1)].at(idx);
                return val - 1; // Convert to 0-based indexing
            }
        };

        const int ic_x = parse_index(0, "x");
        const int ic_y = parse_index(1, "y");
        const int ic_z = parse_index(2, "z");

        // Skip if error in parsing
        if (ic_x == -999 || ic_y == -999 || ic_z == -999) continue;

        // Get IC value and units
        double ic_value = ic_info_i[std::to_string(nICchem + 1)].at(3);
        const std::string ic_units = ic_info_i[std::to_string(nICchem + 1)].at(4);

        // Transform units
        std::vector<std::string> units;
        std::vector<double> unit_multiplers;
        
        OpenWQ_units.Calc_Unit_Multipliers(
            OpenWQ_wqconfig, OpenWQ_output,
            unit_multiplers, ic_units, units, true);

        OpenWQ_units.Convert_Units(ic_value, unit_multiplers);

        const bool is_mass_units = (units[1].compare("EMPTY") == 0);

        /* ########################################
        // Apply IC conditions (parallelized)
        ######################################## */

        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp for collapse(3)
            for (int ix = 0; ix < nx; ix++){
                for (int iy = 0; iy < ny; iy++){
                    for (int iz = 0; iz < nz; iz++){
                        
                        if ((ic_x == ix || ic_x == -1) &&
                            (ic_y == iy || ic_y == -1) &&
                            (ic_z == iz || ic_z == -1)){

                            const double i_volume = is_mass_units ? 1.0 :
                                OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ix, iy, iz);

                            (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi)(ix, iy, iz) = 
                                ic_value * i_volume;
                        }
                    }
                }
            }
        }
    }
}

/* #################################################
// Read and set IC conditions: H5 format
################################################# */
void OpenWQ_initiate::setIC_h5(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    const int icmp,
    const int chemi){

    // Find compartment name
    const std::string CompName_icmp = 
        OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);

    // Get compartment dimensions
    const int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
    const int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
    const int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Get chem name
    const std::string chemname = 
        (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi];

    // Validate JSON structure
    std::string errorMsgIdentifier = "Config file";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config, "BIOGEOCHEMISTRY_CONFIGURATION",
        errorMsgIdentifier, true);

    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"], CompName_icmp,
        errorMsgIdentifier, true);

    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp;
    json json_ic_h5_icmp_SubStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"][CompName_icmp], "INITIAL_CONDITIONS",
        errorMsgIdentifier, true);

    // Get H5 parameters
    errorMsgIdentifier = "Config file > BIOGEOCHEMISTRY_CONFIGURATION > " + CompName_icmp + "INITIAL_CONDITIONS";
    
    std::string ic_h5_folderPath = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_ic_h5_icmp_SubStruct, "FOLDERPATH",
        errorMsgIdentifier, true);

    std::string ic_h5_timestamp = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_ic_h5_icmp_SubStruct, "TIMESTAMP",
        errorMsgIdentifier, true);

    std::string ic_h5_units = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_ic_h5_icmp_SubStruct, "UNITS",
        errorMsgIdentifier, true);

    // Replace "/" with "|" for path compatibility
    std::size_t it = ic_h5_units.find("/");
    if (it != std::string::npos){
        ic_h5_units.replace(it, 1, "|");
    }

    // Generate full filename
    std::string ic_filenamePath = ic_h5_folderPath + "/" + CompName_icmp + 
        "@" + chemname + "#" + ic_h5_units + "-main.h5";

    // Get unit conversion multipliers
    std::vector<std::string> units;
    std::vector<double> unit_multiplers;
    
    OpenWQ_units.Calc_Unit_Multipliers(
        OpenWQ_wqconfig, OpenWQ_output,
        unit_multiplers, ic_h5_units, units, true);

    // Load H5 data
    arma::mat xyzIC_h5;
    arma::mat dataIC_h5;

    xyzIC_h5.load(arma::hdf5_name(ic_filenamePath, "xyz_elements"));
    dataIC_h5.load(arma::hdf5_name(ic_filenamePath, ic_h5_timestamp));

    // Check if data was found
    if (dataIC_h5.is_empty()){
        std::string msg_string = 
            "<OpenWQ> WARNING: IC conc not found for timestamp='" + ic_h5_timestamp +
            "' in h5 file='" + ic_filenamePath +
            "'. Request ignored and IC concentrations defaulted to zero for all chemicals of "
            "compartment=" + CompName_icmp;

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        return;
    }

    const bool is_mass_units = (units[1].compare("EMPTY") == 0);
    const int num_rows = xyzIC_h5.n_rows - 1;

    // Apply IC conditions (parallelized)
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (int irow = 0; irow < num_rows; irow++){

            // Get indices (convert from 1-based to 0-based)
            const int ic_x = xyzIC_h5(irow, 0) - 1;
            const int ic_y = xyzIC_h5(irow, 1) - 1;
            const int ic_z = xyzIC_h5(irow, 2) - 1;

            // Bounds check
            if ((ic_x >= 0 && ic_x < nx) &&
                (ic_y >= 0 && ic_y < ny) &&
                (ic_z >= 0 && ic_z < nz)){

                double ic_value = dataIC_h5(irow);

                // Convert units
                OpenWQ_units.Convert_Units(ic_value, unit_multiplers);

                const double i_volume = is_mass_units ? 1.0 :
                    OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ic_x, ic_y, ic_z);

                (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi)(ic_x, ic_y, ic_z) = 
                    ic_value * i_volume;
            }
            else {
                #pragma omp critical
                {
                    std::string msg_string = 
                        "<OpenWQ> WARNING: IC h5-format data load out of boundaries. "
                        "Requested load ignored: "
                        "Compartment=" + CompName_icmp +
                        ", Chemical=" + chemname +
                        ", ix=" + std::to_string(ic_x) +
                        ", iy=" + std::to_string(ic_y) +
                        ", iz=" + std::to_string(ic_z);

                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                }
            }
        }
    }
}