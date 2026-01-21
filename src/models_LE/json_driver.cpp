

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

#include "readjson/headerfile_RJSON.hpp"

// Set LE_model module options
void OpenWQ_readjson::SetConfigInfo_LEModule(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_hostModelconfig & OpenWQ_hostModelconfig,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    std::string errorMsgIdentifier;
    std::string msg_string;
    std::string LE_method_local;    // module name
    std::string input_filepath;
    // Strings for inputs
    std::string input_direction;    // user input
    std::string input_upper_compartment;    // user input
    std::string input_lower_compartment;    // user input
    // iteractive variables
    std::string icmp_i_name;
    json json_K_Erodib_K;
    json json_IntMov_subStruct;
    json json_BoundMix_subStruct;
    json json_BoundMix_subStruct_subStruct;

    // check if field MODULES exist
    errorMsgIdentifier = "Master file";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "MODULES",
        errorMsgIdentifier,
        true);

    // check if field MODULES > LATERAL_EXCHANGE exist
    errorMsgIdentifier = "Master file in MODULES";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master["MODULES"], "LATERAL_EXCHANGE",
        errorMsgIdentifier,
        true);

    // Get LATERAL_EXCHANGE module name
    errorMsgIdentifier = "Master file in MODULES > LATERAL_EXCHANGE";
    LE_method_local = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master["MODULES"]["LATERAL_EXCHANGE"], "MODULE_NAME",
        errorMsgIdentifier,
        true);
        
    (OpenWQ_wqconfig.LE_model->LE_module).append(LE_method_local);

    // Check if TD_model option not valid, through error
    if ((OpenWQ_wqconfig.LE_model->LE_module).compare("NATIVE_LE_BOUNDMIX") != 0 
        && (OpenWQ_wqconfig.LE_model->LE_module).compare("NONE") != 0){

        // Create Message (Warning Message)
        msg_string = 
            "<OpenWQ> WARNING: BOUNDMIX module - unkown (entry = "
            + OpenWQ_wqconfig.LE_model->LE_module + ")";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true,true); 

        // Abort program
        exit(EXIT_FAILURE);
    }

    // if not NONE, load module configuration file
    if ((OpenWQ_wqconfig.LE_model->LE_module).compare("NONE") != 0){

        // Get configuration file JSON file
        input_filepath = OpenWQ_utils.RequestJsonKeyVal_str(
            OpenWQ_wqconfig, OpenWQ_output,
            OpenWQ_json.Master["MODULES"]["LATERAL_EXCHANGE"], "MODULE_CONFIG_FILEPATH",
            errorMsgIdentifier,
            true);

        // Save configuration file for OpenWQ_json.LE_module
        OpenWQ_readjson::read_JSON_2class(
            OpenWQ_wqconfig,
            OpenWQ_output,
            OpenWQ_utils,
            OpenWQ_json.LE_module,
            false,
            "",
            input_filepath);

        // Load information for the method
        if ((OpenWQ_wqconfig.LE_model->LE_module).compare("NATIVE_LE_BOUNDMIX") == 0){
            
            SetConfigInfo_LEModule_BOUNDMIX(  
                OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);
            
        };
    }
    
}

 