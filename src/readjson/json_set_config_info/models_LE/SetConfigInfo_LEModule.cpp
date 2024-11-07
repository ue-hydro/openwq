

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

#include "readjson/headerfile_readjson.hpp"

// Set LE module options
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
    unsigned int num_entries;
    // Strings for inputs
    std::string input_direction;    // user input
    std::string input_upper_compartment;    // user input
    std::string input_lower_compartment;    // user input
    double input_k_val;                     // user input
    // iteractive variables
    std::string icmp_i_name;
    unsigned int input_direction_index;
    unsigned int input_upper_compartment_index;
    unsigned int input_lower_compartment_index;
    bool input_upper_compartment_index_exist;
    bool input_lower_compartment_index_exist;
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
        
    (OpenWQ_wqconfig.LE_module).append(LE_method_local);

    // Load information for the method
    if ((OpenWQ_wqconfig.LE_module).compare("OPENWQ_NATIVE_BOUNDMIX") == 0){
        
    //    // Get Erodibility coeficients for native IntMob function
    //    errorMsgIdentifier = "TD file";
    //    json_IntMov_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
    //        OpenWQ_wqconfig, OpenWQ_output,
    //        OpenWQ_json.TD_module, "INTMOB_CONFIGURATION",
    //        errorMsgIdentifier,
    //        true);

    //    errorMsgIdentifier = "Master file in MODULES > INTMOB_CONFIGURATION";
    //    json_K_Erodib_K = OpenWQ_utils.RequestJsonKeyVal_json(
    //        OpenWQ_wqconfig, OpenWQ_output,
    //        json_IntMov_subStruct, "K_VAL",
    //        errorMsgIdentifier,
    //        true);

    //    for (unsigned int cmpi=0;cmpi<OpenWQ_hostModelconfig.get_num_HydroComp();cmpi++){
    //        try{
    //            OpenWQ_wqconfig.OpenWQ_TE_native_IntMob_Erodib_K.push_back( 
    //                (double)json_K_Erodib_K.at(cmpi));
    //        }catch(...){
    //            msg_string = 
    //                "<OpenWQ> ERROR: Problem with LEB json > BOUNDMIX module > INTMOB_CONFIGURATION "
    //                " > K_val vector. It must have size "
    //                + std::to_string(OpenWQ_hostModelconfig.get_num_HydroComp())
    //                + " (number of compartments in host_model) and contain only"
    //                " double or integer entries.";
    //            // Print it (Console and/or Log file)
    //            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true,true); 
    //            // Abort program
    //            exit(EXIT_FAILURE);
    //        }
    //    }

        // Get BOUNDMIX JSON file
        input_filepath = OpenWQ_utils.RequestJsonKeyVal_str(
            OpenWQ_wqconfig, OpenWQ_output,
            OpenWQ_json.Master["MODULES"]["LATERAL_EXCHANGE"], "MODULE_CONFIG_FILEPATH",
            errorMsgIdentifier,
            true);

        // Save BOUNDMIX JSON file to OpenWQ_json.LE_module
        OpenWQ_readjson::read_JSON_2class(
            OpenWQ_wqconfig,
            OpenWQ_output,
            OpenWQ_utils,
            OpenWQ_json.LE_module,
            false,
            "",
            input_filepath);

        // Get info for BoundMix function
        // Get number of entries
        errorMsgIdentifier = "LE file";
        json_BoundMix_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            OpenWQ_json.LE_module, "CONFIGURATION",
            errorMsgIdentifier,
            true);

        num_entries = json_BoundMix_subStruct.size();

        for (unsigned int entry_i = 0; entry_i < num_entries; entry_i++){

            // Check if sub-structure exists
            errorMsgIdentifier = "LE file";
            json_BoundMix_subStruct_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
                OpenWQ_wqconfig, OpenWQ_output,
                json_BoundMix_subStruct, std::to_string(entry_i + 1),
                errorMsgIdentifier,
                true);

            // Get entries
            errorMsgIdentifier = "LE file > " + std::to_string(entry_i + 1);
            input_direction = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                json_BoundMix_subStruct_subStruct, "DIRECTION",
                errorMsgIdentifier,
                true); 

            input_upper_compartment = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                json_BoundMix_subStruct_subStruct, "UPPER_COMPARTMENT",
                errorMsgIdentifier,
                true); 

            input_lower_compartment = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                json_BoundMix_subStruct_subStruct, "LOWER_COMPARTMENT",
                errorMsgIdentifier,
                true); 

            input_k_val = OpenWQ_utils.RequestJsonKeyVal_double(
                OpenWQ_wqconfig, OpenWQ_output,
                json_BoundMix_subStruct_subStruct, "K_VAL",
                errorMsgIdentifier,
                true); 

            // Get corresponding indexes
            
            // 1) input_direction
            if (input_direction.compare("X") == 0){
                input_direction_index = 0;
            }else if (input_direction.compare("Y") == 0){
                input_direction_index = 1;
            }else if (input_direction.compare("Z") == 0){
                input_direction_index = 2;
            }else{

                // Create Message (Warning Message)
                msg_string = 
                    "<OpenWQ> WARNING: BOUNDMIX module - unkown 'UPPER_COMPARTMENT'or 'LOWER_COMPARTMENT' in entry = "
                    + std::to_string(entry_i);

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(
                    OpenWQ_wqconfig,    // for Log file name
                    msg_string,         // message
                    true,               // print in console
                    true);              // print in log file

                // Ignore entry
                continue;

            }

            // 2) upper compartment AND lower compartment index
            // Loop through native compartment names
            
            input_upper_compartment_index_exist = false; // reset to false
            input_lower_compartment_index_exist = false; // reset to false

            for (unsigned int icmp_i = 0; icmp_i < OpenWQ_hostModelconfig.get_num_HydroComp(); icmp_i++)
            {
                // Upper compartment: index
                if (input_upper_compartment.compare(OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp_i)) == 0){
                    input_upper_compartment_index = icmp_i;
                    input_upper_compartment_index_exist = true;
                }
                // Lower compartment: index
                if (input_lower_compartment.compare(OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp_i)) == 0){
                    input_lower_compartment_index = icmp_i;
                    input_lower_compartment_index_exist = true;
                }
            }

            // Create Message (Warning Message) if compartments provided don't exist
            if (input_upper_compartment_index_exist == false 
                || input_lower_compartment_index_exist == false){
                 
                msg_string = 
                    "<OpenWQ> WARNING: BOUNDMIX module - unkown 'DIRECTION' = " 
                    + input_direction
                    + "in entry = "
                    + std::to_string(entry_i);

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);              // print in log file

                // Ignore entry
                continue;

            }
                
            // Add values to tuple
            OpenWQ_wqconfig.OpenWQ_TE_native_BoundMix_info.push_back(
                std::tuple<unsigned int,unsigned int,unsigned int,double>(
                    input_direction_index,
                    input_upper_compartment_index,
                    input_lower_compartment_index,
                    input_k_val
                )
            );

        }
    }
    
}
