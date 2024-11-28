

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

#include "readjson/headerfile_RJSON.hpp"


// Set chemModule = BGC_FLEX
void OpenWQ_readjson::SetConfigInfo_TSModule_MMF_hype(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    json json_subStruct;
    std::string errorMsgIdentifier;
    std::string msg_string;
    std::string LE_method_local;    // module name
    std::string input_filepath;
    unsigned int num_entries;
    // Strings for inputs
    std::string input_direction;    // user input
    std::string input_upper_compartment;    // user input
    std::string input_lower_compartment;    // user input
    std::string data_format;
    // iteractive variables
    std::string icmp_i_name;
    unsigned int input_direction_index;
    unsigned int input_upper_compartment_index;
    unsigned int input_lower_compartment_index;
    bool input_upper_compartment_index_exist;
    bool input_lower_compartment_index_exist;
    // Other local variables
    typedef std::tuple<
        unsigned int, 
        unsigned int,
        unsigned int,      
        std::string       // index of chemical in transformation equation (needs to be here for loop reset)
        > HypeMMF_infoVect;          // Tuple with info and expression for BGC cyling

    // Get LATERAL_EXCHANGE module name
    errorMsgIdentifier = "Master file in MODULES > TRANSPORT_SEDIMENTS";

    // Get info for TSModule_MMF_hype function
    // Get number of entries
    errorMsgIdentifier = "TS_model file";
    json_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.TS_module, "CONFIGURATION",
        errorMsgIdentifier,
        true);

    errorMsgIdentifier = "TS_model file > DIRECTION";
    input_direction = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, "DIRECTION",
        errorMsgIdentifier,
        true); 

    errorMsgIdentifier = "TS_model file > UPPER_COMPARTMENT";
    input_upper_compartment = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, "UPPER_COMPARTMENT",
        errorMsgIdentifier,
        true); 

    errorMsgIdentifier = "TS_model file > LOWER_COMPARTMENT";
    input_lower_compartment = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, "LOWER_COMPARTMENT",
        errorMsgIdentifier,
        true); 

    errorMsgIdentifier = "TS_model file > DATA_FORMAT";
    data_format = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, "DATA_FORMAT",
        errorMsgIdentifier,
        true); 

    // #########################
    // Get corresponding indexes for all input values
    
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
            "<OpenWQ> WARNING: TSModule_MMF_hype - unkown 'DIRECTION' =>" + input_direction;

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file

        // Ignore entry
        return;

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
            "<OpenWQ> WARNING: TSModule_MMF_hype module - unkown 'COMPARTMET NAMES'";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);              // print in log file

        // Ignore entry
        return;

    }
        
    // Add values to tuple
    OpenWQ_wqconfig.TS_model->HypeMMF->info_vector = 
        HypeMMF_infoVect(
            input_direction_index,
            input_upper_compartment_index,
            input_lower_compartment_index,
            data_format);

    // The actual cohesion,  erodibility and sreroexp values will be processed in 
    // the function => OpenWQ_TS_model.TS_HYPE_MMF_config()

}