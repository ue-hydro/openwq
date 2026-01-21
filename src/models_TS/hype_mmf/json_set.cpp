

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


// Set chemModule = BGC_FLEX
void OpenWQ_readjson::SetConfigInfo_TSModule_MMF_hype(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output,
    int sediment_cmp_index,
    int transport_cmp_index){

    // Local variables
    json json_subStruct;
    std::string jsonKey;
    std::string errorMsgIdentifier;
    std::string msg_string;
    std::string LE_method_local;    // module name
    std::string input_filepath;
    // Strings for inputs
    std::string input_direction;    // user input
    std::string erodingFlux_cmpt;            // user input
    std::string erodingInhibit_cmpt;   // user input
    std::string data_format;
    // iteractive variables
    std::string icmp_i_name;
    int input_direction_index;
    int erodingFlux_cmpt_index;
    int sedCmp_index;
    int erodingInhibit_cmpt_index;
    bool erodingFlux_cmpt_index_exist;
    bool erodingInhibit_cmpt_index_exist;
    // Other local variables
    typedef std::tuple<
        unsigned int,   // input_direction_index
        unsigned int,   // erodingFlux_cmpt_index
        unsigned int,   // sedCmp_index
        unsigned int,   // erodingInhibit_cmpt_index  
        std::string     // data format
        > HypeMMF_infoVect;          // Tuple with info and expression for BGC cyling

    // Get LATERAL_EXCHANGE module name
    errorMsgIdentifier = "Master file in MODULES > TRANSPORT_SEDIMENTS";

    // Getting eroding (sediment) compartment name
    sedCmp_index = OpenWQ_hostModelconfig.get_HydroComp_index(
        OpenWQ_wqconfig.TS_model->SedCmpt,
        errorMsgIdentifier,
        true);

    // Get info for TSModule_MMF_hype function
    // Get number of entries
    errorMsgIdentifier = "TS_model file";
    json_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.TS_module, "CONFIGURATION",
        errorMsgIdentifier,
        true);

    jsonKey = "DIRECTION";
    errorMsgIdentifier = "TS_model file >" + jsonKey;
    input_direction = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, jsonKey,
        errorMsgIdentifier,
        true); 
    
    jsonKey = "EROSION_INHIBIT_COMPARTMENT";
    errorMsgIdentifier = "TS_model file > " + jsonKey;
    erodingInhibit_cmpt = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, jsonKey,
        errorMsgIdentifier,
        true); 

    jsonKey = "DATA_FORMAT";
    errorMsgIdentifier = "TS_model file >" + jsonKey;
    data_format = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        json_subStruct, jsonKey,
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

    // Get upper and inhibiting compartments
    erodingInhibit_cmpt_index = OpenWQ_hostModelconfig.get_HydroComp_index(
        erodingInhibit_cmpt,
        errorMsgIdentifier,
        true);

    // Add values to tuple
    OpenWQ_wqconfig.TS_model->HypeMMF->info_vector = 
        HypeMMF_infoVect(
            input_direction_index,
            transport_cmp_index,
            erodingInhibit_cmpt_index,
            sediment_cmp_index,
            data_format);

    // The actual cohesion,  erodibility and sreroexp values will be processed in 
    // the function => OpenWQ_TS_model.TS_HYPE_MMF_config()

}