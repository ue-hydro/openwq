
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

// Set chem info
void OpenWQ_readjson::SetConfigInfo_CHModule(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    std::string errorMsgIdentifier;
    json jsonMaster_SubStruct;
    json BGCjson_subStruct;
    json BGCjson_ChemList;
    json BGCjson_mobileSpecies;
    std::string msg_string;           // error/warning message string
    std::string input_module_name;
    std::string input_filepath;

    // check if field MODULES exist
    // this field needs to exist
    errorMsgIdentifier = "Master file inside OPENWQ_INPUT";
    jsonMaster_SubStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "MODULES",
        errorMsgIdentifier,
        true);

    // check if field BIOGEOCHEMISTRY exist 
    errorMsgIdentifier = "Master file inside OPENWQ_INPUT > MODULES";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct, "BIOGEOCHEMISTRY",
        errorMsgIdentifier,
        true);

    errorMsgIdentifier = "Master file inside OPENWQ_INPUT > MODULES > BIOGEOCHEMISTRY";

    // check if field MODULE_NAME exist (compulsary)
    // if not BGQ is desired, then choose MODULE_NAME = NONE
    input_module_name = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct["BIOGEOCHEMISTRY"],"MODULE_NAME",
        errorMsgIdentifier,
        true);

    // save BGC module name
    (OpenWQ_wqconfig.BGC_module).append(input_module_name);

    // Print it (Console and/or Log file)
    msg_string = "<OPENWQ> MODULE LOADED: BIOGEOCHEMISTRY > " 
                + OpenWQ_wqconfig.BGC_module;
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig,    // for Log file name
        msg_string,         // message
        true,               // print in console
        true);              // print in log file

    // Check if CH option not valid, through error
    if ((OpenWQ_wqconfig.BGC_module).compare("NONE") == 0){

        // Create Message (Warning Message)
        msg_string = 
            "<OpenWQ> WARNING: The BIOGEOCHEMISTRY module was set to NONE. "
            "That is NOT ALLOWED because chemical constitutes are not define, so not even the TRANSPORT routines could work.";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true,true); 

        // Abort program
        exit(EXIT_FAILURE);

    }

    // Check if CH option not valid, through error
    if ((OpenWQ_wqconfig.BGC_module).compare("NATIVE_BGC_FLEX") != 0 
        && (OpenWQ_wqconfig.BGC_module).compare("NONE") != 0){

        // Create Message (Warning Message)
        msg_string = 
            "<OpenWQ> WARNING: BIOGEOCHEMISTRY module - unkown (entry = "
            + OpenWQ_wqconfig.BGC_module + ")";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true,true); 

        // Abort program
        exit(EXIT_FAILURE);

    }

    // Get BGC module and extract JSON to class
    input_filepath = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            jsonMaster_SubStruct["BIOGEOCHEMISTRY"],"MODULE_CONFIG_FILEPATH",
            errorMsgIdentifier,
            true);

    OpenWQ_readjson::read_JSON_2class(
        OpenWQ_wqconfig,
        OpenWQ_output,
        OpenWQ_utils,
        OpenWQ_json.BGC_module,
        false,
        "",
        input_filepath);


    // Load information fo the method
    // Native module
    if ((OpenWQ_wqconfig.BGC_module).compare("NATIVE_BGC_FLEX") == 0){
        
        SetConfigInfo_CHModule_BGC_FLEX(  
            OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    }

}