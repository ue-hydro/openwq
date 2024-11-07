
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


/* #################################################
 // Real all configuration files
 ################################################# */
void OpenWQ_readjson::read_all_JSON(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_units &OpenWQ_units,
    OpenWQ_utils & OpenWQ_utils,
    OpenWQ_output& OpenWQ_output)
{

    // Local
    std::string input_module_name;
    std::string input_filepath;
    std::string errorMsgIdentifier;
    std::string msg_string;

    /* ################################################
    // Read JSON files
    // Master, Config, ExtWatFlux, and SinkSource
    ################################################ */

    /* ########################
    // Read with Master file (read)
    ######################## */
    OpenWQ_readjson::read_JSON_2class(
        OpenWQ_wqconfig,
        OpenWQ_output,
        OpenWQ_utils,
        OpenWQ_json.Master,                 // JSON structure TO SAVE TO
        false,                              // Save in subfield of JSON structure? only if multiple files (e.g., source and sinks)
        "",                                 // if above true, provide name of subfield
        OpenWQ_wqconfig.get_OpenWQ_masterjson()); // Name of JSON file

    // Print it (Console and/or Log file)
    msg_string = "<OPENWQ> MASTER file loaded > " + OpenWQ_wqconfig.get_OpenWQ_masterjson();
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig,    // for Log file name
        msg_string,         // message
        true,               // print in console
        true);              // print in log file


    /* ################################################
    // Load all info from Master file
    // Process all relevant config inputs
    // Data Preparation
    // including module-specific JSON files and other config info
    ################################################ */
    OpenWQ_readjson::SetConfigInfo_driver(
        OpenWQ_json,
        OpenWQ_hostModelconfig,
        OpenWQ_wqconfig,
        OpenWQ_utils,
        OpenWQ_units,
        OpenWQ_output);

}


// ########################################
// Change JSON key
// ########################################
void OpenWQ_readjson::change_JSON_key_to_upper_case(
    OpenWQ_utils &OpenWQ_utils,
    json &object, 
    const std::string& old_key,
    std::string& new_key)
{

    try{

        // Convert to upper case
        new_key = OpenWQ_utils.ConvertStringToUpperCase(old_key);

        // get iterator to old key; TODO: error handling if key is not present
        json::iterator it = object.find(old_key);

        // create null value for new key and swap value from old key
        std::swap(object[new_key], it.value());

    }catch(const std::exception &e){}

}


// ########################################
// Extract model configuration info from JSON files
// ########################################
void OpenWQ_readjson::SetConfigInfo_driver(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_hostModelconfig & OpenWQ_hostModelconfig,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units &OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
    

    // #############################
    // Set General info

    // Set logFile
    // First thing because the log file needs to be initiated to 
    // log info about config
    SetConfigInfo_output_logFile(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    // Set metadata
    SetConfigInfo_meta(
        OpenWQ_json, OpenWQ_wqconfig,OpenWQ_output);

    // Set computation settings
    SetConfigInfo_compute(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    // Set solver settings
    SetConfigInfo_solver(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    // #############################
    // Set OPENWQ_INPUT info

    // Set general config data
    SetConfigInfo_INPUT_general(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    // Set external fluxes
    SetConfigInfo_INPUT_EF(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    // set sink and source
    SetConfigInfo_INPUT_SS(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);


    // #############################
    // Set modules

    // Set chem module settings
    SetConfigInfo_chemModule(  
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_output);

    // Set TD module settings (transport dissolved)
    SetConfigInfo_TDModule(
        OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, 
        OpenWQ_utils, OpenWQ_output);

    // Set LE module settings (lateral exchange)
    SetConfigInfo_LEModule(
        OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, 
        OpenWQ_utils, OpenWQ_output);

    // Set ST module settings (sediment transport)
    //SetConfigInfo_STModule(
    //    OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, 
    //    OpenWQ_utils, OpenWQ_output);

    // Set SI module settings (sorption isotherm)
    //SetConfigInfo_SIModule(
    //    OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, 
    //    OpenWQ_utils, OpenWQ_output);

    // #############################
    // Set output options

    SetConfigInfo_output_driver(
        OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig, 
        OpenWQ_utils, OpenWQ_units, OpenWQ_output);

}
