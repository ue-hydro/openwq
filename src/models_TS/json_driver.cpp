
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

// Set LE_model module options
void OpenWQ_readjson::SetConfigInfo_TSModule(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_hostModelconfig & OpenWQ_hostModelconfig,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    json jsonMaster_SubStruct;
    std::string errorMsgIdentifier;
    std::string msg_string;
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

    // check if field TRANSPORT_SEDIMENTS exist
    errorMsgIdentifier = "Master file inside OPENWQ_INPUT > MODULES";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct, "TRANSPORT_SEDIMENTS",
        errorMsgIdentifier,
        true);

    errorMsgIdentifier = "Master file inside OPENWQ_INPUT > MODULES > TRANSPORT_SEDIMENTS";

    // check if field MODULE_NAME exist (compulsary)
    input_module_name = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct["TRANSPORT_SEDIMENTS"],"MODULE_NAME",
        errorMsgIdentifier,
        true);

    // save TS_model module name
    (OpenWQ_wqconfig.TS_model->TS_module).append(input_module_name);

    // Print it (Console and/or Log file)
    msg_string = "<OPENWQ> MODULE LOADED: TRANSPORT_SEDIMENTS > " + OpenWQ_wqconfig.TS_model->TS_module;
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig,    // for Log file name
        msg_string,         // message
        true,               // print in console
        true);              // print in log file

    // If MODULE_NAME not NONE, the get the MODULE_CONFIG_FILEPATH
    if ((OpenWQ_wqconfig.TS_model->TS_module).compare("NONE") != 0){

        input_filepath = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            jsonMaster_SubStruct["TRANSPORT_SEDIMENTS"],"MODULE_CONFIG_FILEPATH",
            errorMsgIdentifier,
            true);

        read_JSON_2class(
            OpenWQ_wqconfig,
            OpenWQ_output,
            OpenWQ_utils,
            OpenWQ_json.TS_module,
            false,
            "",
            input_filepath);
            
    }

    // Load information fo the TS_model model selected
    if ((OpenWQ_wqconfig.TS_model->TS_module).compare("HYPE_MMF") == 0){
        
        SetConfigInfo_TSModule_MMF_hype( 
            OpenWQ_hostModelconfig,
            OpenWQ_json, 
            OpenWQ_wqconfig, 
            OpenWQ_utils, 
            OpenWQ_output);

    }else if ((OpenWQ_wqconfig.TS_model->TS_module).compare("HYPE_HBVSED") == 0){
        
        SetConfigInfo_TSModule_HBVsed_hype(  
            OpenWQ_json, 
            OpenWQ_wqconfig, 
            OpenWQ_utils, 
            OpenWQ_output);

    }


}