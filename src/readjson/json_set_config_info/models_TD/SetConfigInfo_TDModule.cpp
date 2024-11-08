

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


// Set TD module options
void OpenWQ_readjson::SetConfigInfo_TDModule(
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

    // Load OPENWQ_INPUT JSON structure
    errorMsgIdentifier = "Master file";
    jsonMaster_SubStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "MODULES",
        errorMsgIdentifier,
        true);

    // check if field TRANSPORT_DISSOLVED exist
    errorMsgIdentifier = "Master file inside OPENWQ_INPUT > MODULES";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct, "TRANSPORT_DISSOLVED",
        errorMsgIdentifier,
        true);

    errorMsgIdentifier = "Master file inside OPENWQ_INPUT > MODULES > TRANSPORT_DISSOLVED";

    // check if field MODULE_NAME exist (compulsary)
    input_module_name = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct["TRANSPORT_DISSOLVED"],"MODULE_NAME",
        errorMsgIdentifier,
        true);

    // save TD module name
    (OpenWQ_wqconfig.TD_module).append(input_module_name);

    // Print it (Console and/or Log file)
    msg_string = "<OPENWQ> MODULE LOADED: TRANSPORT_DISSOLVED > " + OpenWQ_wqconfig.TD_module;
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig,    // for Log file name
        msg_string,         // message
        true,               // print in console
        true);              // print in log file

    // If MODULE_NAME not NONE, the get the MODULE_CONFIG_FILEPATH
    if (input_module_name.compare("NONE") != 0){

        input_filepath = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            jsonMaster_SubStruct["TRANSPORT_DISSOLVED"],"MODULE_CONFIG_FILEPATH",
            errorMsgIdentifier,
            true);

        OpenWQ_readjson::read_JSON_2class(
            OpenWQ_wqconfig,
            OpenWQ_output,
            OpenWQ_utils,
            OpenWQ_json.TD_module,
            false,
            "",
            input_filepath);
    }

    // Check if TD option not valid, through error
    if ((OpenWQ_wqconfig.TD_module).compare("OPENWQ_NATIVE_TD_ADVDISP") != 0 
        && (OpenWQ_wqconfig.TD_module).compare("NATIVE_TD_ADV") != 0 
        && (OpenWQ_wqconfig.TD_module).compare("NONE") != 0){

        // Create Message (Warning Message)
        msg_string = 
            "<OpenWQ> WARNING: TRANSPORT_DISSOLVED module - unkown (entry = "
            + OpenWQ_wqconfig.TD_module + ")";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true,true); 

        // Abort program
        exit(EXIT_FAILURE);

    }
}