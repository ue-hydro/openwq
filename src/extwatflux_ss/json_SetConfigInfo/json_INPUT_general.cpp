
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

// Set general config info
void OpenWQ_readjson::SetConfigInfo_INPUT_general(
        OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    json jsonMaster_SubStruct;
    std::string errorMsgIdentifier;
    std::string msg_string;           // error/warning message string
    std::string input_filepath;

    // Get OPENWQ_INPUT json substructure (main configuration)
    errorMsgIdentifier = "Master file";
    jsonMaster_SubStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "OPENWQ_INPUT",
        errorMsgIdentifier,
        true);

    // Get configuration json-file filepath
    errorMsgIdentifier = "Master file inside OPENWQ_INPUT";
    input_filepath = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        jsonMaster_SubStruct, "CONFIG_FILEPATH",
        errorMsgIdentifier,
        true);

    // Save CONFIG_FILEPATH json file and save it to OpenWQ_json.Config
    OpenWQ_readjson::read_JSON_2class(
        OpenWQ_wqconfig,
        OpenWQ_output,
        OpenWQ_utils,
        OpenWQ_json.Config,
        false,
        "",
        input_filepath);

    // Print it (Console and/or Log file)
    msg_string = "<OPENWQ> CONFIG file loaded > " + input_filepath;
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig,    // for Log file name
        msg_string,         // message
        true,               // print in console
        true);              // print in log file

}
