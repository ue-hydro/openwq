

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


// Set up logFile
void OpenWQ_readjson::SetConfigInfo_output_logFile(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string errorMsgIdentifier;
    std::string output_format;
    std::string output_folder;
    std::string msg_string;

    // Check if OPENWQ_OUTPUT field exists
    errorMsgIdentifier = "Master file";
    OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "OPENWQ_OUTPUT",
        errorMsgIdentifier,
        true);

    errorMsgIdentifier = "Master file > OPENWQ_OUTPUT";
    output_format = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master["OPENWQ_OUTPUT"],"FORMAT",
        errorMsgIdentifier,
        true);

    output_folder = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master["OPENWQ_OUTPUT"],"RESULTS_FOLDERPATH",
        errorMsgIdentifier,
        true);

    // Create Log file name with full or relative path
    OpenWQ_wqconfig.create_LogFile_name_fullpath(output_format, output_folder);


    // ###############
    // Write header in log file

    // Separator
    msg_string.clear();
    msg_string = "\n ############################## \n";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, false, true); 

}