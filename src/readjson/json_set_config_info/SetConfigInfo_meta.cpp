

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


// Set metadata
void OpenWQ_readjson::SetConfigInfo_meta(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    std::string msg_string;

    // PROJECT_NAME
    msg_string.clear();
    msg_string = "PROJECT_NAME: ";
    try{ msg_string.append(OpenWQ_json.Master["PROJECT_NAME"]);
    }catch(...){ msg_string.append("unkown (not provided in Master file)");};
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,false, true);

    // GEOGRAPHICAL_LOCATION
    msg_string.clear();
    msg_string = "GEOGRAPHICAL_LOCATION: ";
    try{ msg_string.append(OpenWQ_json.Master["GEOGRAPHICAL_LOCATION"]);
    }catch(...){ msg_string.append("unkown (not provided in Master file)");};
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,false, true);
    
    // AUTHORS
    msg_string.clear();
    msg_string = "AUTHORS: ";
    try{ msg_string.append(OpenWQ_json.Master["AUTHORS"]);
    }catch(...){ msg_string.append("unkown (not provided in Master file)");};
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,false, true);

    // DATE
    msg_string.clear();
    msg_string = "DATE: ";
    try{ msg_string.append(OpenWQ_json.Master["DATE"]);
    }catch(...){ msg_string.append("unkown (not provided in Master file)");};
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,false, true);

    // Separator
    msg_string.clear();
    msg_string = "\n ##############################\n ";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,false, true);

}