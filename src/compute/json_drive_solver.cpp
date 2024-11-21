
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

// Set solver info
void OpenWQ_readjson::SetConfigInfo_solver(
        OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string Solver_method_local;
    std::string errorMsgIdentifier;
    std::string msg_string;           // error/warning message string


    // Solver method
    errorMsgIdentifier = "Master file";
    Solver_method_local = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "SOLVER",
        errorMsgIdentifier,
        true);    

    (OpenWQ_wqconfig.SOLVER_module).append(Solver_method_local);

    // Check if SOLVER module exists
    if (OpenWQ_wqconfig.SOLVER_module != "BE" && 
        OpenWQ_wqconfig.SOLVER_module != "SUNDIALS"){

        // Create Message (Warning Message)
        msg_string = 
            "<OpenWQ> WARNING: SOLVER module - unkown (entry = "
            + Solver_method_local + ")";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true,true); 

        // Abort program
        exit(EXIT_FAILURE);
    }

}