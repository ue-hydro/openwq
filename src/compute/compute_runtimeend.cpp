
// Copyright 2020, Diogo Costa (diogo.costa@uevora.pt)
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


#include "compute/headerfile_compute.hpp"
#include "models_chem/openwq_chem_headerfile.hpp"

/* #################################################
// General numerical solver
################################################# */

void OpenWQ_solver::Numerical_Solver(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_chem& OpenWQ_chem){
    
    // Local variables
    std::string msg_string; // error/warning message
    

    if ((OpenWQ_wqconfig.SOLVER_module).compare("SUNDIALS") == 0){

        Solve_with_CVode(
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_vars,
            OpenWQ_json,
            OpenWQ_output,
            OpenWQ_chem);

    } else if ((OpenWQ_wqconfig.SOLVER_module).compare("BE") == 0) {
        
        Solve_with_BE(
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_vars,
            OpenWQ_json,
            OpenWQ_output,
            OpenWQ_chem);

    } else {
        // Create Message
        msg_string = 
            "<OpenWQ> ERROR: Unknown Solver choice";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file
        
        // Abort (Fatal error)
        exit(EXIT_FAILURE);

    }
}
