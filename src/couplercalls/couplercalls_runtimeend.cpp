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


#include "headerfile_couplercalls.hpp"


// ################################################################
// Calls all functions required inside time loop
// But AFTER space loop has been finalized
// ################################################################
void OpenWQ_couplercalls::RunTimeLoopEnd(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
    OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
    OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
    OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
    OpenWQ_readjson& OpenWQ_readjson,            // read json files
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
    OpenWQ_TD_model& OpenWQ_TD_model,      // transport modules
    OpenWQ_LE_model& OpenWQ_LE_model,           // LE model
    OpenWQ_models_CH& OpenWQ_models_CH,                   // biochemistry modules
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
    OpenWQ_compute& OpenWQ_compute,
    OpenWQ_output& OpenWQ_output,
    time_t simtime){


    /* ########################################
    Solver 
    (needs to be inside the time loop, but outside the space loop
    Only place where the state-variables are changed
    ######################################## */ 

    OpenWQ_compute.Numerical_Solver(
        OpenWQ_hostModelconfig,
        OpenWQ_wqconfig,
        OpenWQ_vars, 
        OpenWQ_json,
        OpenWQ_output,
        OpenWQ_models_CH);

    // ########################################
    // Output Results
    // ###########################################

    // Only print if time to print -> Needs to be adapted to host model time conventions
    // Note that OpenWQ_wqconfig.timetep_out converted to seconds in OpenWQ_readjson    

    // Print/Save results
    // Return if not time to print
    if (simtime >= OpenWQ_wqconfig.get_nexttime_out()){

        // Get compartment volumes for calculation of concentration (if requested)
        // neet to convert from mm (CRHM native units) to m3 (OpenWQ native units)

        // Call main printing function
        OpenWQ_output.writeResults(
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_compute,
            simtime);  // needs to be in seconds since 00:00 hours, Jan 1, 1970 UTC

    }

}