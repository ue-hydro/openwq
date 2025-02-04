

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

#include "models_CH/headerfile_CH.hpp"

void OpenWQ_CH_model::CH_driver_config(
        OpenWQ_json& OpenWQ_json,
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_units& OpenWQ_units,
        OpenWQ_output& OpenWQ_output){
    
    // Local variable
    std::string msg_string;

    // NATIVE Bigoeochemical model
    // Parse biogeochemical expressions (and save in global)
    if ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0)
    {
        
        bgc_flex_setBGCexpressions(
            OpenWQ_json,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_vars,
            OpenWQ_units,
            OpenWQ_output);

    }else if ((OpenWQ_wqconfig.CH_model->BGC_module).compare("PHREEQC") == 0){
        phreeqc_setup(
            OpenWQ_json,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_vars,
            OpenWQ_units,
            OpenWQ_output);
    } else{

        // Create Message
        msg_string = 
            "<OpenWQ> ERROR: No BGC_module found or unkown";

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