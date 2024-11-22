
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


/* #################################################
// Compute chemical transformations
################################################# */
void OpenWQ_CH_model::CH_driver_run(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string msg_string; // error/warning message
    
     // NATIVE Bioogeochemical model
    if ((OpenWQ_wqconfig.CH->BGC_module).compare("NATIVE_BGC_FLEX") == 0){

        // Loop over number of compartments
        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){
        
            bgc_flex_transform( // calls exprtk: parsing expression
                OpenWQ_json,
                OpenWQ_vars,
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                icmp);

        }

        // Turn off printing of exception err message (only print on first step)
        OpenWQ_wqconfig.BGC_Transform_print_errmsg = false;
        OpenWQ_wqconfig.invalid_bgc_entry_errmsg = false;

    }else{

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
