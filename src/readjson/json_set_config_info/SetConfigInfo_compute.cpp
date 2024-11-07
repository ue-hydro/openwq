
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

#include "readjson/headerfile_readjson.hpp"

// Set computation options
void OpenWQ_readjson::SetConfigInfo_compute(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    std::string num_cores_input;                            // input of threads info as string
    bool threads_input_err_flag = false;                    // flag for error in threads input (initiate a false, no error)
    std::string msg_string;                                 // error/warning message string
    std::string errorMsgIdentifier;
    json json_computeSettings;

    // Get json-subStructurd related to computational settings
    errorMsgIdentifier = "Master file";
    json_computeSettings = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master, "COMPUTATIONAL_SETTINGS",
        errorMsgIdentifier,
        true);

    // Run model: check if DEBUG model has been requested
    // The model will print all the derivatives
    errorMsgIdentifier = "Master file in COMPUTATIONAL_SETTINGS";
    OpenWQ_wqconfig.debug_mode = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        json_computeSettings, "RUN_MODE_DEBUG",
        errorMsgIdentifier,
        true);

    // Get number of threads in the system
    OpenWQ_wqconfig.set_num_threads_system(sysconf(_SC_NPROCESSORS_ONLN)); 
    
    // Get number of threads requested
    // It can be provided as an int or "all" (string)
    errorMsgIdentifier = "Master file in COMPUTATIONAL_SETTINGS";
    try{
        // Test if provided as integer with number of cores
        // Convert cast to int (if double or flow)
        OpenWQ_wqconfig.set_num_threads_requested((int) json_computeSettings["USE_NUM_THREADS"]);

    }catch(...){
        
        // Check if equal to "ALL"
        try{
            // If not integer, check if string
            num_cores_input = json_computeSettings["USE_NUM_THREADS"];
        
            // If string check if it is a string equal to "ALL"
            if (num_cores_input.compare("ALL") == 0){

                OpenWQ_wqconfig.set_threads_requested_to_system();

            }else{
                // Set error flag = true
                threads_input_err_flag = true;
            }
        }catch(...){
             // Set error flag = true
            threads_input_err_flag = true;
        }

        // Error handling 
        // if not "ALL" get integer. Throw error if not integer either
        if (threads_input_err_flag == true){

            // Create Message
            msg_string =
                 "<OpenWQ> WARNNING: Unrecognized input for COMPUTATIONAL_SETTINGS > USE_NUM_THREADS: '"
                + num_cores_input
                + "' (Only allowed \"all\" or integer). USE_NUM_THREADS has been defaulted to "
                + std::to_string(OpenWQ_wqconfig.get_num_threads_default()) + ".";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, false, true);

            // Default number of threads
            OpenWQ_wqconfig.set_threads_requested_to_default();

        }
    }
    
    // Check if number of threads requested is higher than existing in the system
    // if yes, throw warning message and limited the request
    if (OpenWQ_wqconfig.is_num_threads_requested_valid() == false){

        // Create Message
        msg_string = OpenWQ_wqconfig.get_num_threads_warning();
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,
                                false, true);
        
        // Overwrite OpenWQ_wqconfig.num_threads_requested to max number in the system
        OpenWQ_wqconfig.set_threads_requested_to_system();
        
    }else{

        // Create Message
        msg_string = OpenWQ_wqconfig.get_num_threads_info();

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string,false, true);

    }

}