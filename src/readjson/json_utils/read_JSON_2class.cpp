

// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
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

/* #################################################
// Read JSON file to class
################################################# */
void OpenWQ_readjson::read_JSON_2class(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_utils& OpenWQ_utils,
    json& jsondata,                         // JSON structure to save to
    const bool substruc_flag,               // Save in subfield of JSON structure? only if multiple files (e.g., source and sinks)
    const std::string JsonSubStruct_name,   // if true, name of subfield
    const std::string jsonfile)             // Name of JSON file
{ 

    // Local Variables
    std::string msg_string;             // error/warning message string
    bool readfail = false;

    /* ####################################
    // Save JSON file in JSON data structure
    #################################### */

    try
    {
        // Read json file with ifstrem
        std::ifstream f(jsonfile);

        // Raise exit flag if cannot find file
        if(f.fail()){

            // Error file
            readfail = true;

            // Create Message
            msg_string = 
                "<OpenWQ> ERROR: Could not find file: " 
                + jsonfile;
                
        }

        //* ##########################
        // If substruc_flag = true
        // save JSON inside JsonSubStruct_name (which is actually of int type) field in jsondata
        //* ####################################
        if (substruc_flag == false)
        {
            // Save ifstream to JSON data type
            jsondata = json::parse(     f,
                /* callback */          nullptr,
                /* allow exceptions */  false,
                /* skip_comments */     true
                );

            // Convert all text in JSON to upper case
            OpenWQ_readjson::ConvertJSONtext_2upperCase(
                OpenWQ_utils,
                jsondata);
        }
        else
        {
            // Save ifstream to JSON data type
            jsondata[JsonSubStruct_name] = json::parse(f,
                /* callback */                          nullptr,
                /* allow exceptions */                  false,
                /* skip_comments */                     true
                );
                    
            // Convert all text in JSON to upper case
            OpenWQ_readjson::ConvertJSONtext_2upperCase(
                OpenWQ_utils,
                jsondata[JsonSubStruct_name]);
        }
    }
    /* ####################################
    // Handle exceptions 
    #################################### */
    catch (const std::exception &e)
    {
        // Error flag
        readfail = true;

        // Create Message
        msg_string = 
            "<OpenWQ> ERROR: An exception occurred parsing JSON file: " 
            + jsonfile;

    }

    // Exist if problem with file
    if(readfail){

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
