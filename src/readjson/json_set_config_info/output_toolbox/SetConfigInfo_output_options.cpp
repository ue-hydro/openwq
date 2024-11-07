

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


// Set output options
void OpenWQ_readjson::SetConfigInfo_output_options(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units &OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    std::string output_format;                              // format of output (CSV, VTK)
    std::string output_folder;
    json json_time;
    std::vector<std::string> units;                         // units (numerator and denominator)
    std::vector<double> unit_multiplers;                    // multiplers (numerator and denominator)
    bool volume_unit_flag;                                  // flag to note if denominator is a volume (needed for calculation of concentration in output)
    std::string num_cores_input;                            // input of threads info as string
    std::string msg_string;                                 // error/warning message string
    std::string errorMsgIdentifier;
    json json_output_subStruct; 

    // Output sub-structure
    errorMsgIdentifier = "Master file";
    json_output_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master,"OPENWQ_OUTPUT",
        errorMsgIdentifier,
        true);

    // Output folder
    errorMsgIdentifier = "Master file > OPENWQ_OUTPUT";
    OpenWQ_wqconfig.set_output_dir(OpenWQ_utils.RequestJsonKeyVal_str(
                                OpenWQ_wqconfig, OpenWQ_output,
                                json_output_subStruct,"RESULTS_FOLDERPATH",
                                errorMsgIdentifier,
                                true));

    // Output Units (tuple with all the information needed)
    OpenWQ_wqconfig.set_output_units(
            OpenWQ_utils.RequestJsonKeyVal_str(
            OpenWQ_wqconfig, OpenWQ_output,
            json_output_subStruct,"UNITS",
            errorMsgIdentifier,
            true));

    // Flag for no water
    OpenWQ_wqconfig.noWaterConc = OpenWQ_utils.RequestJsonKeyVal_int(
            OpenWQ_wqconfig, OpenWQ_output,
            json_output_subStruct,"NO_WATER_CONC_FLAG",
            errorMsgIdentifier,
            true);
    
     // Convert time units from host model units to seconds (OpenWQ time units)
    // 1) Calculate unit multiplers
    volume_unit_flag = OpenWQ_units.Calc_Unit_Multipliers(
                OpenWQ_wqconfig, OpenWQ_output,
                unit_multiplers,                            // multiplers (numerator and denominator)
                OpenWQ_wqconfig.get_output_units(),  // input units
                units,
                false);              // direction of the conversion: 
                                     // to native (true) or 
                                     // from native to desired output units (false)
    OpenWQ_wqconfig.set_output_units_numerator(unit_multiplers[0]); // multipler for numerator
    OpenWQ_wqconfig.set_output_units_denominator(unit_multiplers[1]); // multiupler for denominator
    OpenWQ_wqconfig.set_output_units_concentration(volume_unit_flag); // flag if denominator is volume

    // Output format
    output_format = OpenWQ_utils.RequestJsonKeyVal_str(
            OpenWQ_wqconfig, OpenWQ_output,
            json_output_subStruct,"FORMAT",
            errorMsgIdentifier,
            true);

    // ########################################
    // Output format ######

    // Create OpenWQ_wqconfig.output_dir folder if nonexistant
    // OpenWQ_utils.check_mkdir(OpenWQ_wqconfig.get_output_dir());   // Called in set_output_type and set_output_dir

    // Create data type sub-folders if needed
    // CSV format
    if (output_format.compare("CSV") == 0){
        OpenWQ_wqconfig.set_output_type_csv();
    // HDF5 format
    }else if (output_format.compare("HDF5") == 0){
        OpenWQ_wqconfig.set_output_type_hdf5();
    } else {
        // Create Message (Error Message)
        msg_string = 
            "<OpenWQ> ERROR: Output type unkown: " + output_format;
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, false, true); 
        exit(EXIT_FAILURE);
    }
    // OpenWQ_utils.check_mkdir(OpenWQ_wqconfig.output_dir); // Called in set_output_type and set_output_dir

    // ########################################
    // Time ######
    // Get print/output timestep AND convert to seconds
    // Get tuple (value, unit)
    json_time = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        json_output_subStruct,"TIMESTEP",
        errorMsgIdentifier,
        true);
    
    OpenWQ_wqconfig.set_timestep_out(json_time.at(0));
    OpenWQ_wqconfig.set_timestep_out_unit(json_time.at(1));

    // Convert time units from host model units to seconds (OpenWQ time units)
    // 1) Calculate unit multiplers
    units.clear();                           // units (numerator and denominator)
    OpenWQ_units.Calc_Unit_Multipliers(
        OpenWQ_wqconfig, OpenWQ_output,
        unit_multiplers,                    // multiplers (numerator and denominator)
        OpenWQ_wqconfig.get_timestep_out_unit(),  // input units
        units,                              // units (numerator and denominator)
        true);                              // direction of the conversion: 
                                            // to native (true) or 
                                            // from native to desired output units (false)
    
    // Used to be a call to OpenWQ_Units::Convert_Units
    // Changed to encapsulate time within the host model config
    OpenWQ_wqconfig.convert_units_timestep_out(unit_multiplers);

}