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

#include "models_TS/headerfile_TS.hpp"


 void OpenWQ_TS_model::TS_HYPE_MMF_config(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
            
    
    // Local variables
    std::string main_keyName;
    std::string Element_name;
    std::string Chemical_name;                  // chemical name
    std::vector<std::string> elm_list;          // model compartment list
    std::string err_text;                       // iteractive string for text to pass to error messages
    std::string Type;                           // from JSON filec (only used in SS)
    std::string ascii_FilePath;                 // additional information for ASCII data input
   

    std::vector<std::string> headerKeys;        // vector with header keys
    std::string headerWord;                     // interactive header words
    std::string elemName;                       // temporary element name
    json EWF_SS_json_sub_rowi;                  // json row of json-substructure json_subStruct 
    std::string ss_units_json_mass_base;             // units of row data (mass)
    std::vector<double> unit_multiplers;        // multiplers (numerator and denominator)
    
    std::string errorMsgIdentifier;             // error message section identifier
    std::string model_parameter;


    // Local variables
    std::string DataFormat;
    json json_PARAMETER_subStruct;

    // Data format
    DataFormat = std::get<3>(OpenWQ_wqconfig.TS_model->HypeMMF->info_vector);

    // Get PARAMETERS JSON data
    errorMsgIdentifier = "TS_model file";

    json_PARAMETER_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.TS_module, "PARAMETERS",
        errorMsgIdentifier,
        true);

    // #############################
    // Get parameters
    // #############################

    // 1) COHESION
    
    model_parameter = "COHESION";

    OpenWQ_wqconfig.TS_model->HypeMMF
        ->cohesion_entryTuple
            = OpenWQ_utils.LoadModelParameters_asTable_JSONorASCII(
            OpenWQ_wqconfig,
            OpenWQ_output,
            DataFormat,
            model_parameter,
            json_PARAMETER_subStruct,
            errorMsgIdentifier + " > " + model_parameter);

    // 2) ERODIBILITY
    
    model_parameter = "ERODIBILITY";

    OpenWQ_wqconfig.TS_model->HypeMMF
        ->erodibility_entryTuple
            = OpenWQ_utils.LoadModelParameters_asTable_JSONorASCII(
            OpenWQ_wqconfig,
            OpenWQ_output,
            DataFormat,
            model_parameter,
            json_PARAMETER_subStruct,
            errorMsgIdentifier + " > " + model_parameter);

    // 3) SREROEXP
    
    model_parameter = "SREROEXP";

    OpenWQ_wqconfig.TS_model->HypeMMF
        ->sreroexp_entryTuple
            = OpenWQ_utils.LoadModelParameters_asTable_JSONorASCII(
            OpenWQ_wqconfig,
            OpenWQ_output,
            DataFormat,
            model_parameter,
            json_PARAMETER_subStruct,
            errorMsgIdentifier + " > " + model_parameter);
                
}
