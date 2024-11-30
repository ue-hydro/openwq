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
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
            
    
    // Local variables    
    std::string errorMsgIdentifier;             // error message section identifier
    std::string model_parameter;
    std::string DataFormat;
    unsigned int icmp;
    json json_PARAMETER_subStruct;
    json json_PARAMETER_DEFAULTS_subStruct;

    // Data format & affected compartment (lower compartment)
    DataFormat = std::get<3>(OpenWQ_wqconfig.TS_model->HypeMMF->info_vector);
    // Get compartment index of lower compartment
    icmp = std::get<2>(OpenWQ_wqconfig.TS_model->HypeMMF->info_vector);

    // Get PARAMETERS JSON data
    errorMsgIdentifier = "TS_model file";

    json_PARAMETER_DEFAULTS_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.TS_module, "PARAMETER_DEFAULTS",
        errorMsgIdentifier,
        true);

    json_PARAMETER_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.TS_module, "PARAMETERS",
        errorMsgIdentifier,
        true);

    // #############################
    // Get parameters
    // #############################

    // 1) COHESION
    
    model_parameter = "COHESION"; // kPa

    OpenWQ_wqconfig.TS_model->HypeMMF
        ->cohesion_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // 2) ERODIBILITY
    
    model_parameter = "ERODIBILITY"; // g/J

    OpenWQ_wqconfig.TS_model->HypeMMF
        ->erodibility_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // 3) SREROEXP
    
    model_parameter = "SREROEXP"; // [-]

    OpenWQ_wqconfig.TS_model->HypeMMF
        ->sreroexp_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);
                
}
