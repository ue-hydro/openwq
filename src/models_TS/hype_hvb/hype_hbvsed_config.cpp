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

#include "models_TS/headerfile_TS.hpp"


 void OpenWQ_TS_model::TS_HYPE_HBVSED_config(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
            
    
    // Local variables    
    std::string errorMsgIdentifier;             // error message section identifier
    std::string model_parameter;
    std::string DataFormat;
    unsigned int icmp_SedCmpt;
    json json_PARAMETER_subStruct;
    json json_PARAMETER_DEFAULTS_subStruct;


    // Get compartment index of lower compartment
    icmp_SedCmpt = OpenWQ_hostModelconfig.get_HydroComp_index(
        OpenWQ_wqconfig.TS_model->SedCmpt,
        "<OpenWQ> ERROR: TS_model > Unkown sediment compartment name: " 
        + OpenWQ_wqconfig.TS_model->SedCmpt,
        true); 
       
    // Data format & affected compartment (lower compartment)
    DataFormat = OpenWQ_wqconfig.TS_model->HypeHVB->get_data_format();    

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
        false);  // Not compulsory: empty PARAMETERS means use PARAMETER_DEFAULTS for all cells

    // #############################
    // Get parameters
    // #############################

    // #############################
    // 1) SOIL_EROSION_FACTOR_LAND_DEPENDENCE
    
    model_parameter = "SOIL_EROSION_FACTOR_LAND_DEPENDENCE"; //  lusepar       ! soil erosion factor (land use dependence)

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->lusepar_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // #############################
    // 2) SOIL_EROSION_FACTOR_SOIL_DEPENDENCE
    
    model_parameter = "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE"; // soil erosion factor (soil dependence)

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->soilpar_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // #############################
    // 3) SLOPE_EROSION_FACTOR_EXPONENT
    
    model_parameter = "SLOPE_EROSION_FACTOR_EXPONENT"; // slopepar      ! slope erosion factor (exponent) 

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->slopepar_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // #############################
    // 4) PRECIP_EROSION_FACTOR_EXPONENT
    
    model_parameter = "PRECIP_EROSION_FACTOR_EXPONENT"; // precexppar    ! erosion precipitation dependence factor (exponent)

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->precexppar_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // #############################
    // 5) PARAM_SCALING_EROSION_INDEX
    
    model_parameter = "PARAM_SCALING_EROSION_INDEX"; // eroindexpar   ! model parameter for scaling of erosion index

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->eroindexpar_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);
    
    // #############################
    // 6) SLOPE
    
    model_parameter = "SLOPE"; // basin slope 

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->slope_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);
    
    // #############################
    // 7) EROSION_INDEX
    
    model_parameter = "EROSION_INDEX"; // erosion index

    OpenWQ_wqconfig.TS_model->HypeHVB
        ->eroindex_entryArmaCube
            = OpenWQ_utils.
            LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
                OpenWQ_hostModelconfig,
                OpenWQ_wqconfig,
                OpenWQ_output,
                DataFormat,
                model_parameter,
                icmp_SedCmpt,
                json_PARAMETER_subStruct,
                json_PARAMETER_DEFAULTS_subStruct,
                errorMsgIdentifier + " > " + model_parameter);

    // #############################
    // 8) MONTHLY_EROSION_FACTOR (optional, 12 values for Jan-Dec)
    // In HYPE: erodmonth = 1.0 + monthpar(m_erodmon, current_month)
    // This is a top-level key in the TS module JSON (not inside PARAMETERS)

    model_parameter = "MONTHLY_EROSION_FACTOR";

    if (OpenWQ_json.TS_module.contains(model_parameter)){

        json json_monthpar = OpenWQ_json.TS_module[model_parameter];

        if (json_monthpar.is_array() && json_monthpar.size() == 12){
            for (int m = 0; m < 12; m++){
                OpenWQ_wqconfig.TS_model->HypeHVB->monthpar[m] = json_monthpar[m].get<double>();
            }
        } else {

            // Create Message (Warning Message)
            std::string msg_string =
                "<OpenWQ> WARNING: TS_model -> HYPE_HBVSED - \"MONTHLY_EROSION_FACTOR\" "
                "must be an array of 12 values (Jan-Dec). Using defaults (all zeros).";

            // Print it (Console and/or Log file)
            OpenWQ_output.ConsoleLog(
                OpenWQ_wqconfig,    // for Log file name
                msg_string,         // message
                true,               // print in console
                true);              // print in log file
        }
    }
    // If not present, defaults (all zeros) are used => erodmonth = 1.0

}
