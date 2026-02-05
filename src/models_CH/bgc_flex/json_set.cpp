

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


// Set chemModule = BGC_FLEX
void OpenWQ_readjson::SetConfigInfo_CHModule_BGC_FLEX(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string errorMsgIdentifier;
    json BGCjson_subStruct;
    json BGCjson_ChemList;
    json BGCjson_mobileSpecies;
    std::string msg_string;           // error/warning message string

    // Check if BGQ json has CHEMICAL_SPECIES key
    errorMsgIdentifier = "BGQ file";
    BGCjson_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.BGC_module, "CHEMICAL_SPECIES",
        errorMsgIdentifier,
        true); 

    // Check if BGQ json has CHEMICAL_SPECIES key
    errorMsgIdentifier = "BGQ file in CHEMICAL_SPECIES";
    BGCjson_ChemList = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        BGCjson_subStruct, "LIST",
        errorMsgIdentifier,
        true); 

    // Get number of chemical species from BGC_json
    (OpenWQ_wqconfig.CH_model->NativeFlex->num_chem) = BGCjson_ChemList.size();

    // Get chemical species list from BGC_json
    // (must be done before mobile species resolution)
    for (unsigned int chemi = 0; chemi < (OpenWQ_wqconfig.CH_model->NativeFlex->num_chem); chemi++)
    {
        (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list).push_back(
            BGCjson_ChemList[std::to_string(chemi + 1)]);
    }

    // Get mobile species (now specified by name instead of index)
    // Resolve each name to its 0-based index in chem_species_list
    errorMsgIdentifier = "BGQ file in CHEMICAL_SPECIES";
    BGCjson_mobileSpecies = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        BGCjson_subStruct, "BGC_GENERAL_MOBILE_SPECIES",
        errorMsgIdentifier,
        false);  // no abort

    for (unsigned int chemi = 0; chemi < BGCjson_mobileSpecies.size(); chemi++){

        std::string mobile_name = BGCjson_mobileSpecies.at(chemi).get<std::string>();
        bool found = false;

        for (unsigned int si = 0; si < OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list.size(); si++){
            if (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list[si] == mobile_name){
                OpenWQ_wqconfig.CH_model->NativeFlex->mobile_species.push_back((int)si);
                found = true;
                break;
            }
        }

        if (!found){
            msg_string =
                "<OpenWQ> WARNING: BGC_GENERAL_MOBILE_SPECIES name '"
                + mobile_name
                + "' not found in CHEMICAL_SPECIES LIST. Skipping.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }

    }

}
