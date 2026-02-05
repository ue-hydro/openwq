

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


// Parse Freundlich isotherm configuration from JSON
// JSON format:
// {
//     "MODULE_NAME": "FREUNDLICH",
//     "SOIL_PROPERTIES": {
//         "BULK_DENSITY_KG/M3": 1500.0,
//         "LAYER_THICKNESS_M": 1.0
//     },
//     "SPECIES": {
//         "species_name_1": {
//             "KFR": 0.5,
//             "NFR": 0.7,
//             "KADSDES_1/S": 0.001
//         },
//         "species_name_2": { ... }
//     }
// }
void OpenWQ_readjson::SetConfigInfo_SIModule_freundlich(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string errorMsgIdentifier;
    std::string msg_string;
    json si_json;

    errorMsgIdentifier = "SI Module (Freundlich) JSON file";

    // Get the SI module JSON data
    si_json = OpenWQ_json.SI_module;

    // ############################
    // Read soil properties
    // ############################
    json soilProps = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        si_json, "SOIL_PROPERTIES",
        errorMsgIdentifier,
        true);

    OpenWQ_wqconfig.SI_model->FREUNDLICH->bulk_density =
        OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            soilProps, "BULK_DENSITY_KG/M3",
            errorMsgIdentifier + " > SOIL_PROPERTIES",
            true);

    OpenWQ_wqconfig.SI_model->FREUNDLICH->layer_thickness =
        OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            soilProps, "LAYER_THICKNESS_M",
            errorMsgIdentifier + " > SOIL_PROPERTIES",
            true);

    // ############################
    // Read per-species parameters
    // ############################
    json speciesBlock = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        si_json, "SPECIES",
        errorMsgIdentifier,
        true);

    // Get the chemical species list for name-to-index resolution
    const unsigned int num_chem = OpenWQ_wqconfig.cached_num_chem;
    const std::vector<std::string>& chem_list = *OpenWQ_wqconfig.cached_chem_species_list_ptr;

    // Iterate over each species entry in the JSON
    unsigned int si_species_count = 0;
    for (auto it = speciesBlock.begin(); it != speciesBlock.end(); ++it) {

        std::string species_name_upper = it.key();  // Already uppercased by ConvertJSONtext_2upperCase

        // Find matching global species index
        int species_idx = -1;
        for (unsigned int c = 0; c < num_chem; c++) {
            if (chem_list[c] == species_name_upper) {
                species_idx = static_cast<int>(c);
                break;
            }
        }

        if (species_idx < 0) {
            msg_string =
                "<OpenWQ> WARNING (Freundlich SI): Species '"
                + species_name_upper
                + "' not found in chemical species list. Skipping.";
            OpenWQ_output.ConsoleLog(
                OpenWQ_wqconfig, msg_string, true, true);
            continue;
        }

        // Get this species' parameter block
        json specParams = it.value();
        std::string specErrId = errorMsgIdentifier + " > SPECIES > " + species_name_upper;

        double Kfr_val = OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            specParams, "KFR",
            specErrId, true);

        double Nfr_val = OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            specParams, "NFR",
            specErrId, true);

        double Kadsdes_val = OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            specParams, "KADSDES_1/S",
            specErrId, true);

        // Store in the FREUNDLICH data structure
        OpenWQ_wqconfig.SI_model->FREUNDLICH->species_index.push_back(
            static_cast<unsigned int>(species_idx));
        OpenWQ_wqconfig.SI_model->FREUNDLICH->species_name.push_back(species_name_upper);
        OpenWQ_wqconfig.SI_model->FREUNDLICH->Kfr.push_back(Kfr_val);
        OpenWQ_wqconfig.SI_model->FREUNDLICH->Nfr.push_back(Nfr_val);
        OpenWQ_wqconfig.SI_model->FREUNDLICH->Kadsdes.push_back(Kadsdes_val);

        si_species_count++;

        // Log
        msg_string =
            "<OpenWQ> SI (Freundlich): Species='" + species_name_upper
            + "' (idx=" + std::to_string(species_idx)
            + "), Kfr=" + std::to_string(Kfr_val)
            + ", Nfr=" + std::to_string(Nfr_val)
            + ", Kadsdes=" + std::to_string(Kadsdes_val);
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig, msg_string, true, true);
    }

    OpenWQ_wqconfig.SI_model->FREUNDLICH->num_species = si_species_count;

    // Log soil properties
    msg_string =
        "<OpenWQ> SI (Freundlich): bulk_density="
        + std::to_string(OpenWQ_wqconfig.SI_model->FREUNDLICH->bulk_density)
        + " kg/m3, layer_thickness="
        + std::to_string(OpenWQ_wqconfig.SI_model->FREUNDLICH->layer_thickness)
        + " m, num_species=" + std::to_string(si_species_count);
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig, msg_string, true, true);

}
