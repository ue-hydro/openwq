

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


// Parse Langmuir isotherm configuration from JSON
// JSON format:
// {
//     "MODULE_NAME": "LANGMUIR",
//     "SOIL_PROPERTIES": {
//         "BULK_DENSITY_KG/M3": 1500.0,
//         "LAYER_THICKNESS_M": 1.0
//     },
//     "SPECIES": {
//         "species_name_1": {
//             "QMAX_MG/KG": 100.0,
//             "KL_L/MG": 0.01,
//             "KADSDES_1/S": 0.001
//         },
//         "species_name_2": { ... }
//     }
// }
void OpenWQ_readjson::SetConfigInfo_SIModule_langmuir(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string errorMsgIdentifier;
    std::string msg_string;
    json si_json;

    errorMsgIdentifier = "SI Module (Langmuir) JSON file";

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

    OpenWQ_wqconfig.SI_model->LANGMUIR->bulk_density =
        OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            soilProps, "BULK_DENSITY_KG/M3",
            errorMsgIdentifier + " > SOIL_PROPERTIES",
            true);

    OpenWQ_wqconfig.SI_model->LANGMUIR->layer_thickness =
        OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output,
            soilProps, "LAYER_THICKNESS_M",
            errorMsgIdentifier + " > SOIL_PROPERTIES",
            true);

    // ############################
    // Read per-species parameters from the BGC template's
    //   CHEMICAL_SPECIES > BGC_GENERAL_SORBABLE_SPECIES
    // (NOT the SI JSON — that only carries SOIL_PROPERTIES + the compartment).
    // Each entry:  "<dissolved>": { "INTO": "<particulate>",
    //    "LANGMUIR": { "QMAX_MG/KG", "KL_L/MG", "KADSDES_1/S" },
    //    "FREUNDLICH": { ... } }
    // ############################
    errorMsgIdentifier = "BGC module > CHEMICAL_SPECIES > BGC_GENERAL_SORBABLE_SPECIES";

    json chemSpeciesBlock = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.BGC_module, "CHEMICAL_SPECIES",
        "BGC module", false);

    json sorbBlock = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        chemSpeciesBlock, "BGC_GENERAL_SORBABLE_SPECIES",
        errorMsgIdentifier, false);

    // Chemical species list for name-to-index resolution
    const unsigned int num_chem = OpenWQ_wqconfig.cached_num_chem;
    const std::vector<std::string>& chem_list = *OpenWQ_wqconfig.cached_chem_species_list_ptr;

    auto resolve_idx = [&](const std::string& nm) -> int {
        for (unsigned int c = 0; c < num_chem; c++)
            if (chem_list[c] == nm) return static_cast<int>(c);
        return -1;
    };

    // Iterate over each sorbable species entry
    unsigned int si_species_count = 0;
    for (auto it = sorbBlock.begin(); it != sorbBlock.end(); ++it) {

        std::string species_name_upper = it.key();  // dissolved species (uppercased)
        json specObj = it.value();

        // Only configure species that declare a LANGMUIR parameter block
        if (!specObj.contains("LANGMUIR")) continue;

        int species_idx = resolve_idx(species_name_upper);
        if (species_idx < 0) {
            msg_string =
                "<OpenWQ> WARNING (Langmuir SI): sorbable species '"
                + species_name_upper
                + "' not found in CHEMICAL_SPECIES LIST. Skipping.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            continue;
        }

        std::string specErrId = errorMsgIdentifier + " > " + species_name_upper;

        // Particulate ("INTO") phase that receives the sorbed mass
        std::string into_name = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output, specObj, "INTO", specErrId, true);
        int into_idx = resolve_idx(into_name);
        if (into_idx < 0) {
            msg_string =
                "<OpenWQ> WARNING (Langmuir SI): 'INTO' particulate species '"
                + into_name + "' (for '" + species_name_upper
                + "') not found in CHEMICAL_SPECIES LIST. Skipping.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            continue;
        }

        // Langmuir parameter block
        json langParams = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output, specObj, "LANGMUIR", specErrId, true);
        std::string pErr = specErrId + " > LANGMUIR";

        double qmax_val = OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output, langParams, "QMAX_MG/KG", pErr, true);
        double KL_val = OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output, langParams, "KL_L/MG", pErr, true);
        double Kadsdes_val = OpenWQ_utils.RequestJsonKeyVal_double(
            OpenWQ_wqconfig, OpenWQ_output, langParams, "KADSDES_1/S", pErr, true);

        // Store in the LANGMUIR data structure
        OpenWQ_wqconfig.SI_model->LANGMUIR->species_index.push_back(
            static_cast<unsigned int>(species_idx));
        OpenWQ_wqconfig.SI_model->LANGMUIR->species_name.push_back(species_name_upper);
        OpenWQ_wqconfig.SI_model->LANGMUIR->into_index.push_back(
            static_cast<unsigned int>(into_idx));
        OpenWQ_wqconfig.SI_model->LANGMUIR->into_name.push_back(into_name);
        OpenWQ_wqconfig.SI_model->LANGMUIR->qmax.push_back(qmax_val);
        OpenWQ_wqconfig.SI_model->LANGMUIR->KL.push_back(KL_val);
        OpenWQ_wqconfig.SI_model->LANGMUIR->Kadsdes.push_back(Kadsdes_val);

        si_species_count++;

        msg_string =
            "<OpenWQ> SI (Langmuir): '" + species_name_upper
            + "' -> '" + into_name + "' (idx " + std::to_string(species_idx)
            + "->" + std::to_string(into_idx)
            + "), qmax=" + std::to_string(qmax_val)
            + " mg/kg, KL=" + std::to_string(KL_val)
            + " L/mg, Kadsdes=" + std::to_string(Kadsdes_val) + " 1/s";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    OpenWQ_wqconfig.SI_model->LANGMUIR->num_species = si_species_count;

    // Log soil properties
    msg_string =
        "<OpenWQ> SI (Langmuir): bulk_density="
        + std::to_string(OpenWQ_wqconfig.SI_model->LANGMUIR->bulk_density)
        + " kg/m3, layer_thickness="
        + std::to_string(OpenWQ_wqconfig.SI_model->LANGMUIR->layer_thickness)
        + " m, num_species=" + std::to_string(si_species_count);
    OpenWQ_output.ConsoleLog(
        OpenWQ_wqconfig, msg_string, true, true);

}
