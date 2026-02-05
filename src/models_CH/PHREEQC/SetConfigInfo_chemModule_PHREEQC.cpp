

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
#include "PhreeqcRM.h"
#include <memory>

// Set chemModule = PHREEQC
void OpenWQ_readjson::SetConfigInfo_CHModule_PHREEQC(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){

    // Local variables
    std::string errorMsgIdentifier;
    json BGCjson_mobileSpecies;
    std::string PHREEQC_filename;
    std::string database_filename;
    std::string msg_string;           // error/warning message string
    std::vector<std::string> component_names;

    // Check if BGQ json has FILENAME
    errorMsgIdentifier = "BGQ file";
    PHREEQC_filename = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.BGC_module, "FILEPATH",
        errorMsgIdentifier,
        true); 

    errorMsgIdentifier = "database in BGQ file";
    database_filename = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.BGC_module, "DATABASE",
        errorMsgIdentifier,
        true); 

    

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm = std::unique_ptr<PhreeqcRM>(new PhreeqcRM(1, 1));
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->LoadDatabase(database_filename);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunFile(true,true, true, PHREEQC_filename);

    // Discover chemical species via FindComponents
    // FindComponents() returns int: positive = component count, negative = error code.
    // Even on error, GetComponents() may still return valid names from the
    // internal components vector, so we use GetComponents().size() as the
    // authoritative count.
    int find_result = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->FindComponents();

    if (find_result < 0) {
        msg_string =
            "<OpenWQ> WARNING: PHREEQC FindComponents() returned error code "
            + std::to_string(find_result)
            + ". Will attempt to use GetComponents() instead.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    component_names = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponents();

    // Use GetComponents() size as the authoritative component count
    OpenWQ_wqconfig.CH_model->PHREEQC->num_chem = (unsigned int)component_names.size();

    // Validate that components were found
    if (OpenWQ_wqconfig.CH_model->PHREEQC->num_chem == 0) {
        msg_string =
            "<OpenWQ> ERROR: PHREEQC found 0 components. "
            "Check that the PHREEQC input file (" + PHREEQC_filename + ") "
            "and database (" + database_filename + ") are valid.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        exit(EXIT_FAILURE);
    }

    // Log the number and names of all discovered components
    msg_string =
        "<OpenWQ> PHREEQC: Found "
        + std::to_string(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem)
        + " components:";
    for (unsigned int i = 0; i < component_names.size(); i++) {
        msg_string += " [" + std::to_string(i + 1) + "]=" + component_names[i];
    }
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Get chemical species list from component names
    // (must be done before mobile species resolution)
    for (unsigned int chemi = 0; chemi < (OpenWQ_wqconfig.CH_model->PHREEQC->num_chem); chemi++)
    {
        (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list).push_back(
            component_names[chemi]);
    }

    // Get mobile species (now specified by name instead of index)
    // Resolve each name to its 0-based index in chem_species_list
    // Note: For PHREEQC, BGC_GENERAL_MOBILE_SPECIES is a top-level key
    // in the BGC module JSON (unlike NATIVE_BGC_FLEX where it's nested
    // under CHEMICAL_SPECIES)
    errorMsgIdentifier = "BGQ file";
    BGCjson_mobileSpecies = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.BGC_module, "BGC_GENERAL_MOBILE_SPECIES",
        errorMsgIdentifier,
        false);  // no abort

    if (BGCjson_mobileSpecies.empty()) {
        msg_string =
            "<OpenWQ> WARNING: PHREEQC BGC_GENERAL_MOBILE_SPECIES is empty or missing. "
            "No species will be transported.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    for (unsigned int chemi = 0; chemi < BGCjson_mobileSpecies.size(); chemi++){

        std::string mobile_name = BGCjson_mobileSpecies.at(chemi).get<std::string>();
        bool found = false;

        for (unsigned int si = 0; si < OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list.size(); si++){
            if (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list[si] == mobile_name){
                OpenWQ_wqconfig.CH_model->PHREEQC->mobile_species.push_back((int)si);
                found = true;
                break;
            }
        }

        if (!found){
            msg_string =
                "<OpenWQ> WARNING: PHREEQC BGC_GENERAL_MOBILE_SPECIES name '"
                + mobile_name
                + "' not found in PHREEQC component list. Skipping.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }

    }

}
