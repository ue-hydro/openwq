

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
    json BGCjson_subStruct;
    json BGCjson_ChemList;
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

    // Get number of chemical species from BGC_json
    OpenWQ_wqconfig.CH_model->PHREEQC->num_chem = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->FindComponents();

    component_names = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponents();

    // Get mobile species 
    // reset index to start on zero
    errorMsgIdentifier = "BGQ file in CHEMICAL_SPECIES";
    BGCjson_mobileSpecies = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        BGCjson_subStruct, "BGC_GENERAL_MOBILE_SPECIES",
        errorMsgIdentifier,
        false);  // no abort

    for (unsigned int chemi = 0; chemi < BGCjson_mobileSpecies.size(); chemi++){
    
        OpenWQ_wqconfig.CH_model->PHREEQC->mobile_species.push_back(
            (int)BGCjson_mobileSpecies.at(chemi) - 1);

    }

    // Get chemical species list from BGC_json
    for (unsigned int chemi = 0; chemi < (OpenWQ_wqconfig.CH_model->PHREEQC->num_chem); chemi++)
    {
        (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list).push_back(
            component_names[chemi]);
    }

}
