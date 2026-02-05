

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
void OpenWQ_readjson::SetConfigInfo_SIModule_freundlich(
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

    // TODO: Adapt to this module (Freundlich sorption isotherm)
    // Note: BGC_GENERAL_MOBILE_SPECIES now uses species names (strings)
    // instead of integer indices. See bgc_flex/json_set.cpp for the
    // name-based resolution pattern.

}