// Copyright 2026, Diogo Costa
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

#pragma once

#include "global/OpenWQ_json.hpp"
#include "readjson/headerfile_nlohmann.hpp"

using json = nlohmann::json;

/* #################################################
 Project JSON files
################################################# */
class OpenWQ_json
{
    public:

    json Master;
    json Config;
    json ExtWatFlux;
    json SinkSource;
    json TD_module;     // Transport Dissolved module
    json TS_module;     // Transport Sediments module
    json LE_module;     // Lateral Exchange module
    json BGC_module;    // Biogeochemistry module
    json SI_module;     // Sorption Isotherm module

};
