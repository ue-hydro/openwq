
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

// Set output options
void OpenWQ_readjson::SetConfigInfo_OUT_driver(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_hostModelconfig & OpenWQ_hostModelconfig,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units &OpenWQ_units,
    OpenWQ_output& OpenWQ_output){

    // Set options and unit conversions
    SetConfigInfo_output_options(
        OpenWQ_json, OpenWQ_wqconfig, OpenWQ_utils, OpenWQ_units,
        OpenWQ_output);

    // Get compartments and chemicals to print
    SetConfigInfo_output_what2print(OpenWQ_json,
        OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_utils,
        OpenWQ_units,  OpenWQ_output);
    
}