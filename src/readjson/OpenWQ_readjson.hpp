
// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
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

#ifndef OpenWQ_READJSONH_INCLUDED
#define OpenWQ_READJSONH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>

// Linux only (I think)
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include "global/exprtk.hpp"

#include <cstdio>

#include "jnlohmann/json.hpp"
using json = nlohmann::json;

//#include "utility.h"

#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_hostModelconfig.hpp"
#include "units/OpenWQ_units.hpp"
#include "output/OpenWQ_output.hpp"
#include "utils/OpenWQ_utils.hpp"

class OpenWQ_readjson{

    public:

        void read_all(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_utils&  OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);


    private:
        void read_JSON_2class(
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_output& OpenWQ_output,
            OpenWQ_utils& OpenWQ_utils,
            json& jsondata,                         // JSON structure to save to
            const bool substruc_flag,               // Save in subfield of JSON structure? only if multiple files (e.g., source and sinks)        
            const std::string JsonSubStruct_name,   // if true, name of subfield    
            const std::string jsonfile);            // Name of JSON file

        void ConvertJSONtext_2upperCase(
            OpenWQ_utils& OpenWQ_utils,
            json& jsondata);

        void change_JSON_key_to_upper_case(
            OpenWQ_utils& OpenWQ_utils,
            json& object, 
            const std::string& old_key,
            std::string& new_key);

        void change_JSON_value_to_upper_case(
            OpenWQ_utils& OpenWQ_utils,
            json& object,
            std::string new_jsonkey_layer_1);

        void SetConfigInfo_driver(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_output_logFile(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_meta(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_compute(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_chemModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_solver(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);


        void SetConfigInfo_TEModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_output_driver(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_output_options(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_output_what2print(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_output& OpenWQ_output);

        void check_mkdir_openWQ(
            std::string& dirname);

};

#endif

