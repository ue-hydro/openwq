
// Copyright 2020, Diogo Costa, diogo.costa@uevora.pt
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

#include "global/openwq_json.hpp"
using json = nlohmann::json;

//#include "utility.h"

#include "global/openwq_json.hpp"
#include "global/openwq_wqconfig.hpp"
#include "global/openwq_hostmodelconfig.hpp"
#include "units/headerfile_UNITS.hpp"
#include "output/headerfile_OUT.hpp"
#include "utils/headerfile_UTILS.hpp"

class OpenWQ_readjson{

    public:

        void JSON_driver(
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

        void SetConfigInfo_OUT_logFile(
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

        void SetConfigInfo_solver(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_INPUT_general(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_INPUT_EF(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_INPUT_SS(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_OUT_driver(
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

        // ######################
        // Model options
        // ######################

        // ######################
        // BGC model

        // driver
        void SetConfigInfo_CHModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // options
        void SetConfigInfo_CHModule_BGC_FLEX(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // ######################
        // TD_model model

        // driver
        void SetConfigInfo_TDModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // options
        // ????  TODO: need to separate

        // ######################
        // LE_model model

        // driver
        void SetConfigInfo_LEModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // options
        void SetConfigInfo_LEModule_BOUNDMIX(
            OpenWQ_json &OpenWQ_json,
            OpenWQ_hostModelconfig & OpenWQ_hostModelconfig,
            OpenWQ_wqconfig &OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // ######################
        // ST model

        // driver
        void SetConfigInfo_TSModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // options
        void SetConfigInfo_TSModule_HBVsed_hype(
            OpenWQ_json &OpenWQ_json,
            OpenWQ_wqconfig &OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);
             
        void SetConfigInfo_TSModule_MMF_hype(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_json &OpenWQ_json,
            OpenWQ_wqconfig &OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

            

        // SI_model model

        // driver
        void SetConfigInfo_SIModule(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // options
        void SetConfigInfo_SIModule_freundlich(
            OpenWQ_json &OpenWQ_json,
            OpenWQ_wqconfig &OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        void SetConfigInfo_SIModule_langmuir(
            OpenWQ_json &OpenWQ_json,
            OpenWQ_wqconfig &OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

};

#endif

