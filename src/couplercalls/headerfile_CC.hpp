
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


#ifndef OPENWQ_COUPLERCALLSH_INCLUDED
#define OPENWQ_COUPLERCALLSH_INCLUDED

#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_vars.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_hostModelConfig.hpp"

#include "readjson/headerfile_RJSON.hpp"
#include "initiate/headerfile_INIT.hpp"
#include "extwatflux_ss/headerfile_EWF_SS.hpp"
#include "units/headerfile_units.hpp"
#include "compute/headerfile_compute.hpp"
#include "output/headerfile_OUT.hpp"
#include "utils/headerfile_UTILS.hpp"

#include "models_CH/headerfile_CH.hpp"
#include "models_TD/headerfile_TD.hpp"
#include "models_LE/headerfile_LE.hpp"
#include "models_SI/headerfile_SI.hpp"
#include "models_TS/headerfile_TS.hpp"

class OpenWQ_couplercalls{

    public:

    // #######################
    // Calls all functions needed for configuration
    // #######################
    void InitialConfig(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_TD_model& OpenWQ_TD_model,      // transport modules
        OpenWQ_LE_model& OpenWQ_LE_model,            // LE_model model
        OpenWQ_CH_model& OpenWQ_CH_model,                    // biochemistry modules
        OpenWQ_SI_model& OpenWQ_SI_model,
        OpenWQ_TS_model& OpenWQ_TS_model,
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,  // sink and source modules)
        OpenWQ_output& OpenWQ_output);

    // #######################
    // Calls all functions required inside time loop
    // But BEFORE space loop is initiated
    // #######################
    void RunTimeLoopStart(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_TD_model& OpenWQ_TD_model,      // transport modules
        OpenWQ_LE_model& OpenWQ_LE_model,           // LE_model model
        OpenWQ_CH_model& OpenWQ_CH_model,                   // biochemistry modules
        OpenWQ_SI_model& OpenWQ_SI_model,
        OpenWQ_TS_model& OpenWQ_TS_model,
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_compute& OpenWQ_compute,
        OpenWQ_output& OpenWQ_output,
        time_t simtime);                            // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC


    // #######################
    // Calls inside space loop
    // Called for each grid cell
    // there are two variations:
    // 1) Fluxes within the model domain: RunSpaceStep
    //    in this case, if it is an OUT-flux (loss from system), then
    //    set recipient = -1 (in this case ix_s, iy_s, iz_s values will be ignored) 
    // 2) IN-fluxes: RunSpaceStep_IN
    // #######################

    // 1) Fluxes within the model domain or Out-Fluxes (RunSpaceStep)
    void RunSpaceStep(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                       // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                     // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
        OpenWQ_readjson& OpenWQ_readjson,               // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,               // initiate modules
        OpenWQ_TD_model& OpenWQ_TD_model,         // transport modules
        OpenWQ_TS_model& OpenWQ_TS_model,
        OpenWQ_LE_model& OpenWQ_LE_model,               // LE_model model
        OpenWQ_CH_model& OpenWQ_CH_model,                       // biochemistry modules
        OpenWQ_SI_model& OpenWQ_SI_model,
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,     // sink and source modules)
        OpenWQ_compute& OpenWQ_compute,
        OpenWQ_output& OpenWQ_output,
        time_t simtime, // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
        const int source, const int ix_s, const int iy_s, const int iz_s,
        const int recipient, const int ix_r, const int iy_r, const int iz_r,
        const double wflux_s2r, const double wmass_source);

    // 2) IN-fluxes (external water fluxes): RunSpaceStep_IN
    void RunSpaceStep_IN(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                       // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                     // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
        OpenWQ_readjson& OpenWQ_readjson,               // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,               // initiate modules
        OpenWQ_TD_model& OpenWQ_TD_model,
        OpenWQ_CH_model& OpenWQ_CH_model,
        OpenWQ_TS_model& OpenWQ_TS_model,
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,     // sink and source modules)
        OpenWQ_compute& OpenWQ_compute,
        OpenWQ_output& OpenWQ_output,
        time_t simtime, // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
        std::string source_EWF_name,                    // name defined in HydroExtFlux (in couplecalls)
        const int recipient, const int ix_r, const int iy_r, const int iz_r,
        const double wflux_s2r);

    // #######################
    // Calls all functions required inside time loop
    // But AFTER space loop has been finalized
    // #######################
    void RunTimeLoopEnd(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_TD_model& OpenWQ_TD_model,      // transport modules
        OpenWQ_LE_model& OpenWQ_LE_model,           // LE_model model
        OpenWQ_CH_model& OpenWQ_CH_model,                   // biochemistry modules
        OpenWQ_SI_model& OpenWQ_SI_model,
        OpenWQ_TS_model& OpenWQ_TS_model,
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_compute& OpenWQ_compute,
        OpenWQ_output& OpenWQ_output,
        time_t simtime);                            // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC


};

#endif