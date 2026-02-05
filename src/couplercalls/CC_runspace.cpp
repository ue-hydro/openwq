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


#include "headerfile_CC.hpp"


// ################################################################
// Calls inside space loop
// Called for each grid cell
// Applicable to:
// 1) Fluxes between compartments and compartment-domain elements (ix,iy,iz)
// 2) Out-Fluxes (fluxes that exit the system - loss from control volume)
// ################################################################
void OpenWQ_couplercalls::RunSpaceStep(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json& OpenWQ_json,                       // create OpenWQ_json object
    OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
    OpenWQ_units& OpenWQ_units,                     // functions for unit conversion
    OpenWQ_utils& OpenWQ_utils,                       // utility methods/functions
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
    time_t simtime,                                 // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r, const double wmass_source){


    // Local variables
    std::string msg_string;

    // Store simulation time for use by TS models (also accessible via SUNDIALS path)
    OpenWQ_wqconfig.set_simtime(simtime);

    // Return if flux or source_volume is zero
    if (wflux_s2r <= OpenWQ_hostModelconfig.get_waterflux_minlim() || 
    wmass_source <= OpenWQ_hostModelconfig.get_watervol_minlim())
        return;

    // #################################################
    // MODULES

    // #################################################
    // TD_model module
    // Mass transport with water (mobile material only)
    // #################################################

    // OPTIMIZED: use cached bool instead of string comparison
    if (OpenWQ_wqconfig.is_TD_enabled) {

        OpenWQ_TD_model.TD_driver_run(
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_vars,
            OpenWQ_output,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);
    
    }

    // OPTIMIZED: use cached bool instead of string comparison
    if (OpenWQ_wqconfig.is_LE_enabled) {

        OpenWQ_LE_model.LE_driver_run(
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
            OpenWQ_vars,
            OpenWQ_output,                               // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);

    }

    // OPTIMIZED: use cached bool instead of string comparison
    if (OpenWQ_wqconfig.is_TS_enabled) {

        OpenWQ_compute.sediment_calls.push_back(std::make_tuple(source, ix_s, iy_s, iz_s, recipient, ix_r, iy_r, iz_r, wflux_s2r, wmass_source, "TS_type_LE"));

        OpenWQ_TS_model.TS_driver_run(
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
            OpenWQ_vars,
            OpenWQ_utils,
            OpenWQ_output,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source,
            "TS_type_LE");
    }

}

// ################################################################
// Calls inside space loop
// Called for each grid cell
// Applicable to:
// 1) IN-fluxes (external water fluxes)
// ################################################################
void OpenWQ_couplercalls::RunSpaceStep_IN(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json& OpenWQ_json,                       // create OpenWQ_json object
    OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
    OpenWQ_units& OpenWQ_units,                     // functions for unit conversion
    OpenWQ_utils& OpenWQ_utils,                       // utility methods/functions
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
    const std::string& source_EWF_name,             // OPTIMIZED: pass by const ref instead of value
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r){
    
    // Add external loading from precipitation
    OpenWQ_TD_model.EWF_flux_driver_run(
        OpenWQ_hostModelconfig,
        OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
        OpenWQ_vars,
        OpenWQ_utils,
        OpenWQ_output,                               // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
        source_EWF_name,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);

    // OPTIMIZED: use cached bool instead of string comparison
    if (OpenWQ_wqconfig.is_TS_enabled) {

        OpenWQ_compute.sediment_calls.push_back(std::make_tuple(0, 0, 0, 0, recipient, ix_r, iy_r, iz_r, wflux_s2r, 0, "TS_type_EWF"));

        // Run TS_model model
        OpenWQ_TS_model.TS_driver_run(
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
            OpenWQ_vars,
            OpenWQ_utils,
            OpenWQ_output,
            0, 0, 0, 0,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, 0,
            "TS_type_EWF");
    }

}
