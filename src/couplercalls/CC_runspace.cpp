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
    OpenWQ_LE_model& OpenWQ_LE_model,               // LE model
    OpenWQ_CH_model& OpenWQ_CH_model,                       // biochemistry modules
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,     // sink and source modules)
    OpenWQ_compute& OpenWQ_compute,
    OpenWQ_output& OpenWQ_output,
    time_t simtime,                                 // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r, const double wmass_source){

    
    // Local variables
    std::string msg_string;

    // Return if flux or source_volume is zero
    if (wflux_s2r <= OpenWQ_hostModelconfig.get_waterflux_minlim() || 
    wmass_source <= OpenWQ_hostModelconfig.get_watervol_minlim())
        return;

    // #################################################
    // MODULES

    // #################################################
    // TD module
    // Mass transport with water (mobile material only)
    // #################################################

    // if TD_module = OPENWQ_NATIVE_TD_ADVDISP
    if ((OpenWQ_wqconfig.TD_module).compare("OPENWQ_NATIVE_TD_ADVDISP") == 0)
    {
        OpenWQ_TD_model.AdvDisp(
            OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);

    // if TD_module = NATIVE_TD_ADV
    }else if ((OpenWQ_wqconfig.TD_module).compare("NATIVE_TD_ADV") == 0)
    {
        OpenWQ_TD_model.Adv(
            OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);

    // if TD_module != NONE (and any of the others)
    }else if ((OpenWQ_wqconfig.TD_module).compare("NONE")!= 0)
    {

        // Create Message
        msg_string = 
            "<OpenWQ> ERROR: No TD_module found or unkown";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file
        
        // Abort (Fatal error)
        exit(EXIT_FAILURE);

    }

    // Run LE model
    OpenWQ_LE_model.LE_driver_run(
        OpenWQ_hostModelconfig,
        OpenWQ_json,                       // create OpenWQ_json object
        OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
        OpenWQ_units,                     // functions for unit conversion
        OpenWQ_utils,                       // utility methods/functions
        OpenWQ_readjson,               // read json files
        OpenWQ_vars,
        OpenWQ_initiate,               // initiate modules
        OpenWQ_TD_model,         // transport modules
        OpenWQ_LE_model,               // LE model
        OpenWQ_CH_model,                       // biochemistry modules
        OpenWQ_extwatflux_ss,     // sink and source modules)
        OpenWQ_compute,
        OpenWQ_output,
        simtime,                                 // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source);

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
    OpenWQ_TD_model& OpenWQ_TD_model,         // transport modules
    OpenWQ_LE_model& OpenWQ_LE_model,               // LE model
    OpenWQ_CH_model& OpenWQ_CH_model,                       // biochemistry modules
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,     // sink and source modules)
    OpenWQ_compute& OpenWQ_compute,
    OpenWQ_output& OpenWQ_output,
    time_t simtime, // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
    std::string source_EWF_name,                    // name defined in HydroExtFlux (in couplecalls)
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r){

    // Dummy variable
    double ewf_conc_chemi;  // concentraton of chemical at EWF at current time step
    unsigned int iewf;      // index of EWF
    unsigned int ichem_mob; // iterative dummy var for mobile chemical species

    // Find specific EWF index in EWF names
    iewf = OpenWQ_utils.FindStrIndexInVectStr(
        OpenWQ_hostModelconfig.get_HydroExtFlux_names(),
        source_EWF_name);

    unsigned int numspec = OpenWQ_wqconfig.BGC_general_mobile_species.size();

    // Loop over mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        // mobile chemical species index
        ichem_mob = OpenWQ_wqconfig.BGC_general_mobile_species[chemi];

        // Get appropriate chemi concentraton for this ewf and domain ix, iy, iz
        ewf_conc_chemi = (*OpenWQ_vars.ewf_conc)
            (iewf)
            (ichem_mob)
            (ix_r,iy_r,iz_r);

        // Advection and dispersion
        OpenWQ_TD_model.Adv_IN(
            OpenWQ_vars, OpenWQ_wqconfig,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r,                      // external water flux
            ichem_mob,                          // chemical id (from BGQ file and used in state-varibble data structures)
            ewf_conc_chemi);                // concentration taken from EWF json file

    }

}
