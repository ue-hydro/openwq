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

#include "headerfile_EWF_SS.hpp"
#include <algorithm> // for std::upper_bound (OPTIMIZED #11)


/* #################################################
// Check Sources/Sinks and EWF and Apply
// only for JSON and ASCII input
################################################# */
void OpenWQ_extwatflux_ss::CheckApply_EWFandSS_jsonAscii(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    const unsigned int YYYY,                            // current model step: Year
    const unsigned int MM,                              // current model step: month
    const unsigned int DD,                              // current model step: day
    const unsigned int HH,                              // current model step: hour
    const unsigned int MIN,                             // current model step: min
    const unsigned int SEC,                             // current model step: sec
    std::string inputType,                              // flag for SS or EWF
    std::unique_ptr<arma::Mat<double>>& array_FORC){    // array FORC arma (SS or EWF)
    
    /* ########################################
    // Data update/clean-up at 1st timestep
    // Applicable to both SS and EWF
    ######################################## */
    if (OpenWQ_wqconfig.is_tstep1()){

        // Remove requested loads that are prior to the simulation start datetime
        RemoveLoadBeforeSimStart_jsonAscii(
            OpenWQ_wqconfig,
            OpenWQ_units,
            array_FORC,
            YYYY, MM, DD, HH, MIN, SEC);

        // Update time increments for rows with "all" elements
        UpdateAllElemTimeIncremts(
            array_FORC,
            OpenWQ_wqconfig,
            OpenWQ_utils,
            OpenWQ_units,
            YYYY, MM, DD, HH, MIN, SEC);
    }

    // Convert sim time to time_t once
    const time_t simTime = OpenWQ_units.convertTime_ints2time_t(
        OpenWQ_wqconfig, YYYY, MM, DD, HH, MIN, SEC);

    // Get number of rows in FORC
    const unsigned int num_rowdata = (*array_FORC).n_rows;

    // OPTIMIZED: use cached flag instead of string comparison
    const unsigned int num_chem = OpenWQ_wqconfig.is_native_bgc_flex
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    const unsigned int num_ewf = OpenWQ_hostModelconfig.get_num_HydroExtFlux();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();
    const bool is_ss_input = (inputType.compare("ss") == 0);

    // Reset all ewf_conc values to ZERO for new time step (parallelized)
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (unsigned int idx = 0; idx < num_ewf * num_chem; idx++){
            const unsigned int ewfi = idx / num_chem;
            const unsigned int chemi = idx % num_chem;
            (*OpenWQ_vars.ewf_conc)(ewfi)(chemi).zeros();
        }
    }

    /* ########################################
    // Loop over row data in sink-source file
    // Process rows sequentially due to potential time-dependent updates
    ######################################## */

    for (unsigned int ri = 0; ri < num_rowdata; ri++){

        // Skip if row has already been used (not for continuous loads)
        if ((*array_FORC)(ri,14) == -2) continue;

        // ########################################
        // Check if time in array_FORC row ri matches the current model time
        
        // Get requested JSON datetime
        int YYYY_json = (*array_FORC)(ri,3);
        int MM_json   = (*array_FORC)(ri,4);
        int DD_json   = (*array_FORC)(ri,5);
        int HH_json   = (*array_FORC)(ri,6);
        int MIN_json  = (*array_FORC)(ri,7);
        int SEC_json  = (*array_FORC)(ri,8);

        // Track "all" flags
        bool anyAll_flag = false;
        const bool YYYYall_flag = (YYYY_json == -1);
        const bool MMall_flag   = (MM_json == -1);
        const bool DDall_flag   = (DD_json == -1);
        const bool HHall_flag   = (HH_json == -1);
        const bool MINall_flag  = (MIN_json == -1);
        const bool SECall_flag  = (SEC_json == -1);

        anyAll_flag = (YYYYall_flag || MMall_flag || DDall_flag || 
                       HHall_flag || MINall_flag || SECall_flag);

        // Add appropriate increments to "all" flag elements
        if (YYYYall_flag) YYYY_json += (*array_FORC)(ri,14);
        if (MMall_flag)   MM_json   += (*array_FORC)(ri,15);
        if (DDall_flag)   DD_json   += (*array_FORC)(ri,16);
        if (HHall_flag)   HH_json   += (*array_FORC)(ri,17);
        if (MINall_flag)  MIN_json  += (*array_FORC)(ri,18);
        if (SECall_flag)  SEC_json  += (*array_FORC)(ri,19);

        // Convert jsonTime to time_t
        const time_t jsonTime = OpenWQ_units.convertTime_ints2time_t(
            OpenWQ_wqconfig,
            YYYY_json, MM_json, DD_json, HH_json, MIN_json, SEC_json);

        // Skip if not time to load yet
        if (simTime < jsonTime) continue;

        // ########################################
        // Apply source/sink or update EWF concentration
        // ########################################

        // Adjust value based on time increment for continuous/discrete loads
        double value_adjust = (*array_FORC)(ri,12);

        if (SECall_flag && is_ss_input){
            const double timeDiffSecs = difftime(simTime, jsonTime);
            value_adjust *= timeDiffSecs;
        }

        // Get SS type (source=0 or sink=1)
        const long sinksource_flag = (*array_FORC)(ri,2);

        const unsigned int index = (*array_FORC)(ri,1);  // compartment or EWF index
        const unsigned int chemi = (*array_FORC)(ri,0);  // chemical index
        const int ix = (*array_FORC)(ri,9);
        const int iy = (*array_FORC)(ri,10);
        const int iz = (*array_FORC)(ri,11);

        // Apply based on type
        if (is_ss_input && sinksource_flag == 0){
            // Source
            Apply_Source(OpenWQ_vars, OpenWQ_wqconfig, OpenWQ_hostModelconfig,
                        OpenWQ_output, index, chemi, ix, iy, iz, value_adjust);
        }
        else if (is_ss_input && sinksource_flag == 1){
            // Sink
            Apply_Sink(OpenWQ_vars, OpenWQ_wqconfig, OpenWQ_hostModelconfig,
                      OpenWQ_output, index, chemi, ix, iy, iz, value_adjust);
        }
        else if (!is_ss_input){
            // EWF
            Update_EWFconc_jsonAscii(OpenWQ_vars, OpenWQ_wqconfig, OpenWQ_hostModelconfig,
                                     OpenWQ_output, index, chemi, ix, iy, iz, value_adjust);
        }

        // ########################################
        // Prepare time increments for next load
        // ########################################

        // Mark as used if no "all" flags
        if (!anyAll_flag) {
            (*array_FORC)(ri,14) = -2;
            continue;
        }

        // Update increment for rows with "all" flags
        time_t updated_jsonTime = jsonTime;
        bool addAnyIncrem_flag = true;

        // Get max days in month
        const int DD_max = OpenWQ_utils.getNumberOfDaysInMonthYear(YYYY_json, MM_json);

        while (simTime >= updated_jsonTime){

            addAnyIncrem_flag = true;

            // Increment time fields in order of granularity
            if (SECall_flag && addAnyIncrem_flag && SEC_json < 59){
                SEC_json++;
                addAnyIncrem_flag = false;
            }
            else if (MINall_flag && addAnyIncrem_flag && MIN_json < 59){
                MIN_json++;
                if (SECall_flag && SEC_json == 59) SEC_json = 0;
                addAnyIncrem_flag = false;
            }
            else if (HHall_flag && addAnyIncrem_flag && HH_json < 23){
                HH_json++;
                if (SECall_flag && SEC_json == 59) SEC_json = 0;
                if (MINall_flag && MIN_json == 59) MIN_json = 0;
                addAnyIncrem_flag = false;
            }
            else if (DDall_flag && addAnyIncrem_flag && DD_json < DD_max){
                DD_json++;
                if (SECall_flag && SEC_json == 59) SEC_json = 0;
                if (MINall_flag && MIN_json == 59) MIN_json = 0;
                if (HHall_flag && HH_json == 23) HH_json = 0;
                addAnyIncrem_flag = false;
            }
            else if (MMall_flag && addAnyIncrem_flag && MM_json < 12){
                MM_json++;
                if (SECall_flag && SEC_json == 59) SEC_json = 0;
                if (MINall_flag && MIN_json == 59) MIN_json = 0;
                if (HHall_flag && HH_json == 23) HH_json = 0;
                if (DDall_flag && DD_json == DD_max) DD_json = 1;
                addAnyIncrem_flag = false;
            }
            else if (YYYYall_flag && addAnyIncrem_flag){
                YYYY_json++;
                if (SECall_flag && SEC_json == 59) SEC_json = 0;
                if (MINall_flag && MIN_json == 59) MIN_json = 0;
                if (HHall_flag && HH_json == 23) HH_json = 0;
                if (DDall_flag && DD_json == DD_max) DD_json = 1;
                if (MMall_flag && MM_json == 12) MM_json = 1;
                addAnyIncrem_flag = false;
            }

            // Exit if no increment possible
            if (addAnyIncrem_flag) break;

            // Compute new jsonTime
            updated_jsonTime = OpenWQ_units.convertTime_ints2time_t(
                OpenWQ_wqconfig,
                YYYY_json, MM_json, DD_json, HH_json, MIN_json, SEC_json);
        }

        // Update or mark as used
        if (addAnyIncrem_flag){
            (*array_FORC)(ri,14) = -2;
        }
        else {
            // Update increments (+1 because "all" fields use -1 as flag)
            if (YYYYall_flag) (*array_FORC)(ri,14) = YYYY_json + 1;
            if (MMall_flag)   (*array_FORC)(ri,15) = MM_json + 1;
            if (DDall_flag)   (*array_FORC)(ri,16) = DD_json + 1;
            if (HHall_flag)   (*array_FORC)(ri,17) = HH_json + 1;
            if (MINall_flag)  (*array_FORC)(ri,18) = MIN_json + 1;
            if (SECall_flag)  (*array_FORC)(ri,19) = SEC_json + 1;
        }
    }
}

/* #################################################
// Check EWF and apply
// only H5 input
################################################# */
void OpenWQ_extwatflux_ss::CheckApply_EWF_h5(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    const unsigned int YYYY,
    const unsigned int MM,
    const unsigned int DD,
    const unsigned int HH,
    const unsigned int MIN,
    const unsigned int SEC){

    // Convert sim time to time_t once
    const time_t simTime = OpenWQ_units.convertTime_ints2time_t(
        OpenWQ_wqconfig, YYYY, MM, DD, HH, MIN, SEC);

    const unsigned long long num_ewfh5_requests = 
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time).size();

    /* ########################################
    // Data update/clean-up at 1st timestep
    ######################################## */
    if (OpenWQ_wqconfig.is_tstep1()){
        for (unsigned int reqi = 0; reqi < num_ewfh5_requests; reqi++){
            RemoveLoadBeforeSimStart_h5(
                OpenWQ_wqconfig,
                OpenWQ_units,
                OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data,
                OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time,
                reqi,
                YYYY, MM, DD, HH, MIN, SEC);
        }
    }

    // Cache interpolation method check
    const bool is_step_interp = OpenWQ_wqconfig.is_h5EWF_interpMethod("STEP");
    const bool is_linear_interp = OpenWQ_wqconfig.is_h5EWF_interpMethod("LINEAR");
    const bool is_nearest_interp = OpenWQ_wqconfig.is_h5EWF_interpMethod("NEAREST");

    // Loop over all requests
    for (unsigned int reqi = 0; reqi < num_ewfh5_requests; reqi++){

        const unsigned long long num_chems = 
            (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data)[reqi].size();
        const unsigned long long num_timeStamps = 
            (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time)[reqi].size();

        /* ########################################
        // OPTIMIZED #11: Find timestamp at or above current simTime using binary search
        ######################################## */

        const auto& time_vec = (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time)[reqi];
        // std::upper_bound finds first element > simTime (timestamps are sorted)
        auto it = std::upper_bound(time_vec.begin(), time_vec.end(), simTime);

        if (it != time_vec.end()){
            const unsigned long long tStamp = std::distance(time_vec.begin(), it);
            const time_t h5EWF_time_after = *it;

            // Equivalent to original: found first timestamp > simTime
            {

                // Get timestamps for interpolation
                time_t h5EWF_time_before;
                unsigned long long tStamp_before, tStamp_after;

                if (tStamp != 0){
                    tStamp_before = tStamp - 1;
                    tStamp_after = tStamp;
                    h5EWF_time_before = 
                        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time)[reqi][tStamp_before];
                } else {
                    tStamp_before = tStamp;
                    tStamp_after = tStamp;
                    h5EWF_time_before = h5EWF_time_after;
                }

                // Pre-compute interpolation weights for LINEAR method
                double interp_weight = 0.0;
                if (is_linear_interp && h5EWF_time_after != h5EWF_time_before){
                    interp_weight = static_cast<double>(simTime - h5EWF_time_before) / 
                                   static_cast<double>(h5EWF_time_after - h5EWF_time_before);
                }

                // Loop over chemicals
                for (unsigned long long chemi = 0; chemi < num_chems; chemi++){

                    arma::cube h5Conc_chemi_interp;

                    // Get concentrations
                    const arma::cube& h5Conc_chemi_before = 
                        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data)[reqi][chemi][tStamp_before];

                    // Interpolation
                    if (is_step_interp){
                        h5Conc_chemi_interp = h5Conc_chemi_before;
                    }
                    else if (is_linear_interp){
                        const arma::cube& h5Conc_chemi_after = 
                            (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data)[reqi][chemi][tStamp_after];
                        h5Conc_chemi_interp = h5Conc_chemi_before + 
                            (h5Conc_chemi_after - h5Conc_chemi_before) * interp_weight;
                    }
                    else if (is_nearest_interp){
                        const double half_interval = 
                            static_cast<double>(h5EWF_time_after - h5EWF_time_before) / 2.0;
                        if (static_cast<double>(simTime - h5EWF_time_before) > half_interval){
                            h5Conc_chemi_interp = 
                                (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data)[reqi][chemi][tStamp_after];
                        } else {
                            h5Conc_chemi_interp = h5Conc_chemi_before;
                        }
                    }
                    else {
                        // Unknown method - default to STEP
                        OpenWQ_wqconfig.set_h5EWF_interpMethod("STEP");
                        h5Conc_chemi_interp = h5Conc_chemi_before;

                        std::string msg_string = OpenWQ_wqconfig.h5EWF_interp_warning_msg();
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    }

                    // Update concentrations
                    Update_EWFconc_h5(
                        OpenWQ_vars,
                        OpenWQ_wqconfig,
                        OpenWQ_hostModelconfig,
                        OpenWQ_output,
                        reqi,
                        chemi,
                        h5Conc_chemi_interp);
                }

            }
        }
    }
}

/* #################################################
 // Apply Source
 ################################################# */
void OpenWQ_extwatflux_ss::Apply_Source(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output,
    const unsigned int cmpi,
    const unsigned int chemi,
    int ix,
    int iy,
    int iz,
    const double ss_data_json){

    // Get compartment dimensions
    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(cmpi);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(cmpi);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(cmpi);

    // Determine domain region
    const unsigned int spX_min = (ix != -1) ? ix : 0;
    const unsigned int spX_max = (ix != -1) ? ix : nx - 1;
    const unsigned int spY_min = (iy != -1) ? iy : 0;
    const unsigned int spY_max = (iy != -1) ? iy : ny - 1;
    const unsigned int spZ_min = (iz != -1) ? iz : 0;
    const unsigned int spZ_max = (iz != -1) ? iz : nz - 1;

    try{
        (*OpenWQ_vars.d_chemass_ss)(cmpi)(chemi)(
            arma::span(spX_min, spX_max),
            arma::span(spY_min, spY_max),
            arma::span(spZ_min, spZ_max)) += ss_data_json;

    } catch (...) {

        std::string chemname = OpenWQ_wqconfig.is_native_bgc_flex
            ? (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi]
            : (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list)[chemi];

        std::string msg_string =
            "<OpenWQ> WARNING: Sink/Source load out of boundaries. "
            "Requested load ignored: "
            "Compartment=" + OpenWQ_hostModelconfig.get_HydroComp_name_at(cmpi) +
            ", Chemical=" + chemname +
            ", ix=" + std::to_string(ix) +
            ", iy=" + std::to_string(iy) +
            ", iz=" + std::to_string(iz);

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }
}

/* #################################################
 // Apply Sink
 ################################################# */
void OpenWQ_extwatflux_ss::Apply_Sink(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output,
    const unsigned int cmpi,
    const unsigned int chemi,
    int ix,
    int iy,
    int iz,
    const double ss_data_json){

    // Get compartment dimensions
    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(cmpi);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(cmpi);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(cmpi);

    // Determine domain region
    const unsigned int spX_min = (ix != -1) ? ix : 0;
    const unsigned int spX_max = (ix != -1) ? ix : nx - 1;
    const unsigned int spY_min = (iy != -1) ? iy : 0;
    const unsigned int spY_max = (iy != -1) ? iy : ny - 1;
    const unsigned int spZ_min = (iz != -1) ? iz : 0;
    const unsigned int spZ_max = (iz != -1) ? iz : nz - 1;

    try{
        (*OpenWQ_vars.d_chemass_ss)(cmpi)(chemi)(
            arma::span(spX_min, spX_max),
            arma::span(spY_min, spY_max),
            arma::span(spZ_min, spZ_max)) -= ss_data_json;

        // OPTIMIZED #14: Replace negative values with zero only in the modified sub-region
        // instead of transforming the entire cube
        (*OpenWQ_vars.d_chemass_ss)(cmpi)(chemi)(
            arma::span(spX_min, spX_max),
            arma::span(spY_min, spY_max),
            arma::span(spZ_min, spZ_max)).transform(
            [](double val) { return (val < 0.0) ? 0.0 : val; });

    } catch (...) {

        std::string chemname = OpenWQ_wqconfig.is_native_bgc_flex
            ? (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi]
            : (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list)[chemi];

        std::string msg_string =
            "<OpenWQ> WARNING: Sink/Source load out of boundaries. "
            "Requested load ignored: "
            "Compartment=" + OpenWQ_hostModelconfig.get_HydroComp_name_at(cmpi) +
            ", Chemical=" + chemname +
            ", ix=" + std::to_string(ix) +
            ", iy=" + std::to_string(iy) +
            ", iz=" + std::to_string(iz);

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }
}

/* #################################################
 // Update concentration if EWF: JSON and ASCII input
 ################################################# */
void OpenWQ_extwatflux_ss::Update_EWFconc_jsonAscii(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output,
    const unsigned int ewfi,
    const unsigned int chemi,
    int ix,
    int iy,
    int iz,
    const double new_concVal){

    // Get EWF dimensions
    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_x_at(ewfi);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_y_at(ewfi);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_z_at(ewfi);

    // Determine domain region
    const unsigned int spX_min = (ix != -1) ? ix : 0;
    const unsigned int spX_max = (ix != -1) ? ix : nx - 1;
    const unsigned int spY_min = (iy != -1) ? iy : 0;
    const unsigned int spY_max = (iy != -1) ? iy : ny - 1;
    const unsigned int spZ_min = (iz != -1) ? iz : 0;
    const unsigned int spZ_max = (iz != -1) ? iz : nz - 1;

    try{
        (*OpenWQ_vars.ewf_conc)(ewfi)(chemi)(
            arma::span(spX_min, spX_max),
            arma::span(spY_min, spY_max),
            arma::span(spZ_min, spZ_max)) += new_concVal;

        // Replace negative values with zero
        (*OpenWQ_vars.ewf_conc)(ewfi)(chemi).transform(
            [](double val) { return (val < 0.0) ? 0.0 : val; });

    } catch (...) {

        std::string chemname =
            ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0)
            ? (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi]
            : (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list)[chemi];

        std::string msg_string =
            "<OpenWQ> WARNING: EWF conc out of boundaries. "
            "Requested EWF concentration update using JSON/ASCII input ignored and set to zero: "
            "EWF=" + OpenWQ_hostModelconfig.get_HydroComp_name_at(ewfi) +
            ", Chemical=" + chemname +
            ", ix=" + std::to_string(ix) +
            ", iy=" + std::to_string(iy) +
            ", iz=" + std::to_string(iz);

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }
}

/* #################################################
 // Update concentration if EWF: H5
 ################################################# */
void OpenWQ_extwatflux_ss::Update_EWFconc_h5(
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output,
    const unsigned int reqi,
    const unsigned int chemi,
    arma::cube& h5Conc_chemi_interp){

    // Get external water flux id
    const unsigned int ewfi = (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_ewfCompID)[reqi];

    // Get interface dimensions
    const unsigned int nx_interf = 
        OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_x_at(ewfi) - 1;
    const unsigned int ny_interf = 
        OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_y_at(ewfi) - 1;
    const unsigned int nz_interf = 
        OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_z_at(ewfi) - 1;

    try{
        (*OpenWQ_vars.ewf_conc)(ewfi)(chemi)(
            arma::span(0, nx_interf),
            arma::span(0, ny_interf),
            arma::span(0, nz_interf)) = h5Conc_chemi_interp;

        // Replace negative values with zero
        (*OpenWQ_vars.ewf_conc)(ewfi)(chemi).transform(
            [](double val) { return (val < 0.0) ? 0.0 : val; });

    } catch (...) {

        std::string chemname =
            ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0)
            ? (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi]
            : (OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list)[chemi];

        std::string msg_string =
            "<OpenWQ> WARNING: EWF conc out of boundaries. "
            "Requested EWF concentration update using HDF5 input ignored and set to zero: "
            "EWF=" + OpenWQ_hostModelconfig.get_HydroComp_name_at(ewfi) +
            ", Chemical=" + chemname;

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }
}

/* #################################################
// Remove loads prior to simulation start (JSON/ASCII)
#################################################*/
void OpenWQ_extwatflux_ss::RemoveLoadBeforeSimStart_jsonAscii(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_units& OpenWQ_units,
    std::unique_ptr<arma::Mat<double>>& array_FORC,
    const int YYYY,
    const int MM,
    const int DD,
    const int HH,
    const int MIN,
    const int SEC){

    // Convert sim time to time_t once
    const time_t simTime = OpenWQ_units.convertTime_ints2time_t(
        OpenWQ_wqconfig, YYYY, MM, DD, HH, MIN, SEC);

    const unsigned int num_rowdata = (*array_FORC).n_rows;
    std::vector<int> rows2Remove;
    rows2Remove.reserve(num_rowdata / 10); // Reserve some space

    /* ########################################
    // Identify rows to remove
    ######################################## */

    for (unsigned int ri = 0; ri < num_rowdata; ri++){

        // Get requested JSON datetime
        int YYYY_json = (*array_FORC)(ri,3);
        int MM_json   = (*array_FORC)(ri,4);
        int DD_json   = (*array_FORC)(ri,5);
        int HH_json   = (*array_FORC)(ri,6);
        int MIN_json  = (*array_FORC)(ri,7);
        int SEC_json  = (*array_FORC)(ri,8);

        // Check for "all" flags
        const bool allinYYYY_flag = (YYYY_json == -1);
        bool all_flag = allinYYYY_flag;

        // Replace "all" (-1) with current sim time
        if (YYYY_json == -1){ YYYY_json = YYYY; all_flag = true; }
        if (MM_json == -1)  { MM_json = MM; all_flag = true; }
        if (DD_json == -1)  { DD_json = DD; all_flag = true; }
        if (HH_json == -1)  { HH_json = HH; all_flag = true; }
        if (MIN_json == -1) { MIN_json = MIN; all_flag = true; }
        if (SEC_json == -1) { SEC_json = SEC; all_flag = true; }

        // Convert jsonTime to time_t
        const time_t jsonTime = OpenWQ_units.convertTime_ints2time_t(
            OpenWQ_wqconfig,
            YYYY_json, MM_json, DD_json, HH_json, MIN_json, SEC_json);

        // Mark row for removal
        if ((!all_flag && jsonTime < simTime) ||
            (!allinYYYY_flag && jsonTime < simTime && YYYY_json < YYYY)) {
            rows2Remove.push_back(ri);
        }
    }

    // Convert to arma::uvec and remove rows
    if (!rows2Remove.empty()){
        const unsigned int n_rows2remove = rows2Remove.size();
        arma::uvec rows2Remove_uvec(n_rows2remove);

        for (unsigned int ri = 0; ri < n_rows2remove; ri++){
            rows2Remove_uvec(ri) = rows2Remove[ri];
        }

        (*array_FORC).shed_rows(rows2Remove_uvec);
    }
}

/* #################################################
// Remove loads prior to simulation start (H5)
#################################################*/
void OpenWQ_extwatflux_ss::RemoveLoadBeforeSimStart_h5(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_units& OpenWQ_units,
    std::unique_ptr<std::vector<std::vector<std::vector<arma::Cube<double>>>>>& FORC_vec_data,
    std::unique_ptr<std::vector<std::vector<time_t>>>& FORC_vec_time_t,
    const int reqi,
    const int YYYY,
    const int MM,
    const int DD,
    const int HH,
    const int MIN,
    const int SEC){

    // Convert sim time to time_t once
    const time_t simTime = OpenWQ_units.convertTime_ints2time_t(
        OpenWQ_wqconfig, YYYY, MM, DD, HH, MIN, SEC);

    const unsigned long long num_chems = (*FORC_vec_data)[reqi].size();
    const unsigned long long num_timeStamps = (*FORC_vec_time_t)[reqi].size();

    std::vector<int> timStampsIndex2Remove;
    timStampsIndex2Remove.reserve(num_timeStamps / 10);

    /* ########################################
    // Find timestamps before simTime
    ######################################## */

    for (unsigned long long tStamp = 0; tStamp < num_timeStamps; tStamp++){
        const time_t h5EWF_time = (*FORC_vec_time_t)[reqi][tStamp];
        if (h5EWF_time < simTime){
            timStampsIndex2Remove.push_back(tStamp);
        }
    }

    // Sort in descending order for safe removal
    std::sort(timStampsIndex2Remove.rbegin(), timStampsIndex2Remove.rend());

    /* ########################################
    // Remove identified timestamps
    ######################################## */

    for (const int ri2remove : timStampsIndex2Remove){
        (*FORC_vec_time_t)[reqi].erase(
            (*FORC_vec_time_t)[reqi].begin() + ri2remove);

        for (unsigned long long chemi = 0; chemi < num_chems; chemi++){
            (*FORC_vec_data)[reqi][chemi].erase(
                (*FORC_vec_data)[reqi][chemi].begin() + ri2remove);
        }
    }
}

/* #################################################
// Update time increments for rows with "all" elements
#################################################*/
void OpenWQ_extwatflux_ss::UpdateAllElemTimeIncremts(
    std::unique_ptr<arma::Mat<double>>& array_FORC,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    const int YYYY,
    const int MM,
    const int DD,
    const int HH,
    const int MIN,
    const int SEC){

    // Convert sim time to time_t once
    const time_t simTime = OpenWQ_units.convertTime_ints2time_t(
        OpenWQ_wqconfig, YYYY, MM, DD, HH, MIN, SEC);

    const unsigned int num_rowdata = (*array_FORC).n_rows;

    /* ########################################
    // Process each row
    ######################################## */

    for (unsigned int ri = 0; ri < num_rowdata; ri++){

        // Get requested JSON datetime
        int YYYY_json = (*array_FORC)(ri,3);
        int MM_json   = (*array_FORC)(ri,4);
        int DD_json   = (*array_FORC)(ri,5);
        int HH_json   = (*array_FORC)(ri,6);
        int MIN_json  = (*array_FORC)(ri,7);
        int SEC_json  = (*array_FORC)(ri,8);

        // Check for "all" flags
        const bool all_YYYY_flag = (YYYY_json == -1);
        const bool all_MM_flag   = (MM_json == -1);
        const bool all_DD_flag   = (DD_json == -1);
        const bool all_HH_flag   = (HH_json == -1);
        const bool all_MIN_flag  = (MIN_json == -1);
        const bool all_SEC_flag  = (SEC_json == -1);

        // Skip if no "all" elements
        if (!all_YYYY_flag && !all_MM_flag && !all_DD_flag &&
            !all_HH_flag && !all_MIN_flag && !all_SEC_flag){
            (*array_FORC)(ri,14) = -1;
            continue;
        }

        // Initialize increments
        unsigned int increm1 = 0, increm2 = 0, increm3 = 0;
        unsigned int increm4 = 0, increm5 = 0, increm6 = 0;

        // Compute initial increments based on simulation time
        if (!all_YYYY_flag && (YYYY_json > YYYY)){
            if (all_SEC_flag) increm6 = 1;
            if (all_MIN_flag) increm5 = 1;
            if (all_HH_flag)  increm4 = 1;
            if (all_DD_flag)  increm3 = 2;
            if (all_MM_flag)  increm2 = 2;
        }
        else if (!all_YYYY_flag && (YYYY_json == YYYY) &&
                 !all_MM_flag && (MM_json > MM)){
            if (all_SEC_flag) increm6 = 1;
            if (all_MIN_flag) increm5 = 1;
            if (all_HH_flag)  increm4 = 1;
            if (all_DD_flag)  increm3 = 2;
        }
        else if (!all_YYYY_flag && (YYYY_json == YYYY) &&
                 !all_MM_flag && (MM_json == MM) &&
                 !all_DD_flag && (DD_json > DD)){
            if (all_SEC_flag) increm6 = 1;
            if (all_MIN_flag) increm5 = 1;
            if (all_HH_flag)  increm4 = 1;
        }
        else if (!all_YYYY_flag && (YYYY_json == YYYY) &&
                 !all_MM_flag && (MM_json == MM) &&
                 !all_DD_flag && (DD_json == DD) &&
                 !all_HH_flag && (HH_json > HH)){
            if (all_SEC_flag) increm6 = 1;
            if (all_MIN_flag) increm5 = 1;
        }
        else if (!all_YYYY_flag && (YYYY_json == YYYY) &&
                 !all_MM_flag && (MM_json == MM) &&
                 !all_DD_flag && (DD_json == DD) &&
                 !all_HH_flag && (HH_json == HH) &&
                 !all_MIN_flag && (MIN_json > MIN)){
            if (all_SEC_flag) increm6 = 1;
        }
        else {
            if (all_SEC_flag)  increm6 = SEC - SEC_json;
            if (all_MIN_flag)  increm5 = MIN - MIN_json;
            if (all_HH_flag)   increm4 = HH - HH_json;
            if (all_DD_flag)   increm3 = DD - DD_json;
            if (all_MM_flag)   increm2 = MM - MM_json;
            if (all_YYYY_flag) increm1 = YYYY - YYYY_json;
        }

        // Compute initial jsonTime
        time_t jsonTime = OpenWQ_units.convertTime_ints2time_t(
            OpenWQ_wqconfig,
            YYYY_json + increm1,
            MM_json + increm2,
            DD_json + increm3,
            HH_json + increm4,
            MIN_json + increm5,
            SEC_json + increm6);

        // Refine increments to reach or exceed simTime
        if (jsonTime < simTime){

            // Try SEC increments
            if (all_SEC_flag){
                while (jsonTime < simTime && (SEC_json + increm6) < 59){
                    increm6++;
                    jsonTime = OpenWQ_units.convertTime_ints2time_t(
                        OpenWQ_wqconfig,
                        YYYY_json + increm1, MM_json + increm2, DD_json + increm3,
                        HH_json + increm4, MIN_json + increm5, SEC_json + increm6);
                }
                if (jsonTime < simTime) increm6 = 1;
            } else {
                increm6 = 1;
            }

            // Try MIN increments
            if (all_MIN_flag){
                while (jsonTime < simTime && (MIN_json + increm5) < 59){
                    increm5++;
                    jsonTime = OpenWQ_units.convertTime_ints2time_t(
                        OpenWQ_wqconfig,
                        YYYY_json + increm1, MM_json + increm2, DD_json + increm3,
                        HH_json + increm4, MIN_json + increm5, SEC_json + increm6);
                }
                if (jsonTime < simTime) increm5 = 1;
            } else {
                increm5 = 1;
            }

            // Try HH increments
            if (all_HH_flag){
                while (jsonTime < simTime && (HH_json + increm4) < 23){
                    increm4++;
                    jsonTime = OpenWQ_units.convertTime_ints2time_t(
                        OpenWQ_wqconfig,
                        YYYY_json + increm1, MM_json + increm2, DD_json + increm3,
                        HH_json + increm4, MIN_json + increm5, SEC_json + increm6);
                }
                if (jsonTime < simTime) increm4 = 1;
            } else {
                increm4 = 1;
            }

            // Try DD increments
            const int DD_max = OpenWQ_utils.getNumberOfDaysInMonthYear(YYYY_json, MM_json);
            if (all_DD_flag){
                while (jsonTime < simTime && (DD_json + increm3) < DD_max){
                    increm3++;
                    jsonTime = OpenWQ_units.convertTime_ints2time_t(
                        OpenWQ_wqconfig,
                        YYYY_json + increm1, MM_json + increm2, DD_json + increm3,
                        HH_json + increm4, MIN_json + increm5, SEC_json + increm6);
                }
                if (jsonTime < simTime) increm3 = 2;
            } else {
                increm3 = 1;
            }

            // Try MM increments
            if (all_MM_flag){
                while (jsonTime < simTime && (MM_json + increm2) < 12){
                    increm2++;
                    jsonTime = OpenWQ_units.convertTime_ints2time_t(
                        OpenWQ_wqconfig,
                        YYYY_json + increm1, MM_json + increm2, DD_json + increm3,
                        HH_json + increm4, MIN_json + increm5, SEC_json + increm6);
                }
                if (jsonTime < simTime) increm2 = 2;
            } else {
                increm2 = 1;
            }

            // Try YYYY increments
            if (all_YYYY_flag){
                while (jsonTime < simTime){
                    increm1++;
                    jsonTime = OpenWQ_units.convertTime_ints2time_t(
                        OpenWQ_wqconfig,
                        YYYY_json + increm1, MM_json + increm2, DD_json + increm3,
                        HH_json + increm4, MIN_json + increm5, SEC_json + increm6);
                }
            }
        }

        // Store computed increments
        (*array_FORC)(ri,14) = increm1;
        (*array_FORC)(ri,15) = increm2;
        (*array_FORC)(ri,16) = increm3;
        (*array_FORC)(ri,17) = increm4;
        (*array_FORC)(ri,18) = increm5;
        (*array_FORC)(ri,19) = increm6;
    }
}