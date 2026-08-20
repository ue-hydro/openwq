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

#include "output/headerfile_OUT.hpp"

// ############################################################################
// [REDESIGN] Time-series output layout.
//
// Previously openWQ created ONE new dataset per timestep (named by timestamp),
// so a multi-decade hourly run produced hundreds of thousands of datasets in a
// single group -> a huge link B-tree that corrupts on weakly-consistent
// filesystems (Docker bind mount, Lustre/NFS) and is slow.
//
// New layout (created ONCE per file, then appended to each step):
//   /concentrations : 2D double, shape [time x cells], chunked, time-unlimited
//   /timestamps     : 1D variable-length strings, shape [time], time-unlimited
// plus the existing one-shot metadata datasets (xyz_elements, xyz_elements_size,
// and the host-id label dataset). Only ~5 datasets per file regardless of run
// length, so the group metadata stays tiny and robust on any filesystem.
static void create_timeseries_datasets(hid_t fh, unsigned int num_cells,
                                       OpenWQ_wqconfig& wqcfg,
                                       const std::string& filename){
    if (num_cells == 0) return;

    // ---- /concentrations : [time x cells], unlimited in time, chunked ----
    hsize_t dims[2]    = {0, (hsize_t)num_cells};
    hsize_t maxdims[2] = {H5S_UNLIMITED, (hsize_t)num_cells};
    hid_t   space = H5Screate_simple(2, dims, maxdims);
    hid_t   dcpl  = H5Pcreate(H5P_DATASET_CREATE);
    // Target ~1 MB per chunk so per-step row appends accumulate in the chunk
    // cache and flush a chunk at a time (scales with domain size).
    hsize_t chunk_time = (hsize_t)(1048576UL / (num_cells * sizeof(double)));
    if (chunk_time < 1)    chunk_time = 1;
    if (chunk_time > 8192) chunk_time = 8192;
    hsize_t chunk[2] = {chunk_time, (hsize_t)num_cells};
    H5Pset_chunk(dcpl, 2, chunk);
    hid_t ds = H5Dcreate(fh, "concentrations", H5T_NATIVE_DOUBLE, space,
                         H5P_DEFAULT, dcpl, H5P_DEFAULT);
    wqcfg.conc_dsets[filename] = ds;   // keep open (closed in OpenWQ_wqconfig dtor)
    H5Pclose(dcpl);
    H5Sclose(space);

    // ---- /timestamps : 1D variable-length strings, unlimited, chunked ----
    hsize_t tdims[1] = {0};
    hsize_t tmax[1]  = {H5S_UNLIMITED};
    hid_t   tspace = H5Screate_simple(1, tdims, tmax);
    hid_t   tdcpl  = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t tchunk[1] = {1024};
    H5Pset_chunk(tdcpl, 1, tchunk);
    hid_t   stype = H5Tcopy(H5T_C_S1);
    H5Tset_size(stype, H5T_VARIABLE);
    hid_t   tds = H5Dcreate(fh, "timestamps", stype, tspace,
                            H5P_DEFAULT, tdcpl, H5P_DEFAULT);
    wqcfg.time_dsets[filename] = tds;  // keep open (closed in OpenWQ_wqconfig dtor)
    H5Tclose(stype);
    H5Pclose(tdcpl);
    H5Sclose(tspace);
}

int OpenWQ_output::writeResults(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_compute& OpenWQ_compute,
    time_t simtime){

    // Create time string
    struct tm *tm_simtime = gmtime(&simtime);
    char timechar[30];
    strftime(timechar, 30, "%Y%b%d-%H:%M:%S", tm_simtime);
    std::string timestr(timechar);

    // Convert to upper case
    try {
        std::transform(timestr.begin(), timestr.end(), timestr.begin(),
                      [](unsigned char c) { return std::toupper(c); });
    } catch (...) {}

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
    const bool is_csv_output = OpenWQ_wqconfig.is_output_type_csv();
    const bool is_hdf5_output = OpenWQ_wqconfig.is_output_type_hdf5();
    const bool debug_mode = OpenWQ_wqconfig.debug_mode;
    const bool print_once = OpenWQ_wqconfig.print_oneStep;

    std::string outputfile_type = is_csv_output ? "CSV" : "HDF5";

    /* ########################################
    // Loop over compartments
    ######################################## */

    for (unsigned int icmp = 0; icmp < num_comps; icmp++){

        // Check if this compartment should be printed
        auto it = std::find(OpenWQ_wqconfig.compt2print.begin(),
                           OpenWQ_wqconfig.compt2print.end(), icmp);

        if (it == OpenWQ_wqconfig.compt2print.end()) continue;

        // ################################################
        // CSV Output
        if (is_csv_output){
            
            // Main output
            std::string output_label = "main";
            writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                    OpenWQ_vars.chemass, output_label, timestr, icmp);

            // Debug outputs
            if (debug_mode){
                
                output_label = "d_output_dt_chemistry";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_dt_chem_out, output_label, timestr, icmp);

                output_label = "d_output_dt_sorption";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_dt_sorpt_out, output_label, timestr, icmp);

                output_label = "d_output_dt_transport_diss";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_dt_transp_diss_out, output_label, timestr, icmp);

                output_label = "d_output_dt_transport_part";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_dt_transp_part_out, output_label, timestr, icmp);

                output_label = "d_output_ss";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_ss_out, output_label, timestr, icmp);

                output_label = "d_output_ewf";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_ewf_out, output_label, timestr, icmp);

                if (print_once){
                    output_label = "d_output_ic";
                    writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                            OpenWQ_vars.d_chemass_ic, output_label, timestr, icmp);
                }
            }
        }
        // ################################################
        // HDF5 Output
        else if (is_hdf5_output){

            // Main output
            std::string output_label = "main";
            writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                     OpenWQ_vars.chemass, output_label, timestr, icmp);

            // Sediment output if applicable
            if (export_sediment &&
                (OpenWQ_wqconfig.TS_model->ErodTranspCmpt).compare(
                    OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)) == 0){
                writeHDF5_Sediment(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                                  OpenWQ_vars.sedmass, output_label, timestr, icmp);

                // Sediment derivative channels (debug mode)
                if (debug_mode){
                    output_label = "d_output_sed_transport";
                    writeHDF5_Sediment(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                                      OpenWQ_vars.d_sedmass_transport_dt_out, output_label, timestr, icmp);
                    output_label = "main";
                }
            }
            // Erosion source channel lives on the SEDIMENT compartment (may
            // differ from the transport compartment)
            if (export_sediment && debug_mode &&
                (OpenWQ_wqconfig.TS_model->SedCmpt).compare(
                    OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)) == 0){
                output_label = "d_output_sed_mobilized";
                writeHDF5_Sediment(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                                  OpenWQ_vars.d_sedmass_mobilized_dt_out, output_label, timestr, icmp);
                output_label = "main";
            }

            // Debug outputs
            if (debug_mode){

                output_label = "d_output_dt_chemistry";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_dt_chem_out, output_label, timestr, icmp);

                output_label = "d_output_dt_sorption";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_dt_sorpt_out, output_label, timestr, icmp);

                output_label = "d_output_dt_transport_diss";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_dt_transp_diss_out, output_label, timestr, icmp);

                output_label = "d_output_dt_transport_part";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_dt_transp_part_out, output_label, timestr, icmp);

                output_label = "d_output_ss";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_ss_out, output_label, timestr, icmp);

                output_label = "d_output_ewf";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_ewf_out, output_label, timestr, icmp);

                if (print_once){
                    output_label = "d_output_ic";
                    writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                             OpenWQ_vars.d_chemass_ic, output_label, timestr, icmp);
                }
            }
        }
    }

    // Turn off one-step printing
    OpenWQ_wqconfig.print_oneStep = false;

    // Update next output time
    OpenWQ_wqconfig.update_nexttime_out();

    // Reset cumulative derivatives
    OpenWQ_compute.Reset_Deriv(
        OpenWQ_hostModelconfig,
        OpenWQ_wqconfig,
        OpenWQ_vars,
        false,  // reset inst derivatives
        true);  // reset cumulative derivatives

    // Log success message
    std::string msg_string = "<OpenWQ> Output export successful (" +
                            outputfile_type + "): " + timestr;
    ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    return EXIT_SUCCESS;
}

/* ########################################
// CSV output
######################################## */
int OpenWQ_output::writeCSV(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    std::unique_ptr<arma::field<arma::field<arma::cube>>>& OpenWQ_var2print,
    std::string& output_file_label,
    std::string timestr,
    int icmp){

    // Get number of cells to print
    const unsigned int num_cells2print = OpenWQ_wqconfig.cells2print_vec.at(icmp).n_rows;
    const bool printflag = OpenWQ_wqconfig.cells2print_bool.at(icmp);

    // Early exit if nothing to print
    if (!printflag || num_cells2print == 0) return EXIT_SUCCESS;

    // Cache frequently accessed values
    const std::string CompName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
    const unsigned int num_chem2print = OpenWQ_wqconfig.chem2print.size();
    const bool is_conc_requested = OpenWQ_wqconfig.is_concentration_requested();
    const double watervol_minlim = OpenWQ_hostModelconfig.get_watervol_minlim();
    const double noWaterConc = OpenWQ_wqconfig.noWaterConc;
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Get unit multipliers
    const double unit_mult_num = OpenWQ_wqconfig.get_output_units_numerator();
    const double unit_mult_den = OpenWQ_wqconfig.get_output_units_denominator();
    const std::string output_units = OpenWQ_wqconfig.get_output_units();

    // Build filename
    std::string filename = OpenWQ_wqconfig.get_output_dir() + "/" +
                          CompName_icmp + "@" + timestr + "-" +
                          output_file_label + ".txt";

    // Initialize output matrix
    arma::dmat filedata(num_cells2print, num_chem2print + 3);

    // Define header
    arma::field<std::string> header(num_chem2print + 3);
    header(0) = "ix [-]";
    header(1) = "iy [-]";
    header(2) = "iz [-]";

    // Get chemical names for header
    for (unsigned int ichem = 0; ichem < num_chem2print; ichem++){
        const std::string chem_name = (*OpenWQ_wqconfig.cached_chem_species_list_ptr)[OpenWQ_wqconfig.chem2print[ichem]];
        
        header(ichem + 3) = chem_name + "#" + output_units;
    }

    // Parallelize data extraction
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (unsigned int ixyz = 0; ixyz < num_cells2print; ixyz++){

            // Get indices (convert to 1-based for output)
            const unsigned int ix = OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 0);
            const unsigned int iy = OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 1);
            const unsigned int iz = OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 2);

            filedata(ixyz, 0) = ix + 1;
            filedata(ixyz, 1) = iy + 1;
            filedata(ixyz, 2) = iz + 1;

            // Get water volume once per cell
            const double water_vol = is_conc_requested ?
                OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ix, iy, iz) : 1.0;

            // Process all chemicals for this cell
            for (unsigned int ichem = 0; ichem < num_chem2print; ichem++){

                if (water_vol > watervol_minlim){
                    filedata(ixyz, ichem + 3) =
                        (*OpenWQ_var2print)(icmp)(OpenWQ_wqconfig.chem2print[ichem])(ix, iy, iz) *
                        unit_mult_num / (water_vol * unit_mult_den);
                } else {
                    filedata(ixyz, ichem + 3) = noWaterConc;
                }
            }
        }
    }

    // Write to file
    filedata.save(arma::csv_name(filename, header));

    return EXIT_SUCCESS;
}

/* ########################################
// HDF5 format
######################################## */
int OpenWQ_output::writeHDF5(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    std::unique_ptr<arma::field<arma::field<arma::cube>>>& OpenWQ_var2print,
    std::string& output_file_label,
    std::string timestr,
    int icmp){

    // Get number of cells to print
    const unsigned int num_cells2print = OpenWQ_wqconfig.cells2print_vec.at(icmp).n_rows;
    const bool printflag = OpenWQ_wqconfig.cells2print_bool.at(icmp);

    // Early exit if nothing to print
    if (!printflag || num_cells2print == 0) return EXIT_SUCCESS;

    // Cache frequently accessed values
    const std::string CompName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
    const unsigned int num_chem2print = OpenWQ_wqconfig.chem2print.size();
    const bool is_conc_requested = OpenWQ_wqconfig.is_concentration_requested();
    const double watervol_minlim = OpenWQ_hostModelconfig.get_watervol_minlim();
    const double noWaterConc = OpenWQ_wqconfig.noWaterConc;
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Get unit multipliers
    const double unit_mult_num = OpenWQ_wqconfig.get_output_units_numerator();
    const double unit_mult_den = OpenWQ_wqconfig.get_output_units_denominator();

    // Units string (replace "/" with "|")
    std::string units_string = OpenWQ_wqconfig.get_output_units();
    std::size_t slash_pos = units_string.find("/");
    if (slash_pos != std::string::npos){
        units_string.replace(slash_pos, 1, "|");
    }

    // Prepare compartment size info
    arma::mat cells2print_xyzElements_size(1, 3);
    cells2print_xyzElements_size(0, 0) = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
    cells2print_xyzElements_size(0, 1) = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
    cells2print_xyzElements_size(0, 2) = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

    /* ########################################
    // Save metadata & open files on first timestep.
    // Files stay open in OpenWQ_wqconfig.files for the
    // entire simulation (closed in destructor).
    // No close-reopen cycle = no Docker volume sync issues.
    ######################################## */

    if (OpenWQ_wqconfig.print_oneStep){

        // Prepare xyz elements (convert to 1-based indexing)
        arma::mat cells2print_xyzElements = OpenWQ_wqconfig.cells2print_vec[icmp];
        cells2print_xyzElements.for_each([](arma::mat::elem_type& val) { val += 1.0; });

        // Prepare host IDs
        std::vector<std::string> host_ids;
        host_ids.reserve(num_cells2print);
        for (unsigned int ixyz = 0; ixyz < num_cells2print; ixyz++){
            const int ix = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 0));
            const int iy = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 1));
            const int iz = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 2));
            host_ids.push_back(OpenWQ_hostModelconfig.get_cellid_to_wq_at(icmp, ix, iy, iz));
        }
        const std::string internal_db_name = OpenWQ_hostModelconfig.get_cellid_to_wqlabel();
        std::vector<const char*> cstrs(host_ids.size());
        for (size_t i = 0; i < host_ids.size(); ++i) cstrs[i] = host_ids[i].c_str();

        // File properties: latest format (fractal heaps), no locking
        hid_t fcpl = H5Pcreate(H5P_FILE_CREATE);
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        // Locking off: required on bind-mount / network filesystems that don't
        // support file locking. Otherwise use HDF5's DEFAULT format and DEFAULT
        // metadata cache — the new few-dataset layout (one extensible dataset per
        // file, created below) keeps metadata tiny, so the LATEST/fractal-heap
        // format and the no-evict metadata cache are no longer needed. (The
        // LATEST-format extensible-array chunk index combined with frequent
        // H5Dset_extent was an unnecessary risk.)
        H5Pset_file_locking(fapl, 0, 0);

        for (unsigned int ichem = 0; ichem < num_chem2print; ichem++){

            const std::string chem_name = (*OpenWQ_wqconfig.cached_chem_species_list_ptr)[OpenWQ_wqconfig.chem2print[ichem]];
            std::string filename = OpenWQ_wqconfig.get_output_dir() + "/" +
                                  CompName_icmp + "@" + chem_name + "#" +
                                  units_string + "-" + output_file_label + ".h5";

            hid_t fh = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, fcpl, fapl);
            if (fh < 0) {
                std::string msg = "<OpenWQ> ERROR: Failed to create HDF5 file: " + filename;
                ConsoleLog(OpenWQ_wqconfig, msg, true, true);
                continue;
            }

            // Write xyz_elements
            {
                hsize_t dims[2] = {cells2print_xyzElements.n_cols,
                                   cells2print_xyzElements.n_rows};
                hid_t sp = H5Screate_simple(2, dims, NULL);
                hid_t ds = H5Dcreate(fh, "xyz_elements", H5T_NATIVE_DOUBLE,
                                     sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (ds >= 0) {
                    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, cells2print_xyzElements.memptr());
                    H5Dclose(ds);
                }
                H5Sclose(sp);
            }

            // Write xyz_elements_size
            {
                hsize_t dims[2] = {cells2print_xyzElements_size.n_cols,
                                   cells2print_xyzElements_size.n_rows};
                hid_t sp = H5Screate_simple(2, dims, NULL);
                hid_t ds = H5Dcreate(fh, "xyz_elements_size", H5T_NATIVE_DOUBLE,
                                     sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (ds >= 0) {
                    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, cells2print_xyzElements_size.memptr());
                    H5Dclose(ds);
                }
                H5Sclose(sp);
            }

            // Write host IDs as variable-length strings
            if (!host_ids.empty()){
                hsize_t dims[1] = {host_ids.size()};
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dtype = H5Tcopy(H5T_C_S1);
                H5Tset_size(dtype, H5T_VARIABLE);
                hid_t dset = H5Dcreate(fh, internal_db_name.c_str(), dtype, space,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset >= 0){
                    H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cstrs.data());
                    H5Dclose(dset);
                }
                H5Tclose(dtype);
                H5Sclose(space);
            }

            // [REDESIGN] create the time-extensible datasets that each timestep
            // appends to (/concentrations + /timestamps)
            create_timeseries_datasets(fh, num_cells2print, OpenWQ_wqconfig, filename);

            // Store handle — file stays open until destructor
            OpenWQ_wqconfig.files[filename] = fh;

            std::string msg_string = "<OpenWQ> Created output file " + filename;
            ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }

        H5Pclose(fcpl);
        H5Pclose(fapl);
    }

    /* ########################################
    // Write data for current timestep (parallel).
    // Files are already open in OpenWQ_wqconfig.files.
    ######################################## */

    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (unsigned int ichem = 0; ichem < num_chem2print; ichem++){

            const std::string chem_name = (*OpenWQ_wqconfig.cached_chem_species_list_ptr)[OpenWQ_wqconfig.chem2print[ichem]];
            std::string filename = OpenWQ_wqconfig.get_output_dir() + "/" +
                                  CompName_icmp + "@" + chem_name + "#" +
                                  units_string + "-" + output_file_label + ".h5";

            // Look up persistent handle (read-only map access, thread-safe)
            auto it = OpenWQ_wqconfig.files.find(filename);
            if (it == OpenWQ_wqconfig.files.end() || it->second < 0) {
                #pragma omp critical
                {
                    std::string msg = "<OpenWQ> WARNING: HDF5 handle not found, skipping: " + filename;
                    ConsoleLog(OpenWQ_wqconfig, msg, true, true);
                }
                continue;
            }

            // Prepare data matrix
            arma::mat data2print(num_cells2print, 1);

            // Extract data
            for (unsigned int celli = 0; celli < num_cells2print; celli++){

                const unsigned int ix = OpenWQ_wqconfig.cells2print_vec[icmp](celli, 0);
                const unsigned int iy = OpenWQ_wqconfig.cells2print_vec[icmp](celli, 1);
                const unsigned int iz = OpenWQ_wqconfig.cells2print_vec[icmp](celli, 2);

                const double water_vol = is_conc_requested ?
                    OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ix, iy, iz) : 1.0;

                if ((is_conc_requested && water_vol > watervol_minlim) || !is_conc_requested){
                    data2print(celli, 0) =
                        (*OpenWQ_var2print)(icmp)(OpenWQ_wqconfig.chem2print[ichem])(ix, iy, iz) *
                        unit_mult_num / (water_vol * unit_mult_den);
                } else {
                    data2print(celli, 0) = noWaterConc;
                }
            }

            // [FIX] libhdf5_serial (1.10.10) is NOT thread-safe. This loop runs
            // under "#pragma omp parallel for", so the actual HDF5 write must be
            // serialized — otherwise concurrent threads corrupt HDF5's global
            // metadata cache (HDF5 error: "Target already protected") and the
            // heap, producing non-deterministic munmap_chunk aborts at random
            // simulation times. The data extraction above stays parallel.
            #pragma omp critical (openwq_hdf5_write)
            {
                auto _cit = OpenWQ_wqconfig.conc_dsets.find(filename);
                auto _tit = OpenWQ_wqconfig.time_dsets.find(filename);
                hid_t _cds = (_cit != OpenWQ_wqconfig.conc_dsets.end()) ? _cit->second : -1;
                hid_t _tds = (_tit != OpenWQ_wqconfig.time_dsets.end()) ? _tit->second : -1;
                appendData_to_HDF5_file(_cds, _tds, data2print, timestr);
            }
        }
    }

    return EXIT_SUCCESS;
}

int OpenWQ_output::writeHDF5_Sediment(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    std::unique_ptr<arma::cube>& OpenWQ_var2print,
    std::string& output_file_label,
    std::string timestr,
    int icmp){

    // Get number of cells to print
    const unsigned int num_cells2print = OpenWQ_wqconfig.cells2print_vec.at(icmp).n_rows;
    const bool printflag = OpenWQ_wqconfig.cells2print_bool.at(icmp);

    // Early exit if nothing to print
    if (!printflag || num_cells2print == 0) return EXIT_SUCCESS;

    // Cache values
    // Sediment is exported as a bulk MASS (sedmass) — unlike the chemistry output
    // it is NOT divided by water volume and NOT no-water-flagged, so the concentration
    // helpers (is_conc_requested / watervol_minlim / noWaterConc / unit multipliers)
    // are intentionally not used here.
    const std::string CompName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);

    // Build filename
    std::string filename = OpenWQ_wqconfig.get_output_dir() + "/" +
                          CompName_icmp + "@Sediment-" + output_file_label + ".h5";

    /* ########################################
    // Save metadata on first timestep
    ######################################## */
    if (OpenWQ_wqconfig.print_oneStep){

        // Prepare xyz elements
        arma::mat cells2print_xyzElements = OpenWQ_wqconfig.cells2print_vec[icmp];
        cells2print_xyzElements.for_each([](arma::mat::elem_type& val) { val += 1.0; });

        // Prepare domain size
        arma::mat cells2print_xyzElements_size(1, 3);
        cells2print_xyzElements_size(0, 0) = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        cells2print_xyzElements_size(0, 1) = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        cells2print_xyzElements_size(0, 2) = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

        // Write host IDs
        std::vector<std::string> host_ids;
        host_ids.reserve(num_cells2print);

        for (unsigned int ixyz = 0; ixyz < num_cells2print; ixyz++){
            const int ix = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 0));
            const int iy = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 1));
            const int iz = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 2));
            host_ids.push_back(OpenWQ_hostModelconfig.get_cellid_to_wq_at(icmp, ix, iy, iz));
        }

        std::vector<const char*> cstrs(host_ids.size());
        for (size_t i = 0; i < host_ids.size(); ++i) {
            cstrs[i] = host_ids[i].c_str();
        }
        const std::string internal_db_name = OpenWQ_hostModelconfig.get_cellid_to_wqlabel();

        // Create file with latest format (fractal heaps, not B-trees)
        hid_t fcpl_s = H5Pcreate(H5P_FILE_CREATE);
        hid_t fapl_s = H5Pcreate(H5P_FILE_ACCESS);
        // Locking off (see note in writeHDF5); default format + metadata cache.
        H5Pset_file_locking(fapl_s, 0, 0);

        hid_t file_h = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, fcpl_s, fapl_s);
        if (file_h >= 0){

            // Write xyz_elements
            {
                hsize_t dims[2] = {cells2print_xyzElements.n_cols,
                                   cells2print_xyzElements.n_rows};
                hid_t sp = H5Screate_simple(2, dims, NULL);
                hid_t ds = H5Dcreate(file_h, "xyz_elements", H5T_NATIVE_DOUBLE,
                                     sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (ds >= 0) {
                    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, cells2print_xyzElements.memptr());
                    H5Dclose(ds);
                }
                H5Sclose(sp);
            }

            // Write xyz_elements_size
            {
                hsize_t dims[2] = {cells2print_xyzElements_size.n_cols,
                                   cells2print_xyzElements_size.n_rows};
                hid_t sp = H5Screate_simple(2, dims, NULL);
                hid_t ds = H5Dcreate(file_h, "xyz_elements_size", H5T_NATIVE_DOUBLE,
                                     sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (ds >= 0) {
                    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, cells2print_xyzElements_size.memptr());
                    H5Dclose(ds);
                }
                H5Sclose(sp);
            }

            // Write host IDs as variable-length strings
            if (!host_ids.empty()){
                hsize_t dims[1] = {host_ids.size()};
                hid_t space = H5Screate_simple(1, dims, NULL);
                hid_t dtype = H5Tcopy(H5T_C_S1);
                H5Tset_size(dtype, H5T_VARIABLE);
                hid_t dset = H5Dcreate(file_h, internal_db_name.c_str(), dtype, space,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset >= 0){
                    H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cstrs.data());
                    H5Dclose(dset);
                }
                H5Tclose(dtype);
                H5Sclose(space);
            }

            // [REDESIGN] create the time-extensible datasets that each timestep
            // appends to (/concentrations + /timestamps)
            create_timeseries_datasets(file_h, num_cells2print, OpenWQ_wqconfig, filename);

            // Store handle — file stays open until destructor
            OpenWQ_wqconfig.files[filename] = file_h;
        }

        H5Pclose(fcpl_s);
        H5Pclose(fapl_s);

        std::string msg_string = "<OpenWQ> Created sediment output file " + filename;
        ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    /* ########################################
    // Write data for current timestep.
    // File is already open in OpenWQ_wqconfig.files.
    ######################################## */

    // Look up persistent handle (read-only map access)
    auto it = OpenWQ_wqconfig.files.find(filename);
    if (it == OpenWQ_wqconfig.files.end() || it->second < 0) {
        std::string msg = "<OpenWQ> WARNING: Sediment HDF5 handle not found, skipping: " + filename;
        ConsoleLog(OpenWQ_wqconfig, msg, true, true);
        return EXIT_SUCCESS;
    }

    arma::mat data2print(num_cells2print, 1);

    // Extract data
    for (unsigned int celli = 0; celli < num_cells2print; celli++){

        const unsigned int ix = OpenWQ_wqconfig.cells2print_vec[icmp](celli, 0);
        const unsigned int iy = OpenWQ_wqconfig.cells2print_vec[icmp](celli, 1);
        const unsigned int iz = OpenWQ_wqconfig.cells2print_vec[icmp](celli, 2);

        // Sediment is a bulk MASS (sedmass) — write it directly. Do NOT divide by
        // water volume or apply the no-water flag (the concentration logic used by
        // the chemistry output): that would mask sediment to the no-water value
        // whenever the transport compartment (e.g. RUNOFF) holds no stored water.
        data2print(celli, 0) = (*OpenWQ_var2print)(ix, iy, iz);
    }

    {
        auto _cit = OpenWQ_wqconfig.conc_dsets.find(filename);
        auto _tit = OpenWQ_wqconfig.time_dsets.find(filename);
        hid_t _cds = (_cit != OpenWQ_wqconfig.conc_dsets.end()) ? _cit->second : -1;
        hid_t _tds = (_tit != OpenWQ_wqconfig.time_dsets.end()) ? _tit->second : -1;
        appendData_to_HDF5_file(_cds, _tds, data2print, timestr);
    }

    return EXIT_SUCCESS;
}

/* ########################################
// Append one timestep to the HDF5 file.
// [REDESIGN] Instead of creating a NEW dataset per timestep (which exploded the
// group's link metadata for long runs and corrupted on weak filesystems), this
// appends one row to the pre-created extensible datasets:
//   data (num_cells x 1) -> new row of /concentrations
//   name (timestamp str) -> new element of /timestamps
######################################## */
bool OpenWQ_output::appendData_to_HDF5_file(
    hid_t conc_ds,
    hid_t time_ds,
    arma::mat& data,
    std::string name){

    // conc_ds / time_ds are the per-file dataset handles kept open for the whole
    // run (no per-step H5Dopen/H5Dclose), so the active chunk stays in the chunk
    // cache and is flushed once per full chunk instead of read-modify-written
    // every timestep.
    if (conc_ds < 0) return false;
    const hsize_t ncells = (hsize_t) data.n_rows;   // data is (num_cells x 1)

    // -------- append a row to /concentrations --------
    hid_t fspace = H5Dget_space(conc_ds);
    hsize_t cur[2] = {0, 0};
    H5Sget_simple_extent_dims(fspace, cur, NULL);
    H5Sclose(fspace);

    hsize_t newdims[2] = {cur[0] + 1, ncells};
    if (H5Dset_extent(conc_ds, newdims) < 0) return false;

    fspace = H5Dget_space(conc_ds);
    hsize_t start[2] = {cur[0], 0};
    hsize_t count[2] = {1, ncells};
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t mspace = H5Screate_simple(2, count, NULL);
    H5Dwrite(conc_ds, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, data.mem);
    H5Sclose(mspace);
    H5Sclose(fspace);

    // -------- append the timestamp to /timestamps --------
    if (time_ds >= 0){
        hid_t tfs = H5Dget_space(time_ds);
        hsize_t tcur[1] = {0};
        H5Sget_simple_extent_dims(tfs, tcur, NULL);
        H5Sclose(tfs);

        hsize_t tnew[1] = {tcur[0] + 1};
        H5Dset_extent(time_ds, tnew);

        tfs = H5Dget_space(time_ds);
        hsize_t tstart[1] = {tcur[0]};
        hsize_t tcount[1] = {1};
        H5Sselect_hyperslab(tfs, H5S_SELECT_SET, tstart, NULL, tcount, NULL);
        hid_t tms = H5Screate_simple(1, tcount, NULL);
        hid_t stype = H5Tcopy(H5T_C_S1);
        H5Tset_size(stype, H5T_VARIABLE);
        const char* sptr = name.c_str();
        H5Dwrite(time_ds, stype, tms, tfs, H5P_DEFAULT, &sptr);
        H5Tclose(stype);
        H5Sclose(tms);
        H5Sclose(tfs);
    }

    return true;
}