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

                output_label = "d_output_dt_transport";
                writeCSV(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                        OpenWQ_vars.d_chemass_dt_transp_out, output_label, timestr, icmp);

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
            }

            // Debug outputs
            if (debug_mode){

                output_label = "d_output_dt_chemistry";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_dt_chem_out, output_label, timestr, icmp);

                output_label = "d_output_dt_transport";
                writeHDF5(OpenWQ_json, OpenWQ_hostModelconfig, OpenWQ_wqconfig,
                         OpenWQ_vars.d_chemass_dt_transp_out, output_label, timestr, icmp);

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
    const bool is_native_flex = 
        (OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0;
    
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
    const bool is_native_flex = 
        (OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0;

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
    // Save metadata on first timestep
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

        // Write metadata for each chemical
        for (unsigned int ichem = 0; ichem < num_chem2print; ichem++){

            const std::string chem_name = (*OpenWQ_wqconfig.cached_chem_species_list_ptr)[OpenWQ_wqconfig.chem2print[ichem]];

            // Build filename
            std::string filename = OpenWQ_wqconfig.get_output_dir() + "/" +
                                  CompName_icmp + "@" + chem_name + "#" +
                                  units_string + "-" + output_file_label + ".h5";

            // Save xyz elements
            cells2print_xyzElements.save(arma::hdf5_name(filename, "xyz_elements"));

            // Save domain size
            cells2print_xyzElements_size.save(
                arma::hdf5_name(filename, "xyz_elements_size", arma::hdf5_opts::append));

            // Write host IDs as variable-length strings
            if (!host_ids.empty()){
                const std::string internal_db_name = OpenWQ_hostModelconfig.get_cellid_to_wqlabel();
                
                std::vector<const char*> cstrs(host_ids.size());
                for (size_t i = 0; i < host_ids.size(); ++i) {
                    cstrs[i] = host_ids[i].c_str();
                }

                hid_t file_h = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                if (file_h >= 0){
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
                    H5Fclose(file_h);
                }
            }

            // Open file for subsequent writes
            OpenWQ_wqconfig.files[filename] = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

            std::string msg_string = "<OpenWQ> Opening output file " + filename +
                                    " " + std::to_string(OpenWQ_wqconfig.files[filename]);
            ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }

    /* ########################################
    // Write data for current timestep
    ######################################## */

    // Parallelize over chemicals
    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (unsigned int ichem = 0; ichem < num_chem2print; ichem++){

            const std::string chem_name = (*OpenWQ_wqconfig.cached_chem_species_list_ptr)[OpenWQ_wqconfig.chem2print[ichem]];

            // Build filename
            std::string filename = OpenWQ_wqconfig.get_output_dir() + "/" +
                                  CompName_icmp + "@" + chem_name + "#" +
                                  units_string + "-" + output_file_label + ".h5";

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

            // Write to HDF5 (thread-safe as each thread writes to different file)
            hid_t file = OpenWQ_wqconfig.files[filename];
            appendData_to_HDF5_file(file, data2print, timestr);
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
    const std::string CompName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
    const bool is_conc_requested = OpenWQ_wqconfig.is_concentration_requested();
    const double watervol_minlim = OpenWQ_hostModelconfig.get_watervol_minlim();
    const double noWaterConc = OpenWQ_wqconfig.noWaterConc;
    const double unit_mult_num = OpenWQ_wqconfig.get_output_units_numerator();
    const double unit_mult_den = OpenWQ_wqconfig.get_output_units_denominator();

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

        // Save metadata
        cells2print_xyzElements.save(arma::hdf5_name(filename, "xyz_elements"));
        cells2print_xyzElements_size.save(
            arma::hdf5_name(filename, "xyz_elements_size", arma::hdf5_opts::append));

        // Write host IDs
        std::vector<std::string> host_ids;
        host_ids.reserve(num_cells2print);

        for (unsigned int ixyz = 0; ixyz < num_cells2print; ixyz++){
            const int ix = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 0));
            const int iy = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 1));
            const int iz = static_cast<int>(OpenWQ_wqconfig.cells2print_vec[icmp](ixyz, 2));
            host_ids.push_back(OpenWQ_hostModelconfig.get_cellid_to_wq_at(icmp, ix, iy, iz));
        }

        if (!host_ids.empty()){
            const std::string internal_db_name = OpenWQ_hostModelconfig.get_cellid_to_wqlabel();
            
            std::vector<const char*> cstrs(host_ids.size());
            for (size_t i = 0; i < host_ids.size(); ++i) {
                cstrs[i] = host_ids[i].c_str();
            }

            hid_t file_h = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (file_h >= 0){
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
                H5Fclose(file_h);
            }
        }

        OpenWQ_wqconfig.files[filename] = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    /* ########################################
    // Write data for current timestep
    ######################################## */

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
                (*OpenWQ_var2print)(ix, iy, iz) *
                unit_mult_num / (water_vol * unit_mult_den);
        } else {
            data2print(celli, 0) = noWaterConc;
        }
    }

    // Write to HDF5
    hid_t file = OpenWQ_wqconfig.files[filename];
    appendData_to_HDF5_file(file, data2print, timestr);

    return EXIT_SUCCESS;
}

/* ########################################
// Append data to HDF5 file
######################################## */
bool OpenWQ_output::appendData_to_HDF5_file(
    hid_t file,
    arma::mat& data,
    std::string name){

    hsize_t dims[2] = {data.n_cols, data.n_rows};

    hid_t dataspace = H5Screate_simple(2, dims, NULL);
    hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    hid_t dataset = H5Dcreate(file, name.c_str(), datatype, dataspace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dataset < 0){
        H5Sclose(dataspace);
        H5Tclose(datatype);
        return false;
    }

    // Write data
    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.mem);

    // Cleanup
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Tclose(datatype);

    return true;
}