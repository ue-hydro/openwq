
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

#include "models_CH/headerfile_CH.hpp"
#include <unistd.h>  // for getcwd
#include <climits>   // for PATH_MAX

// Helper function to resolve relative paths
static std::string resolve_path(const std::string& path) {
    // If path is already absolute, return as-is
    if (!path.empty() && path[0] == '/') {
        return path;
    }
    // Get current working directory and prepend to relative path
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        return std::string(cwd) + "/" + path;
    }
    // Fallback: return original path
    return path;
}

/* #################################################
// Parse biogeochemical expressions
################################################# */
void OpenWQ_CH_model::phreeqc_setup(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
    int nx, ny, nz, nxyz;
    bool component_h2o;
    std::string  database, PHREEQC_filename;

    nxyz = 0;
    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++) {
        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
        nxyz += nx*ny*nz;
    }

    // Use single thread for PhreeqcRM to avoid threading issues in calc_dens()
    // TODO: investigate multi-threaded PhreeqcRM crashes
    int nthreads = 1;  // OpenWQ_wqconfig.get_num_threads_requested();
    std::string msg_string = "<OpenWQ> PHREEQC: Creating PhreeqcRM with nxyz=" + std::to_string(nxyz)
        + ", nthreads=" + std::to_string(nthreads);
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm = std::unique_ptr<PhreeqcRM>(new PhreeqcRM(nxyz, nthreads));


    try {
        component_h2o = OpenWQ_json.BGC_module["COMPONENT_H2O"];
    } catch (...){
        component_h2o = true;
    }

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetComponentH2O(component_h2o);

    // Set error handler mode (0 = return error codes - safest for debugging)
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetErrorHandlerMode(0);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetFilePrefix("phreeqc");
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->OpenFiles();

    database = OpenWQ_json.BGC_module["DATABASE"];
    std::string database_fullpath = resolve_path(database);
    msg_string = "<OpenWQ> PHREEQC: Loading database from: " + database_fullpath;
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    IRM_RESULT status = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->LoadDatabase(database_fullpath);
    if (status != IRM_OK) {
        msg_string = "<OpenWQ> ERROR: PHREEQC LoadDatabase failed for: " + database_fullpath
            + " (error code: " + std::to_string(status) + ")";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        exit(EXIT_FAILURE);
    }
    msg_string = "<OpenWQ> PHREEQC: Database loaded successfully";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    PHREEQC_filename = OpenWQ_json.BGC_module["FILEPATH"];
    std::string phreeqc_fullpath = resolve_path(PHREEQC_filename);
    msg_string = "<OpenWQ> PHREEQC: Running input file: " + phreeqc_fullpath;
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    status = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunFile(true, true, true, phreeqc_fullpath);
    if (status != IRM_OK) {
        msg_string = "<OpenWQ> ERROR: PHREEQC RunFile failed for: " + PHREEQC_filename
            + " (error code: " + std::to_string(status) + ")";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        exit(EXIT_FAILURE);
    }
    msg_string = "<OpenWQ> PHREEQC: Input file executed successfully";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Set units BEFORE FindComponents
    // Units: 1 = mg/kgw, 2 = mol/kgw, 3 = kg/kgs
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetUnitsSolution(2);  // 2 = mol/kgw (PHREEQC default)

    // FindComponents must be called on the production PhreeqcRM instance
    int ncomp = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->FindComponents();
    msg_string = "<OpenWQ> PHREEQC FindComponents returned: " + std::to_string(ncomp);
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    int num_chem_check = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponentCount();
    msg_string = "<OpenWQ> PHREEQC GetComponentCount: " + std::to_string(num_chem_check);
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    if (num_chem_check > 0 && (unsigned int)num_chem_check != OpenWQ_wqconfig.CH_model->PHREEQC->num_chem) {
        msg_string =
            "<OpenWQ> WARNING: PHREEQC production grid has "
            + std::to_string(num_chem_check)
            + " components, but "
            + std::to_string(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem)
            + " were found during config. Updating to production value.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        OpenWQ_wqconfig.CH_model->PHREEQC->num_chem = (unsigned int)num_chem_check;
    }

    // Set representative volume and saturation for all cells
    std::vector<double> rv(nxyz, 1.0);  // 1 liter per cell
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetRepresentativeVolume(rv);

    std::vector<double> sat(nxyz, 1.0);  // Fully saturated
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetSaturationUser(sat);

    // Set fixed density to avoid calc_dens() issues with some databases
    std::vector<double> density(nxyz, 1.0);  // 1.0 kg/L (water density)
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetDensityUser(density);

    // CRITICAL: Disable solution density/volume calculations to avoid calc_dens() crash
    // false = use values from SetDensityUser instead of calculating internally
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->UseSolutionDensityVolume(false);

    // Log component names
    std::vector<std::string> components = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponents();
    msg_string = "<OpenWQ> PHREEQC components:";
    for (size_t i = 0; i < components.size(); i++) {
        msg_string += " [" + std::to_string(i) + "]=" + components[i];
    }
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

}
