
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

// ###########################################################################
// Resolve a cell spec -> explicit (ix,iy,iz) cell list, by HOST cell_id.
// Shared by COMPARTMENTS_AND_CELLS and FLUXES_CONC_TO_PRINT so BOTH select
// cells the SAME way SS/EWF do: via find_indices_from_cellid (reachID/hruId),
// NOT by raw array position.
//   spec = ["all", ...]               -> every cell of the compartment
//        = ["<id>", ...]              -> that reach/HRU by its host id
//        = [ ["<id1>","<id2>"], ... ] -> several, by host id
// (iy/iz come from the resolution; for 1-D reach/HRU compartments they're "all").
// ###########################################################################
static arma::mat _resolve_cellspec(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    int comp_index,
    int nx, int ny, int nz,
    const nlohmann::json& spec,
    const std::string& ctx_label)
{
    auto _as_id = [](const nlohmann::json& v, std::string& s)->bool{
        if (v.is_string()){
            s = v.get<std::string>();
            std::string u = s; std::transform(u.begin(), u.end(), u.begin(), ::toupper);
            return u.compare("ALL") != 0;                 // false if "all"
        }
        if (v.is_number_integer()){ s = std::to_string(v.get<long long>()); return true; }
        if (v.is_number()){ s = std::to_string(v.get<double>()); return true; }
        return false;
    };

    std::vector<std::string> ids;
    bool use_all = true;
    if (spec.is_array() && spec.size() >= 1){
        const nlohmann::json& slot0 = spec.at(0);
        if (slot0.is_array()){
            for (auto& e : slot0){ std::string s; if (_as_id(e, s)) ids.push_back(s); }
        } else { std::string s; if (_as_id(slot0, s)) ids.push_back(s); }
    }
    if (!ids.empty()) use_all = false;

    std::vector<std::array<int,3>> sel;
    if (use_all){
        for (int ix=0; ix<nx; ix++)
        for (int iy=0; iy<ny; iy++)
        for (int iz=0; iz<nz; iz++)
            sel.push_back({ix, iy, iz});
    } else {
        std::string msg;
        for (auto& id0 : ids){
            int rx=-1, ry=-1, rz=-1; bool pm=false;
            if (OpenWQ_hostModelconfig.find_indices_from_cellid(comp_index, id0, rx, ry, rz, pm)){
                sel.push_back({rx, ry, rz});
            } else {
                msg = "<OpenWQ> WARNING: cell_id '" + id0 + "' (from "
                    + OpenWQ_hostModelconfig.get_cellid_to_wqlabel() + ") not found for "
                    + ctx_label + " (entry skipped).";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg, true, true);
            }
        }
    }
    arma::mat cells((arma::uword)sel.size(), 3);
    for (arma::uword k=0; k<sel.size(); k++){
        cells(k,0)=sel[k][0]; cells(k,1)=sel[k][1]; cells(k,2)=sel[k][2];
    }
    return cells;
}

// Set what to print
void OpenWQ_readjson::SetConfigInfo_output_what2print(
    OpenWQ_json &OpenWQ_json,
    OpenWQ_hostModelconfig & OpenWQ_hostModelconfig,
    OpenWQ_wqconfig &OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units &OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    int nx;                                                 // interactive nx inforationfor each compartment
    int ny;                                                 // interactive ny inforationfor each compartment
    int nz;                                                 // interactive nz inforationfor each compartment
    int ix_json;                                            // iteractive ix info for COMPARTMENTS_AND_CELLS cell data 
    int iy_json;                                            // iteractive iy info for COMPARTMENTS_AND_CELLS cell data  
    int iz_json;                                            // iteractive iz info for COMPARTMENTS_AND_CELLS cell data  
    int spX_min, spX_max, spY_min, spY_max, spZ_min, spZ_max; // range of ix, iy, iz of cells to print (if selected "all")
    long num_compt2print;                                   // number of compartments to print
    std::string chem_name2print;                            // iteractive chem name from openWQ_OUTPUT
    std::string chem_namelist;                              // iteractive chem name from BGC list
    long num_chem2print;                                    // number of chemicals to print based on openWQ_OUTPUT               
    long num_cells2print;                                   // iteractive number of cells to print for each compartment
    std::vector<std::string> compt_names_vec;               // compartment names to print  
    std::string CompName_icmp;                              // iteractive compartment name from list
    std::string compt_name2print;                           // iteractive compartment name to print
    std::string cells_input;                                // interactive string for cells input for each compartment
    arma::mat cells2print_cmpt;                             // cumulative vector with all cells to print for each compartment
    bool noValPrintRequest_flag;                            // flag to through message if no valid compartments/cells have been selected for printing
    std::vector<int>::iterator it;                          // iteractor for flagging if compartment i has been selected for printing
    unsigned long long iRow;                                // dummy variable to track the number of rows for cell2print
    unsigned long long nRows;                               // dummy variable to store the total number of rows for cell2print
    std::string msg_string;
    json json_output_subStruct;
    json json_output_subStruct_CmpCells;
    std::string errorMsgIdentifier;

    // Chemicals to print
    // num of chem to print
    errorMsgIdentifier = "Master file";
    json_output_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.Master,"OPENWQ_OUTPUT",
        errorMsgIdentifier,
        true);

    errorMsgIdentifier = "Master file > OPENWQ_OUTPUT";
    json json_chemSpecies = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        json_output_subStruct,"CHEMICAL_SPECIES",
        errorMsgIdentifier,
        true);
    num_chem2print = json_chemSpecies.size();

    // Get indexes for the list of chemicals requested
    // CHEMICAL_SPECIES is now an array of species name strings (e.g., ["NO3-N", "NH4-N"])
    // Resolve each name to its 0-based index in the internal species list
    const auto& chem_list = *OpenWQ_wqconfig.cached_chem_species_list_ptr;
    const unsigned int num_chem_total = OpenWQ_wqconfig.cached_num_chem;
    for (unsigned int chemi = 0; chemi < num_chem2print; chemi++){
        // Chemical name from the JSON array
        chem_name2print = json_chemSpecies.at(chemi).get<std::string>();
        bool found = false;
        for (unsigned int chemlisti = 0; chemlisti < num_chem_total; chemlisti++){
            chem_namelist = chem_list[chemlisti];
            // Case-insensitive match: the internal species list can hold
            // mixed-case names (e.g. PHREEQC components "Amm", "Ca", "Charge")
            // while the requested output name is upper-cased by the JSON loader,
            // so a case-sensitive compare would silently drop them.
            bool namematch = (chem_namelist.size() == chem_name2print.size());
            for (size_t _ci = 0; namematch && _ci < chem_namelist.size(); ++_ci){
                namematch = (std::tolower((unsigned char)chem_namelist[_ci])
                          == std::tolower((unsigned char)chem_name2print[_ci]));
            }
            if (namematch){
                OpenWQ_wqconfig.chem2print.push_back(chemlisti);
                found = true;
                break;
            }
        }
        if (!found){
            msg_string =
                "<OpenWQ> WARNING: CHEMICAL_SPECIES name '"
                + chem_name2print
                + "' not found in internal species list. Skipping.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        }
    }

    // ########################################
    // Compartments to print

    errorMsgIdentifier = "Master file > OUTPUT";
    json_output_subStruct_CmpCells = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        json_output_subStruct,"COMPARTMENTS_AND_CELLS",
        errorMsgIdentifier,
        true);

    // num of comartments to print
    num_compt2print = json_output_subStruct_CmpCells.size();
    
    // Get list of compartments to print from COMPARTMENTS_AND_CELLS
    for (auto& x1 : json_output_subStruct_CmpCells.items()){
                
        try{
            compt_names_vec.push_back(x1.key());
        }catch (...){}
    }

    // Get indexes for the list of compartments requested
    for (unsigned int cmpti = 0; cmpti < num_compt2print; cmpti++){
        compt_name2print = compt_names_vec[cmpti]; 
        for (unsigned int icmp = 0; icmp < OpenWQ_hostModelconfig.get_num_HydroComp(); icmp++){     
            CompName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);  // name
            // Check if compartments listed match internal compartment names
            if (CompName_icmp.compare(compt_name2print) == 0){                               
                OpenWQ_wqconfig.compt2print.push_back(icmp);
                break;
            }
        }
    }

    // ########################################
    // Flux-concentration exports to print (OPTIONAL -> backward compatible).
    // SAME structure as COMPARTMENTS_AND_CELLS: an object keyed by flux name,
    // each value a numbered set of cell entries (by host cell_id):
    //   "FLUXES_CONC_TO_PRINT": {
    //       "RUNOFF_TO_STREAM": { "1": ["all","all","all"] },
    //       "REACH_OUTFLOW":    { "1": ["78830302","all","all"],
    //                             "2": ["78830325","all","all"] }
    //   }
    // Each entry is resolved by the shared _resolve_cellspec (reachID/hruId).
    // Also accepts an array shorthand ["RUNOFF_TO_STREAM", ...] = all cells.
    // Names matched case-insensitively against the host-registered flux exports.
    if (json_output_subStruct.contains("FLUXES_CONC_TO_PRINT")){

        // Reset first so a re-parse of the config cannot double-populate the
        // list (the flux write loop iterates this vector directly, unlike the
        // compartment loop which iterates unique compartments).
        OpenWQ_wqconfig.fluxconc2print.clear();
        OpenWQ_wqconfig.fluxconc_cells2print.clear();

        nlohmann::json& _fxnode = json_output_subStruct["FLUXES_CONC_TO_PRINT"];

        auto _process_flux = [&](std::string flux_name2print,
                                 const nlohmann::json& _entries){
            std::transform(flux_name2print.begin(), flux_name2print.end(),
                           flux_name2print.begin(), ::toupper);
            // Match to a registered flux export CASE-INSENSITIVELY: flux exports
            // are named after host-model variables, which may be mixed-case
            // (e.g. "averageRoutedRunoff"), so compare uppercased-vs-uppercased.
            // The registered (mixed-case) name is still what's used for the h5
            // filename, so the host variable name is preserved verbatim on disk.
            int _ifx = -1;
            for (unsigned int ifx = 0;
                 ifx < OpenWQ_hostModelconfig.get_num_FluxConcExport(); ifx++){
                std::string _reg_uc = OpenWQ_hostModelconfig.get_FluxConcExport_name_at(ifx);
                std::transform(_reg_uc.begin(), _reg_uc.end(), _reg_uc.begin(), ::toupper);
                if (_reg_uc.compare(flux_name2print) == 0){ _ifx = (int)ifx; break; }
            }
            if (_ifx < 0){
                msg_string = "<OpenWQ> WARNING: FLUXES_CONC_TO_PRINT name '"
                    + flux_name2print + "' not found in host flux-export list. Skipping.";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                return;
            }
            const int _src = OpenWQ_hostModelconfig.get_FluxConcExport_srcComp_at(_ifx);
            const int _nx = (int)OpenWQ_hostModelconfig.get_FluxConcExport_num_cells_x_at(_ifx);
            const int _ny = (int)OpenWQ_hostModelconfig.get_FluxConcExport_num_cells_y_at(_ifx);
            const int _nz = (int)OpenWQ_hostModelconfig.get_FluxConcExport_num_cells_z_at(_ifx);
            const std::string _ctx = "FLUXES_CONC_TO_PRINT > " + flux_name2print;

            // Accumulate the cells from each numbered entry (via the shared
            // resolver -- same pattern as COMPARTMENTS_AND_CELLS + SS/EWF).
            arma::mat _cells;
            if (_entries.is_object()){
                for (auto& _ent : _entries.items()){
                    arma::mat _c = _resolve_cellspec(
                        OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_output,
                        _src, _nx, _ny, _nz, _ent.value(), _ctx);
                    if (!_c.is_empty()) _cells.insert_rows(_cells.n_rows, _c);
                }
            } else {
                // array-shorthand name -> all cells
                _cells = _resolve_cellspec(
                    OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_output,
                    _src, _nx, _ny, _nz, nlohmann::json(nullptr), _ctx);
            }
            OpenWQ_wqconfig.fluxconc2print.push_back(_ifx);
            OpenWQ_wqconfig.fluxconc_cells2print[_ifx] = _cells;
        };

        if (_fxnode.is_array()){
            for (auto& fx : _fxnode){
                try { _process_flux(fx.get<std::string>(), nlohmann::json(nullptr)); }
                catch (...) {}
            }
        } else if (_fxnode.is_object()){
            for (auto& kv : _fxnode.items()) _process_flux(kv.key(), kv.value());
        }
    }

    // ########################################
    // Cells to print

    // Set noValPrintRequest_flag to default false
    // if no valid request are found, then this will remain false and a warning message will be sent
    // otherwise, this will change to true and all is good
    noValPrintRequest_flag = false;

    // Check what cells to print for each viable (compt2print) compartment requested
    for (unsigned int icmp = 0; icmp < OpenWQ_hostModelconfig.get_num_HydroComp(); icmp++){

        // Check if icmp is included in compt2print
        it = find (
            OpenWQ_wqconfig.compt2print.begin(), 
            OpenWQ_wqconfig.compt2print.end(),
            icmp);

        // if icmp is not in compt2print, then print dummy_cell and skip
        if (it == OpenWQ_wqconfig.compt2print.end()){
            // set cells2print_vec to be null because compartment has not been selected
            OpenWQ_wqconfig.cells2print_bool.push_back(false);
            OpenWQ_wqconfig.cells2print_vec.push_back(cells2print_cmpt); 
            cells2print_cmpt.clear();
            continue;
        }

        // compartnment icmp domain
        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);
        
        // first set default no print
        // this gets true if loading of json input is sucessfull
        OpenWQ_wqconfig.cells2print_bool.push_back(false);

        // Get number of cell entries provided
        num_cells2print = json_output_subStruct_CmpCells
            [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)].size();

        // Loop over cell entries; resolve each by host cell_id via the shared
        // resolver (same pattern as FLUXES_CONC_TO_PRINT + SS/EWF).
        for (unsigned int celli = 0; celli < num_cells2print; celli++){
            const nlohmann::json& _entry = json_output_subStruct_CmpCells
                [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                [std::to_string(celli + 1)];
            arma::mat _entry_cells = _resolve_cellspec(
                OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_output,
                (int)icmp, nx, ny, nz, _entry,
                "COMPARTMENTS_AND_CELLS > "
                    + OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp));
            if (!_entry_cells.is_empty())
                cells2print_cmpt.insert_rows(cells2print_cmpt.n_rows, _entry_cells);
        }

        // Add the complete list of cells to print
        OpenWQ_wqconfig.cells2print_vec.push_back(cells2print_cmpt);

         // all went well in "try", so set print to true
         // and noValPrintRequest_flag to true
        if(cells2print_cmpt.is_empty() == false){
            OpenWQ_wqconfig.cells2print_bool[icmp] = true;
            noValPrintRequest_flag = true;
        }

        // Clear cells2print_cmpt for reuse
        cells2print_cmpt.clear();

    }

    // If no valid output request are found, then through a warning message
    if (noValPrintRequest_flag == false){

        // Create Message (Error - locate problematic cell)
        msg_string = 
            "<OpenWQ> WARNING: No 'COMPARTMENTS_AND_CELLS' requests have been found. "
            "The model will not print any results. "
            "Revise the 'OPENWQ_OUTPUT' keys in the master-json file and try again.";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);

    }

}