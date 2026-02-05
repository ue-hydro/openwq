
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
            // Check if species name matches internal species list
            if (chem_namelist.compare(chem_name2print) == 0){
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

        // Loop over all cells
        for (unsigned int celli = 0; celli < num_cells2print; celli++){
            
            // Get ix value (as integer)
            try{
                ix_json = json_output_subStruct_CmpCells
                    [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                    [std::to_string(celli + 1)]
                    .at(0);
                ix_json --; // remove 1 to match c++ convention to start in zero
            
            }catch(...){
                // if not integer, then check if "all"
                cells_input = json_output_subStruct_CmpCells
                    [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                    [std::to_string(celli + 1)]
                    .at(0);

                if (cells_input.compare("ALL")==0){
                    // then, use "all" flag (=-1)
                    ix_json = -1;

                }else{
                    msg_string = 
                        "<OpenWQ> WARNING: Unkown entry for ix (" 
                        + cells_input
                        + ") for OPENWQ_OUTPUT > COMPARTMENTS_AND_CELLS > " 
                        + OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)
                        + "(entry skipped)";

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);

                    // not a valid entry
                    continue;
                    
                }
            }

            // Get iy value
            try{
                iy_json = json_output_subStruct_CmpCells
                    [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                    [std::to_string(celli + 1)]
                    .at(1);
                iy_json--;  // remove 1 to match c++ convention to start in zero
            
            }catch(...){

                // if not integer, then check if "all"
                cells_input = json_output_subStruct_CmpCells
                    [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                    [std::to_string(celli + 1)]
                    .at(1);

                if (cells_input.compare("ALL")==0){
                    // then, use "all" flag (=-1)
                    iy_json = -1;

                }else{
                    msg_string = 
                        "<OpenWQ> WARNING: Unkown entry for iy (" 
                        + cells_input
                        + ") for OPENWQ_OUTPUT > COMPARTMENTS_AND_CELLS > " 
                        + OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)
                        + "(entry skipped)";

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);

                    // not a valid entry
                    continue;

                }
            }

            // Get iz value
            try{
                iz_json = json_output_subStruct_CmpCells
                    [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                    [std::to_string(celli + 1)]
                    .at(2);
                iz_json --; // remove 1 to match c++ convention to start in zero

            }catch(...){

                // if not integer, then check if "all"
                cells_input = json_output_subStruct_CmpCells
                    [OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)]
                    [std::to_string(celli + 1)]
                    .at(2);

                if (cells_input.compare("ALL")==0){
                    // then, use "all" flag (=-1)
                    iz_json = -1;

                }else{
                    msg_string = 
                        "<OpenWQ> WARNING: Unkown entry for iz (" 
                        + cells_input
                        + ") for OPENWQ_OUTPUT > COMPARTMENTS_AND_CELLS > " 
                        + OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)
                        + "(entry skipped)";

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);

                    // not a valid entry
                    continue;

                }
            }

            // Check if cell requested is witin the boundaries of the spatial domain
            // If yes, add the cell selected to cells2print_vec
            // Otherwise, write warning message (skip entry)
            if (ix_json <= nx
                && iy_json <= ny
                && iz_json <= nz){

                // ix
                if(ix_json != -1){spX_min = ix_json; spX_max = ix_json;}
                else{spX_min = 0; spX_max = nx - 1;}
                // iy
                if(iy_json != -1){spY_min = iy_json; spY_max = iy_json;}
                else{spY_min = 0; spY_max = ny - 1;}
                // iz
                if(iz_json != -1){spZ_min = iz_json; spZ_max = iz_json;}
                else{spZ_min = 0; spZ_max = nz - 1;}

                // Create arma mat to append
                nRows = (spX_max - spX_min + 1) * (spY_max - spY_min + 1) * (spZ_max - spZ_min + 1);
                arma::mat cells2print_row(nRows,3,arma::fill::zeros);       // iteractive vector with cell x, y and z indexes

                // Reset dummy variables
                iRow = 0;

                // add cell requested to list to print for each compartment
                // first create the vector cells2print_row with x, y and z values
                // loop is to account for "all" entries
                for (int ix_j = spX_min; ix_j <= spX_max; ix_j++){
                    for (int iy_j = spY_min; iy_j <= spY_max; iy_j++){
                        for (int iz_j = spZ_min; iz_j <= spZ_max; iz_j++){
                            
                            cells2print_row(iRow,0) = ix_j;
                            cells2print_row(iRow,1) = iy_j;
                            cells2print_row(iRow,2) = iz_j;
                            iRow++;

                        }
                    }
                }

                // add that cells2print_row to cells2print_cmpt
                cells2print_cmpt.insert_rows(
                    cells2print_cmpt.n_rows,
                    cells2print_row);

                // Reset dummy variables
                cells2print_row.zeros();

            }else{
                
                // Create Message (Error - locate problematic cell)
                msg_string = 
                    "<OpenWQ> ERROR: Cell entry provided out of domain"
                    " in OPENWQ_OUTPUT > COMPARTMENTS_AND_CELLS > "
                    + OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)
                    + "> '" + std::to_string(celli + 1) + "'";

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig,msg_string,true,true);

                // not a valid entry
                continue;

            }
                
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