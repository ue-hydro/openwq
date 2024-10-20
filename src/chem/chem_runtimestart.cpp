
// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
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

#include "OpenWQ_chem.hpp"


/* #################################################
// Compute chemical transformations
################################################# */
void OpenWQ_chem::Run(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){
    

    // Loop over number of compartments
    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){
    
        BGC_Transform( // calls exprtk: parsing expression
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_output,
            icmp);

    }

    // Turn off printing of exception err message (only print on first step)
    OpenWQ_wqconfig.BGC_Transform_print_errmsg = false;
    OpenWQ_wqconfig.invalid_bgc_entry_errmsg = false;
}

/* #################################################
// Compute each chemical transformation
################################################# */
void OpenWQ_chem::BGC_Transform(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    unsigned int icmp){
    
    // Local variables
    std::string chemname; // chemical name
    unsigned int index_cons,index_prod; // indexed for consumed and produced chemical
    unsigned int num_transf; // number of transformation in biogeochemical cycle bgci
    unsigned int num_BGCcycles; // number of BGC frameworks in comparment icmp
    unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // number of cell in x-direction
    unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // number of cell in y-direction
    unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // number of cell in z-direction
    std::string BGCcycles_name_icmp,BGCcycles_name_list; // BGC Cycling name i of compartment icmp  
    double transf_mass;         // interactive mass to transfer between species
    std::string msg_string;             // error/warning message string

    // Local variables for expression evaluator: exprtk
    std::string expression_string; // expression string


    // Find compartment icmp name from code (host hydrological model)
    std::string CompName_icmp = 
        OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp); // Compartment icomp name

    // Get cycling_frameworks for compartment icomp (name = CompName_icmp) 
    // (compartment names need to match)
    try{
        std::vector<std::string> BGCcycles_icmp =             
            OpenWQ_json.Config["BIOGEOCHEMISTRY_CONFIGURATION"]
                [CompName_icmp]["CYCLING_FRAMEWORK"];

        // Get number of BGC frameworks in comparment icmp
        num_BGCcycles = BGCcycles_icmp.size();
            
        /* ########################################
        // Loop over biogeochemical cycling frameworks
        ######################################## */

        for (unsigned int bgci=0;bgci<num_BGCcycles;bgci++){

            // Get framework name i of compartment icmp
            BGCcycles_name_icmp = BGCcycles_icmp.at(bgci);

            // Get number transformations in biogeochemical cycle BGCcycles_name_icmp
            std::vector<unsigned int> transf_index;
            for (unsigned int index_j=0;index_j<OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_info.size();index_j++){
                BGCcycles_name_list = std::get<0>(OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_info.at(index_j)); // get BGC cycle from OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_info
                if(BGCcycles_name_list.compare(BGCcycles_name_icmp) == 0){
                    transf_index.push_back(index_j);
                }
            }
            
            num_transf = transf_index.size(); // number of transformation for this BGC cycle (as identified in openWQ_BGCnative_BGCexpressions_info)

            // Check if found any valid BBC cycling framework (based on the json BCG)
            if (num_transf == 0 && OpenWQ_wqconfig.invalid_bgc_entry_errmsg == true){

                // Create Message
                msg_string = "<OpenWQ> Unkown CYCLING_FRAMEWORK with name '"
                    + BGCcycles_name_icmp
                    + "' defined for compartment: " 
                    + CompName_icmp;

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(
                    OpenWQ_wqconfig,    // for Log file name
                    msg_string,         // message
                    true,               // print in console
                    true);              // print in log file

            }
                        
            /* ########################################
            // Loop over transformations in biogeochemical cycle bgci
            ######################################## */
            
            for (unsigned int transi=0;transi<num_transf;transi++){

                index_cons = std::get<3>(OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_info.at(transf_index[transi]));
                index_prod = std::get<4>(OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_info.at(transf_index[transi]));

                // Get indexes of chemicals in transformation equation (needs to be here for loop reset)
                std::vector<unsigned int> index_chemtransf = 
                    std::get<5>(
                        OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_info.at(transf_index[transi]));


                /* ########################################
                // Loop over space: nx, ny, nz
                ######################################## */
                for (unsigned int ix=0;ix<nx;ix++){
                    for (unsigned int iy=0;iy<ny;iy++){
                        for (unsigned int iz=0;iz<nz;iz++){                    
                            
                            //std::vector<double> openWQ_BGCnative_chemass_InTransfEq; // chemical mass involved in transformation (needs to be here for loop reset)
                            OpenWQ_wqconfig.openWQ_BGCnative_chemass_InTransfEq.clear();
                            // loop to get all the variables inside the expression
                            for (unsigned int chem=0;chem<index_chemtransf.size();chem++){
                                OpenWQ_wqconfig.openWQ_BGCnative_chemass_InTransfEq.push_back(
                                    (*OpenWQ_vars.chemass)
                                    (icmp)
                                    (index_chemtransf[chem])
                                    (ix,iy,iz));}

                            // Update current dependencies for current x, y and z
                            // dependVar_scalar needed for exportk
                            for (unsigned int depi=0;depi<OpenWQ_hostModelconfig.get_num_HydroDepend();depi++){
                                OpenWQ_hostModelconfig.set_dependVar_scalar_at(depi, OpenWQ_hostModelconfig.get_dependVar_at(depi,ix,iy,iz));}

                            // Mass transfered: Consumed -> Produced (using exprtk)
                            // Make sure that transf_mass is positive, otherwise ignore transformation
                            transf_mass = std::fmax(
                                OpenWQ_wqconfig.openWQ_BGCnative_BGCexpressions_eq[transf_index[transi]].value(),
                                0.0f); // Needs to be positive (from consumed to produced)

                            // Guarantee that removed mass is not larger than existing mass
                            if (!OpenWQ_hostModelconfig.is_time_step_0()){
                                transf_mass = std::fmin(
                                        (*OpenWQ_vars.chemass)(icmp)(index_cons)(ix,iy,iz),
                                        transf_mass * OpenWQ_hostModelconfig.get_time_step());    // need to multiply by the time step
                                                                                            // because that's the total mass that is 
                                                                                            // going to be added/substracted in the 
                                                                                            // solver
                            }else{
                                transf_mass = 0.0f;
                            }
                            // New mass of consumed chemical
                            (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(index_cons)(ix,iy,iz) -= transf_mass;

                            // New mass of produced chemical
                            (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(index_prod)(ix,iy,iz) += transf_mass;

                        }
                    }
                }
            }
        }
    }catch(json::exception& e){

        // Print error message if 1st time step
        if (OpenWQ_wqconfig.BGC_Transform_print_errmsg == true){

            // Create Message
            msg_string = 
                "<OpenWQ> No CYCLING_FRAMEWORK defined for compartment: "
                + CompName_icmp;

            // Print it (Console and/or Log file)
            OpenWQ_output.ConsoleLog(
                OpenWQ_wqconfig,    // for Log file name
                msg_string,         // message
                true,               // print in console
                true);              // print in log file
        }
        
        return;
    } 
}