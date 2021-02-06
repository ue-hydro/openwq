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

#include <iostream>
#include <fstream>
#include <armadillo>
#include <string>

// #include "utility.h"

#include "OpenWQ_global.h"
#include "OpenWQ_load.h"
#include "OpenWQ_initiate.h"
#include "OpenWQ_chem.h"
#include "OpenWQ_watertransp.h"
//#include "OpenWQ_print.h"

int main(int argc, char* argv[]) 
{   
    std::string vtufilename;

    // Create Object: OpenWQ_hostModelconfig (link to host hydrological model)
    OpenWQ_hostModelconfig OpenWQ_hostModelconfig;
    // Link openWQ data strucuture indexes to hydrological model compartments 
    typedef std::tuple<int,std::string,int, int, int> hydroTuple;
    /* hydroTuple: 
    (1) index in openWQ variable, 
    (2) name identifier in openWQ_config.json, 
    (3) number of cell in x-direction
    (4) number of cell in y-direction
    (5) number of cell in z-direction
    */
    OpenWQ_hostModelconfig.HydroComp.push_back(hydroTuple(0,"Snow",100,50,1));
    OpenWQ_hostModelconfig.HydroComp.push_back(hydroTuple(1,"Soil",100,50,1));
    OpenWQ_hostModelconfig.HydroComp.push_back(hydroTuple(2,"Groundwater",100,50,1));
    OpenWQ_hostModelconfig.HydroComp.push_back(hydroTuple(3,"Streams",100,50,1));
    // (add other compartments as needed)...

    int num_HydroComp = OpenWQ_hostModelconfig.HydroComp.size(); // number of hydrological compartments in host model


    // Create Object: OpenWQ_json (Input JSON files)
    OpenWQ_json OpenWQ_json; // create object
    
    // Load input data
    OpenWQ_load OpenWQ_load; // create object: json files load modules
    OpenWQ_load.loadinit(
        OpenWQ_json);
    
    // Create Object: OpenWQ_vars (openWQ variables)
    OpenWQ_vars OpenWQ_vars(
        num_HydroComp);
   
    // OpenWQ_initiate
    OpenWQ_initiate OpenWQ_initiate; // create object: start modules e.g., initiate
    
    // Allocate memory: set variable sizes
    OpenWQ_initiate.initmemory(
        OpenWQ_json,
        OpenWQ_vars,
        OpenWQ_hostModelconfig);
    
    // Read and Set Initial Conditions 
    // Needs loop in host model because it requires grid-cell volume (mass IC) or water volume (concentration IC))
    int icmp = 0;                    // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    int ix = 1;                      // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    int iy = 1;                      // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    int iz = 0;                      // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    double igridcell_volume = 1;     // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    double iwater_volume = 0.5;      // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    OpenWQ_initiate.readSetIC(
        OpenWQ_json,
        OpenWQ_vars,
        OpenWQ_hostModelconfig,
        icmp,
        ix,
        iy,
        iz,
        igridcell_volume,   // asssuming units = m3
        iwater_volume);     // asssuming units = m3

    // OpenWQ_watertransp
    OpenWQ_watertransp OpenWQ_watertransp;   // create object: transport modules
    OpenWQ_chem OpenWQ_chem;                 // create object: biochemistry modules
    //OpenWQ_print OpenWQ_print;             // print modules
    
    int ts_hosthydromod = 1000; // (timesteps) TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
    
    // Loop: Time 
    // Space loop is inside OpenWQ_chem.Run function
    // No space loop for OpenWQ_watertransp.Adv: it needs to be called throughout the host model code
    for (int ts=0;ts<ts_hosthydromod;ts++){ // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL

        int source = 1;                 // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int ix_s = 1;                   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int iy_s = 1;                   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int iz_s = 0;                   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int recipient = 2;              // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int ix_r = 1;                   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int iy_r = 1;                   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        int iz_r = 0;                   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        double wflux_s2r = 0.01;        // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL
        double wmass_recipient = 0.1;   // TO REMOVE/REPLACE IN HOST HYDROLOGICAL MODEL


        /* ###################################################################################
        Transport with water fluxes
        ################################################################################### */ 

        // FUNCTION 1: Just Advection (Land Surface model)
        OpenWQ_watertransp.Adv(
            OpenWQ_json,
            OpenWQ_vars,
            source,
            ix_s, 
            iy_s,
            iz_s,
            recipient,
            ix_r,
            iy_r,
            iz_r,
            wflux_s2r,
            wmass_recipient
            );

        // FUNCTION 2: Advection and Dispersion (Aquatic Systems)
        //OpenWQ_watertransp.AdvDisp(
        //    OpenWQ_json,
        //    OpenWQ_vars);

        //
        //void OpenWQ_watertransp::ChemCompExchange(
        //    OpenWQ_json& OpenWQ_json, 
        //    OpenWQ_vars& OpenWQ_vars, 
        //    int source, std::string kinetics, 
        //    std::vector<std::string> parameter_names, 
        //    std::vector<double> parameter_values,
        //    std::array<double,7> & linedata, 
        //    int & index_chem){


        /* ###################################################################################
        Biogeochemistry
        ################################################################################### */ 

        OpenWQ_chem.Run(
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_hostModelconfig);

        /* ###################################################################################
        Sources and Sink
        ################################################################################### */ 

        // Need to create functions for this

        /* ###################################################################################
        Print Results
        ################################################################################### */ 
        //OpenWQ_print.writeVTU(OpenWQ_json,OpenWQ_vars,num_HydroComp,0); 

    }
    

}

