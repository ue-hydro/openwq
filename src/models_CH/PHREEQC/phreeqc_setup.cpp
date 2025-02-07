
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

#include "models_CH/headerfile_CH.hpp"


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

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm = std::unique_ptr<PhreeqcRM>(new PhreeqcRM(nxyz, OpenWQ_wqconfig.get_num_threads_requested()));


    try {
        component_h2o = OpenWQ_json.BGC_module["COMPONENT_H2O"];
    } catch (...){
        component_h2o = true;
    }

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetComponentH2O(component_h2o);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetFilePrefix("phreeqc");
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->OpenFiles();


    database = OpenWQ_json.BGC_module["DATABASE"];
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->LoadDatabase(database);
    PHREEQC_filename = OpenWQ_json.BGC_module["FILEPATH"];
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunFile(true,true, true, PHREEQC_filename);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->FindComponents();

}
