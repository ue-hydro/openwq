#include "models_CH/headerfile_CH.hpp"

/* #################################################
// Compute each chemical transformation
################################################# */
void OpenWQ_CH_model::phreeqc_run(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output){

    std::vector<double> c, volumes, temperatures, pressures;
    int nxyz;
    int nx, ny, nz;
    int indx;
    std::string dependancy;
    std::string compName_icmp;
    unsigned int ndep;
    int temperature_indx;
    int pressure_indx;
    std::string msg_string;

    nxyz = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetGridCellCount();


    c.resize(nxyz * (OpenWQ_wqconfig.CH_model->PHREEQC->num_chem));
    volumes.resize(nxyz);
    temperatures.resize(nxyz);
    pressures.resize(nxyz);
    indx = 0;
    pressures = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetPressure();
    temperatures = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetTemperature();
    for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){
           

        nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
        ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
        nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
        compName_icmp = OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp);
        try {
            dependancy = OpenWQ_json.BGC_module["TEMPERATURE"][compName_icmp];
            temperature_indx = -1;
            for (ndep = 0; ndep < OpenWQ_hostModelconfig.get_num_HydroDepend(); ndep ++) {
                if (OpenWQ_hostModelconfig.get_HydroDepend_name_at(ndep) == dependancy) temperature_indx = ndep;
            }
            if (temperature_indx == -1) {
                // Create Message
                msg_string = "<OpenWQ> Unknown TEMPERATURE dependancy with name '"
                    + dependancy
                    + "' defined for compartment: " 
                    + compName_icmp;

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(
                    OpenWQ_wqconfig,    // for Log file name
                    msg_string,         // message
                    true,               // print in console
                    true);              // print in log file

            }
        } catch(...) {
            temperature_indx = -1;
        }

        try {
            dependancy = OpenWQ_json.BGC_module["PRESSURE"][compName_icmp];
            pressure_indx = -1;
            for (ndep = 0; ndep < OpenWQ_hostModelconfig.get_num_HydroDepend(); ndep ++) {
                if (OpenWQ_hostModelconfig.get_HydroDepend_name_at(ndep) == dependancy) pressure_indx = ndep;
            }
            if (pressure_indx == -1) {
                // Create Message
                msg_string = "<OpenWQ> Unknown PRESSURE dependancy with name '"
                    + dependancy
                    + "' defined for compartment: " 
                    + compName_icmp;

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(
                    OpenWQ_wqconfig,    // for Log file name
                    msg_string,         // message
                    true,               // print in console
                    true);              // print in log file

            }
        } catch(...) {
            pressure_indx = -1;
        }
        for (int ix=0; ix<nx; ix++) {
            for (int iy = 0; iy<ny;iy++) {
                for (int iz=0; iz<nz; iz++) {
                    volumes[indx] = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(icmp, ix, iy, iz);
                    if (volumes[indx] == 0) volumes[indx] = 1;

                    if (temperature_indx!=-1) {
                        temperatures[indx] = OpenWQ_hostModelconfig.get_dependVar_at(temperature_indx, ix, iy, iz);
                    }
                    if (pressure_indx != -1) {
                        pressures[indx] = OpenWQ_hostModelconfig.get_dependVar_at(pressure_indx, ix, iy, iz);
                    }
                    indx++;

                }
            }
        }
    }


    for (unsigned int chemi=0;chemi<(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);chemi++){
        indx = 0;
        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){
           

            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        c[chemi*nxyz+indx] = (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz)/volumes[indx];
                        indx++;
                    }
                }
            }
        }
    }

    indx = 0;

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetPressure(pressures);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTemperature(temperatures);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetConcentrations(c);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetRepresentativeVolume(volumes);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTimeStep(OpenWQ_hostModelconfig.get_time_step());
    if (OpenWQ_hostModelconfig.get_time_step() != 0)    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunCells();
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetConcentrations(c);
    for (unsigned int chemi=0;chemi<(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);chemi++){
        indx = 0;
        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){
           

            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(chemi)(ix,iy,iz) = c[chemi*nxyz+indx]*volumes[indx] - (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz);
                        indx++;
                    }
                }
            }
        }
    }


}
