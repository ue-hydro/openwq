#include "models_CH/headerfile_CH.hpp"
#include <cmath>  // for isnan, isinf

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

    // Get the actual component count from PhreeqcRM
    int nc = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetComponentCount();
    msg_string = "<OpenWQ> PHREEQC phreeqc_run: nxyz=" + std::to_string(nxyz)
        + ", nc=" + std::to_string(nc)
        + ", num_chem=" + std::to_string(OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Use nc (actual PhreeqcRM component count) for sizing
    c.resize(nxyz * nc);
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
                        double temp_val = OpenWQ_hostModelconfig.get_dependVar_at(temperature_indx, ix, iy, iz);
                        // Auto-detect if temperature is in Kelvin (>100 suggests Kelvin)
                        // PHREEQC expects Celsius, so convert if needed
                        // Typical water temps: 0-40°C or 273-313K
                        if (temp_val > 100.0) {
                            // Likely Kelvin - convert to Celsius
                            temp_val = temp_val - 273.15;
                        }
                        temperatures[indx] = temp_val;
                    }
                    if (pressure_indx != -1) {
                        pressures[indx] = OpenWQ_hostModelconfig.get_dependVar_at(pressure_indx, ix, iy, iz);
                    }
                    indx++;

                }
            }
        }
    }


    // Initialize concentration array - use nc from PhreeqcRM
    // Note: num_chem from config might differ from actual PhreeqcRM component count
    unsigned int num_chem_to_use = std::min((unsigned int)nc, OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);
    for (unsigned int chemi=0;chemi<num_chem_to_use;chemi++){
        indx = 0;
        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){


            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        double val = (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz)/volumes[indx];
                        // Validate - replace NaN/Inf with 0
                        if (std::isnan(val) || std::isinf(val) || val < 0) val = 0.0;
                        c[chemi*nxyz+indx] = val;
                        indx++;
                    }
                }
            }
        }
    }

    indx = 0;

    // Validate and fix temperatures - PHREEQC needs reasonable values (0-100°C typical)
    for (int i = 0; i < nxyz; i++) {
        if (std::isnan(temperatures[i]) || std::isinf(temperatures[i]) || temperatures[i] < -50 || temperatures[i] > 200) {
            temperatures[i] = 25.0;  // Default to 25°C
        }
        if (std::isnan(pressures[i]) || std::isinf(pressures[i]) || pressures[i] <= 0) {
            pressures[i] = 1.0;  // Default to 1 atm
        }
    }

    msg_string = "<OpenWQ> PHREEQC: Setting T[0]=" + std::to_string(temperatures[0])
        + ", P[0]=" + std::to_string(pressures[0])
        + ", c[0]=" + std::to_string(c[0])
        + ", vol[0]=" + std::to_string(volumes[0]);
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetPressure(pressures);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTemperature(temperatures);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetConcentrations(c);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetRepresentativeVolume(volumes);

    // Set fixed density and saturation to avoid calc_dens() crash
    std::vector<double> density(nxyz, 1.0);  // 1.0 kg/L
    std::vector<double> saturation(nxyz, 1.0);  // Fully saturated
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetDensityUser(density);
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetSaturationUser(saturation);

    // CRITICAL: Tell PhreeqcRM to use user-provided density, not calculate internally
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->UseSolutionDensityVolume(false);

    msg_string = "<OpenWQ> PHREEQC: About to call RunCells with timestep="
        + std::to_string(OpenWQ_hostModelconfig.get_time_step());
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTimeStep(OpenWQ_hostModelconfig.get_time_step());
    if (OpenWQ_hostModelconfig.get_time_step() != 0)    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunCells();
    msg_string = "<OpenWQ> PHREEQC: RunCells completed, getting concentrations";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetConcentrations(c);
    for (unsigned int chemi=0;chemi<num_chem_to_use;chemi++){
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
