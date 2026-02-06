#include "models_CH/headerfile_CH.hpp"
#include <cmath>  // for isnan, isinf

/* #################################################
// Compute each chemical transformation using PHREEQC
//
// UNIT CONVERSION NOTES:
// ----------------------
// OpenWQ internal units:
//   - chemass: stored in grams (g)
//   - concentration: mg/L (calculated as chemass/volume * 1000)
//   - volumes: liters (L)
//
// PHREEQC PhreeqcRM units (with SetUnitsSolution(2)):
//   - concentration: mol/kgw (moles per kilogram water)
//   - For dilute solutions: mol/kgw ≈ mol/L (since water density ~1 kg/L)
//
// Conversion formulas:
//   mg/L to mol/kgw:  conc_mol = conc_mg_L / (GFW * 1000)
//   mol/kgw to mg/L:  conc_mg_L = conc_mol * GFW * 1000
//
// where GFW = Gram Formula Weight in g/mol
//
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

    // Get GFW for unit conversion
    const std::vector<double>& gfw = OpenWQ_wqconfig.CH_model->PHREEQC->gfw;
    bool has_gfw = (gfw.size() >= (size_t)nc);

    if (!has_gfw) {
        msg_string = "<OpenWQ> WARNING: GFW not available for all components. "
            "Unit conversion disabled (assuming units already in mol/kgw).";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

    // Use nc (actual PhreeqcRM component count) for sizing
    c.resize(nxyz * nc);
    volumes.resize(nxyz);
    temperatures.resize(nxyz);
    pressures.resize(nxyz);
    indx = 0;
    pressures = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetPressure();
    temperatures = OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetTemperature();

    // First pass: collect volumes, temperatures, pressures
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


    // ########################################################
    // UNIT CONVERSION: OpenWQ (mg/L) -> PHREEQC (mol/kgw)
    // ########################################################
    // OpenWQ chemass is in grams (g)
    // concentration = chemass / volume gives g/L
    // To get mg/L: multiply by 1000
    // To get mol/kgw: divide mg/L by (GFW * 1000)
    // Combined: conc_mol = (chemass_g / volume_L) / GFW
    //         = chemass_g / (volume_L * GFW)
    // ########################################################

    unsigned int num_chem_to_use = std::min((unsigned int)nc, OpenWQ_wqconfig.CH_model->PHREEQC->num_chem);
    for (unsigned int chemi=0;chemi<num_chem_to_use;chemi++){
        indx = 0;

        // Get GFW for this component (default to 1.0 if not available)
        double component_gfw = (has_gfw && chemi < gfw.size()) ? gfw[chemi] : 1.0;

        // Protect against invalid GFW
        if (component_gfw <= 0 || std::isnan(component_gfw) || std::isinf(component_gfw)) {
            component_gfw = 1.0;
        }

        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){

            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        // OpenWQ: chemass is in grams, volume in liters
                        // conc_g_L = chemass_g / volume_L
                        double conc_g_L = (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz) / volumes[indx];

                        // Convert g/L to mol/kgw (mol/L for dilute solutions)
                        // mol/L = g/L / GFW(g/mol)
                        double conc_mol_kgw = conc_g_L / component_gfw;

                        // Validate - replace NaN/Inf with 0
                        if (std::isnan(conc_mol_kgw) || std::isinf(conc_mol_kgw) || conc_mol_kgw < 0) {
                            conc_mol_kgw = 0.0;
                        }

                        c[chemi*nxyz+indx] = conc_mol_kgw;
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

    // Log first cell values for debugging (with units)
    if (num_chem_to_use > 0 && has_gfw) {
        double first_gfw = (gfw.size() > 0) ? gfw[0] : 1.0;
        msg_string = "<OpenWQ> PHREEQC: Setting T[0]=" + std::to_string(temperatures[0]) + " C"
            + ", P[0]=" + std::to_string(pressures[0]) + " atm"
            + ", c[0]=" + std::to_string(c[0]) + " mol/kgw"
            + " (GFW=" + std::to_string(first_gfw) + " g/mol)"
            + ", vol[0]=" + std::to_string(volumes[0]) + " L";
    } else {
        msg_string = "<OpenWQ> PHREEQC: Setting T[0]=" + std::to_string(temperatures[0])
            + ", P[0]=" + std::to_string(pressures[0])
            + ", c[0]=" + std::to_string(c[0])
            + ", vol[0]=" + std::to_string(volumes[0]);
    }
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
        + std::to_string(OpenWQ_hostModelconfig.get_time_step()) + " s";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->SetTimeStep(OpenWQ_hostModelconfig.get_time_step());
    if (OpenWQ_hostModelconfig.get_time_step() != 0) {
        OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->RunCells();
    }
    msg_string = "<OpenWQ> PHREEQC: RunCells succeeded";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // Get updated concentrations from PHREEQC (in mol/kgw)
    OpenWQ_wqconfig.CH_model->PHREEQC->phreeqcrm->GetConcentrations(c);
    msg_string = "<OpenWQ> PHREEQC: GetConcentrations succeeded";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    // ########################################################
    // UNIT CONVERSION: PHREEQC (mol/kgw) -> OpenWQ (grams)
    // ########################################################
    // PHREEQC returns concentration in mol/kgw
    // Convert to g/L: conc_g_L = conc_mol_kgw * GFW
    // Convert to mass: mass_g = conc_g_L * volume_L
    // Delta mass = new_mass - old_mass
    // ########################################################

    for (unsigned int chemi=0;chemi<num_chem_to_use;chemi++){
        indx = 0;

        // Get GFW for this component
        double component_gfw = (has_gfw && chemi < gfw.size()) ? gfw[chemi] : 1.0;

        // Protect against invalid GFW
        if (component_gfw <= 0 || std::isnan(component_gfw) || std::isinf(component_gfw)) {
            component_gfw = 1.0;
        }

        for (unsigned int icmp=0;icmp<OpenWQ_hostModelconfig.get_num_HydroComp();icmp++){

            nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp); // num of x elements
            ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp); // num of y elements
            nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp); // num of z elements
            for (int ix=0; ix<nx; ix++) {
                for (int iy = 0; iy<ny;iy++) {
                    for (int iz=0; iz<nz; iz++) {
                        // PHREEQC concentration in mol/kgw
                        double conc_mol_kgw = c[chemi*nxyz+indx];

                        // Validate concentration
                        if (std::isnan(conc_mol_kgw) || std::isinf(conc_mol_kgw)) {
                            conc_mol_kgw = 0.0;
                        }
                        if (conc_mol_kgw < 0) {
                            conc_mol_kgw = 0.0;  // PHREEQC shouldn't return negative, but protect anyway
                        }

                        // Convert mol/kgw to g/L
                        // g/L = mol/kgw * GFW(g/mol)
                        double conc_g_L = conc_mol_kgw * component_gfw;

                        // Convert to mass: mass_g = conc_g_L * volume_L
                        double new_mass_g = conc_g_L * volumes[indx];

                        // Current mass in OpenWQ
                        double old_mass_g = (*OpenWQ_vars.chemass)(icmp)(chemi)(ix,iy,iz);

                        // Delta mass for the chemistry derivative
                        (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(chemi)(ix,iy,iz) = new_mass_g - old_mass_g;

                        indx++;
                    }
                }
            }
        }
    }

    // Log summary of first component change for debugging
    if (num_chem_to_use > 0 && nxyz > 0) {
        double first_gfw = (has_gfw && gfw.size() > 0) ? gfw[0] : 1.0;
        double conc_after_mol = c[0];
        double conc_after_g_L = conc_after_mol * first_gfw;
        msg_string = "<OpenWQ> PHREEQC: After reaction c[0]=" + std::to_string(conc_after_mol) + " mol/kgw"
            + " = " + std::to_string(conc_after_g_L) + " g/L"
            + " = " + std::to_string(conc_after_g_L * 1000) + " mg/L";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
    }

}
