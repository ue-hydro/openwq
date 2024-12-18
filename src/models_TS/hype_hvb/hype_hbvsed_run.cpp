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

#include "models_TS/headerfile_TS.hpp"

/* ########################################
// Calculate eroded particles from soil by HBV-sed based model
#########################################*/
void OpenWQ_TS_model::hbvsed_hype_erosion_run(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient,  // no need for ix_r, iy_r, iz_r (because mobilization occurs at the adjacent cell of the upper compartment)
    double wflux_s2r, 
    double wmass_source,
    std::string TS_type){

  /*
	// ###################
	// CHECK IF RUN APPLICABLE

	// Return if no flux: wflux_s2r == 0
	if(wflux_s2r == 0.0f){return;}

	// Return if flux across compartments
	// This lateral mixing only occurs when fluxes occur when fluxes along the interface of compartments or if leaving the system
	if (source != recipient && recipient != -1)
		return;
	
	// Return if upper compartment not that defined in HypeMMF
	int erodFlux_icmp = OpenWQ_wqconfig.TS_model->HypeMMF->get_eroding_flux_compartment();
	if (source != erodFlux_icmp)
		return;

	//OpenWQ_hostModelconfig.get_HydroComp_index(erodFlux_icmp)

	// Local variables
	int exchange_direction = OpenWQ_wqconfig.TS_model->HypeMMF->get_exchange_direction();
	std::vector<int> xyz_erodingFlux_compartment;
	std::vector<unsigned int> xyz_lowerComp;

	// Num of cells in upper_compartment
	xyz_erodingFlux_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(erodFlux_icmp));
	xyz_erodingFlux_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(erodFlux_icmp));
	xyz_erodingFlux_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(erodFlux_icmp));
	
	// Get number of cells in input_direction_index
	unsigned int index_lower_cell = xyz_erodingFlux_compartment[exchange_direction] - 1;
	std::vector<unsigned int> xyz_source;
	xyz_source.push_back(ix_s);
	xyz_source.push_back(iy_s);
	xyz_source.push_back(iz_s);         // get x,y,z from source comparment in vector

	// Ignore if current source cell is not the lower adjacent cell where mixing may occur
	if (xyz_source[exchange_direction] != index_lower_cell)
		return;

	// Return if there is snow
  // TODO: if (frostdepth<0. .AND. frostdepth>-9999.) then return
	double snow = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
				OpenWQ_wqconfig.TS_model->HypeMMF->get_erosion_inhibit_compartment(),
				ix_s,iy_s,iz_s);
	if (snow != 0.0)
		return;

	// Now get the coordenates of the lower compartment from which the mass exchange will take place
	// which will be the boundary cell: x or y or z = 0
	xyz_lowerComp = xyz_source;
	xyz_lowerComp[exchange_direction] = 0; 

	// ##################
	// CODE

	// Local variables
	
	// PARAMETERS IN
	double cohesion = OpenWQ_wqconfig.TS_model->HypeMMF->cohesion_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double erodibility = OpenWQ_wqconfig.TS_model->HypeMMF->erodibility_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double sreroexp = OpenWQ_wqconfig.TS_model->HypeMMF->sreroexp_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double cropcover = OpenWQ_wqconfig.TS_model->HypeMMF->cropcover_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double groundcover = OpenWQ_wqconfig.TS_model->HypeMMF->groundcover_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double slope = OpenWQ_wqconfig.TS_model->HypeMMF->slope_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double trans1 = OpenWQ_wqconfig.TS_model->HypeMMF->trans_1_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	double trans2 = OpenWQ_wqconfig.TS_model->HypeMMF->trans_2_entryArmaCube(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);
	// DUMMY
	double prec_erosion_threshold;
	double rainfall_energy;
	double transportfactor;
	// OUT
	double mobilisedsed_EWF = 0.0; // mobilised suspended sediment from rainfall (g/m2)
	double mobilisedsed_LE = 0.0;  // mobilised suspended sediment from surface runoff (g/m2)
	double erodedsed = 0.0;        // <eroded (transported) sediment (kg/km2)
	
	// Check if cropcover + groundcover is equal to 1
	// Otherwise, scale it back
	std::string msg_string;
	if (cropcover + groundcover != 1.0){

		// Create Message (Warning Message)
		msg_string = 
			"<OpenWQ> WARNING: TS_model -> MMF_HYPE - \"CROPCOVER + GROUNDCOVER\" needs to be equal to 1 (entry skiped): \\
			COMPARTMENT: " + OpenWQ_wqconfig.TS_model->SedCmpt 
			+ "(" + std::to_string(xyz_lowerComp[0]) + " ," 
			+ std::to_string(xyz_lowerComp[1]) + " ," 
			+ std::to_string(xyz_lowerComp[2]) + ")";

		// Print it (Console and/or Log file)
		OpenWQ_output.ConsoleLog(
			OpenWQ_wqconfig,    // for Log file name
			msg_string,         // message
			true,               // print in console
			true);              // print in log file

	}

	// #################################
	// Get relevant parameter vals at ix, iy, iz
	// #################################

	// There is no erosion if any of the following conditions are met
	if(cohesion == 0.0 || erodibility == 0.0){
		return;
	};       

	// ###############################################
	// Calculate particles that are eroded by rain splash detachment and by overland flow (mobilised sediment)
	// ###############################################

	// ###############################################
	// If erosion driven by EWF (e.g., precipitation)
	if (TS_type.compare("TS_type_EWF") == 0){
		
		// Hype uses the threshold of 5 mm per day for erosion to occur
		// Here, because OpenWQ can have different timesteps, this threshold is adjusted (simple rule of three)
		prec_erosion_threshold = 5 / (24 * 60 * 60) * OpenWQ_hostModelconfig.get_time_step();
		
		if(wflux_s2r > prec_erosion_threshold){     // shorter timestep, other threshold?, holds for all over the world?, reference?
			
			rainfall_energy = 8.95 + 8.44 * log10(
				wflux_s2r * (0.257 + sin(
					2 * 3.14 * ((OpenWQ_wqconfig.TS_model->HypeMMF->pdayno - 70.0)/365)) * 0.09) * 2.0);

		}else{

			rainfall_energy = 0.0;

		}

		rainfall_energy = wflux_s2r * rainfall_energy;        // J/m2
		mobilisedsed_EWF = rainfall_energy * (1.0 - cropcover) * erodibility;  //g/m2

		// Adding rain mobilization to potential
		OpenWQ_wqconfig.TS_model->HypeMMF
			->mobilisedsed_rain_potential(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]) 
			= mobilisedsed_EWF;

	
	} else if (TS_type.compare("TS_type_LE") == 0){
	// ###############################################
	// If erosion driven by LE_model (e.g., runoff)
	
		// TODO 
		// wflux_s2r needs to be in mm
		// (((surfacerunoff * 365.) ** sreroexp) * (1. - common_groundcover) * (1./(0.5 * cohesion)) * SIN(basin(i)%slope / 100.)) / 365. !g/m2 
		mobilisedsed_LE = ( pow(wflux_s2r * 365, sreroexp) 
			* (1.0 - groundcover) * (1.0/(0.5f * cohesion)) * 
			sin(slope / 100.0)) / 365; // g/m2   

		mobilisedsed_LE = std::max(mobilisedsed_LE, 0.0);

		// Transport capacity of fast flowing water may limit transport of sediment
		// transportfactor = MIN(1.,(flow / trans1)**trans2) 
		transportfactor = std::min(1.0d, pow(wflux_s2r / trans1, trans2));

		// #####################################
		// Eroded sediment calculated from mobilised sediment, possibly limited by the transport capacity
		
		// Mobilization by rainfall if transport/flow exists
		mobilisedsed_EWF = OpenWQ_wqconfig.TS_model->HypeMMF->mobilisedsed_rain_potential(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]);

		(*OpenWQ_vars.d_sedmass_dt)(xyz_lowerComp[0], xyz_lowerComp[1], xyz_lowerComp[2]) 
			= 1000.0 * (
				mobilisedsed_EWF    // mobilization by rainfall
				+ mobilisedsed_LE   // mobilization by runoff/flow
				) * transportfactor;  // kg/km2

	}

  */
   
 
}

/*
void hbvsed_hype_erosion_run(i,prec,lusepar,soilpar,slopepar, &
                    precexppar,eroindexpar,snow,frostdepth,erodedsed){

    USE HYPEVARIABLES, ONLY : m_erodmon
    USE MODVAR, ONLY : basin,monthpar,current_time

    // Argument declarations
    int i             // index of current subbasin
    double prec          // precipitation (rainfall only)    
    double lusepar       // soil erosion factor (land use dependence)
    double soilpar       // soil erosion factor (soil dependence)
    double slopepar      // slope erosion factor (exponent) 
    double precexppar    // erosion precipitation dependence factor (exponent)
    double eroindexpar   // model parameter for scaling of erosion index
    double snow          // snow water (mm)
    double frostdepth    // frost depth (cm)
    double erodedsed     // eroded (transported) sediment (kg/km2)

    // Local variables
    double erosionindex;
    double mobilisedsed;   //mobilised suspended sediment from rainfall (g/m2)
    double erodmonth;      // monthly erosion factor (-)
    
    // Algorithm
    // now and soil frost turn off particle mobilization
    IF(snow>0. .OR. (frostdepth<0. .AND. frostdepth>-9999.))THEN
      erodedsed = 0.
      RETURN
    ENDIF

    // Calculate mobilised sediments from erosion index, slope and rainfall
    erosionindex = basin(i)%eroindex / eroindexpar    
    mobilisedsed = (basin(i)%slope / 5.)**slopepar * lusepar * soilpar * erosionindex * prec**precexppar       !tonnes/km2 = g/m2
      
    // Eroded sediment calculated from mobilised sediment with a monthly factor
    erodmonth = 1. + monthpar(m_erodmon,current_time%date%month)
    erodedsed = 1000. * mobilisedsed * erodmonth  // kg/km2

}
*/