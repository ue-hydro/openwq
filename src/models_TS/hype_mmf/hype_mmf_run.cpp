
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

#include "models_TS/headerfile_TS.hpp"
#include<cmath>

// Calculate eroded particles from soil with a method based on the Morgan-Morgan-Finney erosion model. (from HYPE)
//
// Reference: R. P. C. Morgan, D. D. V. Morgan, and H. J. Finney, 1984. 
// A predictive model for the assessment of soil erosion risk, Journal 
// of Agricultural Engineering Research, vol. 30, pp. 245�253
// 

/*
TODO: 
1) make if statement where this function only runs if flow parallel to soil (see LE appraoch)
3) update snow status for every timestep (if statement for if upper comparment = icmp_inhibit_erosion): 
		if input_compartment_inhibitErosion_index = -1 (it mean no inhibiting compartment)
		// Get compartment that inhibits erosion
		icmp_inhibit_erosion = std::get<3>(OpenWQ_wqconfig.TS_model->HypeMMF->info_vector);

		// Reset to zero
		OpenWQ_wqconfig.TS_model->HypeMMF->snow_entryArmaCube(ix, iy, iz).zeros();

		// Loop over SUMMA results to save snow status
		OpenWQ_wqconfig.TS_model->HypeMMF->snow_entryArmaCube(ix, iy, iz) = snow data from summa

4) inhibit erosion if snow in icmp_inhibit_erosion > 0
*/


void OpenWQ_TS_model::mmf_hype_erosion_run(
	OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
	OpenWQ_vars& OpenWQ_vars,
	OpenWQ_wqconfig& OpenWQ_wqconfig,
	OpenWQ_utils& OpenWQ_utils,
	OpenWQ_output& OpenWQ_output,
	const int source, const int ix_s, const int iy_s, const int iz_s,
	const int recipient, const int ix_r, const int iy_r, const int iz_r,
	double wflux_s2r, 
	std::string TS_type){
	

	// ###################
	// CHECK IF RUN APPLICABLE

	// Return if no flux: wflux_s2r == 0
	if(wflux_s2r == 0.0f){return;}


	// Get TS model elements
	int exchange_direction = OpenWQ_wqconfig.TS_model->HypeMMF->get_exchange_direction();
  int erodFlux_icmp = OpenWQ_wqconfig.TS_model->HypeMMF->get_eroding_transp_compartment();
  int erodInhibut_icmp = OpenWQ_wqconfig.TS_model->HypeMMF->get_erosion_inhibit_compartment();
	int sediment_icmp = OpenWQ_wqconfig.TS_model->HypeMMF->get_sediment_compartment();

	// Return if there is snow
	double snow = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
				erodInhibut_icmp,
				ix_s,iy_s,iz_s);
	if (snow != 0.0)
		return;

	// Compute pseudo day-of-year from simulation time
	// Used for seasonal variation in rainfall energy calculation
	{
		time_t simtime = OpenWQ_wqconfig.get_simtime();
		struct tm tm_sim;
		gmtime_r(&simtime, &tm_sim);
		OpenWQ_wqconfig.TS_model->HypeMMF->pdayno = tm_sim.tm_yday + 1; // tm_yday is 0-based, pdayno is 1-based
	}

	// ######################
	// If runoff mobilization: TS_type_LE
	if (TS_type.compare("TS_type_LE") == 0){
		
		// Return if flux across compartments
		// This lateral mixing only occurs when fluxes occur when fluxes along the interface of compartments or if leaving the system
		// Only applicable for runoff erosion but not to precipitation monilization
		if (source != recipient && recipient != -1)
			return;

		// Return if upper compartment not that defined in HypeMMF
		if (source != erodFlux_icmp)
			return;

	}
	
	// ################################
	// Eroding / Transport compartment (upper compartment)
	
	// info for eroding/transport compartment 
	std::vector<int> n_xyz_erodingTransp_compartment;

	n_xyz_erodingTransp_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(erodFlux_icmp));
	n_xyz_erodingTransp_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(erodFlux_icmp));
	n_xyz_erodingTransp_compartment.push_back(OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(erodFlux_icmp));
	
	unsigned int index_erodTranp_interface = n_xyz_erodingTransp_compartment[exchange_direction] - 1;

	// Get number of cells in input_direction_index
	std::vector<unsigned int> xyz_source;
	xyz_source.push_back(ix_s);
	xyz_source.push_back(iy_s);
	xyz_source.push_back(iz_s);         // get x,y,z from source comparment in vector

	// Ignore if current source cell is not the lower adjacent cell where mixing may occur
	if (xyz_source[exchange_direction] != index_erodTranp_interface)
		return;

	// ################################
	// Sediment compartment (lower compartment)

	// Now get the coordenates of the lower compartment from which the mass exchange will take place
	// which will be the boundary cell: x or y or z = 0
	std::vector<unsigned int> xyz_SedCmpt_interface;
	xyz_SedCmpt_interface = xyz_source;
	xyz_SedCmpt_interface[exchange_direction] = 0; 

	// ##################
	// CODE

	// Local variables
	
	// PARAMETERS IN
	double cohesion = OpenWQ_wqconfig.TS_model->HypeMMF->cohesion_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double erodibility = OpenWQ_wqconfig.TS_model->HypeMMF->erodibility_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double sreroexp = OpenWQ_wqconfig.TS_model->HypeMMF->sreroexp_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double cropcover = OpenWQ_wqconfig.TS_model->HypeMMF->cropcover_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double groundcover = OpenWQ_wqconfig.TS_model->HypeMMF->groundcover_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double slope = OpenWQ_wqconfig.TS_model->HypeMMF->slope_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double trans1 = OpenWQ_wqconfig.TS_model->HypeMMF->trans_1_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
	double trans2 = OpenWQ_wqconfig.TS_model->HypeMMF->trans_2_entryArmaCube(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);
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
			COMPARTMENT: " + OpenWQ_hostModelconfig.get_HydroComp_name_at(sediment_icmp) 
			+ "(" + std::to_string(xyz_SedCmpt_interface[0]) + " ," 
			+ std::to_string(xyz_SedCmpt_interface[1]) + " ," 
			+ std::to_string(xyz_SedCmpt_interface[2]) + ")";

		// Print it (Console and/or Log file)
		OpenWQ_output.ConsoleLog(
			OpenWQ_wqconfig,    // for Log file name
			msg_string,         // message
			true,               // print in console
			true);              // print in log file

		// Skip this entry as cover fractions are invalid
		return;

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
			->mobilisedsed_rain_potential(ix_r, iy_r, iz_r) 
			= mobilisedsed_EWF;

	
	} else if (TS_type.compare("TS_type_LE") == 0){
	// ###############################################
	// If erosion driven by LE_model (e.g., runoff)

	// 1)
	// Mobilization from SedCompt (sediment compartment)
	
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
	mobilisedsed_EWF = OpenWQ_wqconfig.TS_model->HypeMMF->mobilisedsed_rain_potential(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);

	(*OpenWQ_vars.d_sedmass_mobilized_dt)(xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]) 
		= 1000.0 * (
			mobilisedsed_EWF    // mobilization by rainfall
			+ mobilisedsed_LE   // mobilization by runoff/flow
			) * transportfactor;  // kg/km2

	}

	// 2)
	// Mobilization with flow
	// Assuming here that, since d_sedmass_mobilized_dt is the mobilized sediment,
	// then it will all move with flow regardless of the flow/runoff intensity
	if(ix_r!=-1 && iy_r!=-1 && iz_r!=-1){

		(*OpenWQ_vars.d_sedmass_dt)(ix_r, iy_r, iz_r) 
			+= (*OpenWQ_vars.d_sedmass_mobilized_dt)(
					xyz_SedCmpt_interface[0], xyz_SedCmpt_interface[1], xyz_SedCmpt_interface[2]);

	// Move sediments with flow: from source to recipient
	// removing from source
	(*OpenWQ_vars.d_sedmass_dt)(ix_s, iy_s, iz_s) -= (*OpenWQ_vars.sedmass)(ix_s, iy_s, iz_s);

	// adding to recipient
	(*OpenWQ_vars.d_sedmass_dt)(ix_r, iy_r, iz_r) += (*OpenWQ_vars.sedmass)(ix_s, iy_s, iz_s);
	
	}

}

/*
	USE MODVAR, ONLY : basin, &
											pi                              

	!Argument declarations
	INTEGER, INTENT(IN) :: i             !<index of current subbasin
	INTEGER, INTENT(IN) :: j             !<index of current class
	INTEGER, INTENT(IN) :: pdayno        !<current pseudo day number of the year
	REAL, INTENT(IN)    :: prec          !<precipitation (rainfall only)    
	REAL, INTENT(IN)    :: surfacerunoff !<saturated overland flow and excess infiltration (mm)
	REAL, INTENT(IN)    :: cohesion      !<(kPa)
	REAL, INTENT(IN)    :: erodibility   !<(g/J)
	REAL, INTENT(IN)    :: snow          !<snow water (mm)
	REAL, INTENT(IN)    :: sreroexp      !<surface runoff erosion exponent 
	REAL, INTENT(IN)    :: flow          !<fast flow
	REAL, INTENT(OUT)   :: erodedsed     !<eroded (transported) sediment (kg/km2)

	!Local variables
	REAL Rainfall_energy
	REAL intensity
	REAL common_cropcover, common_groundcover
	REAL transportfactor
	REAL mobilisedsed         !mobilised suspended sediment from rainfall and surface runoff (g/m2)
	
	!Local parameters
	REAL, PARAMETER :: trans1 = 4.0
	REAL, PARAMETER :: trans2 = 1.3  

	!>\b Algorithm
	erodedsed = 0.
	IF(cohesion==0.OR.erodibility==0) RETURN      !no parameter values -> no erosion

	!>Calculate current cropcover and groundcover, will limit erosion
	CALL get_current_crop_and_ground_cover(i,j,common_cropcover,common_groundcover)

	!Check for snow limiting erosion
	intensity = 1. !intenspar
	IF(snow>0.) intensity = 0.    !snow

	!>Calculate particles that is eroded by rain splash detachment and by overland flow (mobilised sediment)
	mobilisedsed = 0.
	IF(prec > 0.) THEN
		IF(intensity > 0.) THEN
			IF(prec>5.0) THEN     !TODO: shorter timestep, other threshold?, holds for all over the world?, reference?
				Rainfall_energy = 8.95+8.44*LOG10(prec*(0.257+sin(2*3.14*((pdayno-70.)/365.))*0.09)*2.)
			ELSE
				Rainfall_energy = 0.
			ENDIF
			Rainfall_energy = prec * Rainfall_energy        !J/m2
			mobilisedsed = Rainfall_energy * (1. - common_cropcover) * erodibility  !g/m2
		ENDIF
	ENDIF
	IF(surfacerunoff > 0.) THEN   
		mobilisedsed = mobilisedsed + (((surfacerunoff * 365.) ** sreroexp) * (1. - common_groundcover) * (1./(0.5 * cohesion)) * SIN(basin(i)%slope / 100.)) / 365. !g/m2   
	ENDIF

	!>Transport capacity of fast flowing water may limit transport of sediment
	IF(flow>0.)THEN
		transportfactor = MIN(1.,(flow / trans1)**trans2)   
	ELSE
		transportfactor = 1.
	ENDIF
		
	!>Eroded sediment calculated from mobilised sediment, possibly limited by the transport capacity
	erodedsed = 1000. * mobilisedsed * transportfactor  !kg/km2

END SUBROUTINE calculate_MMF_erosion
*/