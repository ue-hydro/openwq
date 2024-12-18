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
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    std::string TS_type){


	// ###################
	// CHECK IF RUN APPLICABLE

  // only applicable with precipitation input
  if (TS_type.compare("TS_type_LE") != 0)
    return;

	// Return if no flux: wflux_s2r == 0
	if(wflux_s2r == 0.0f){return;}

	// Local variables
	int exchange_direction = OpenWQ_wqconfig.TS_model->HypeHVB->get_exchange_direction();
	std::vector<unsigned int> xyz_lowerComp;

	// Get number of cells in input_direction_index
	std::vector<unsigned int> xyz_SedCompt;
	xyz_SedCompt.push_back(ix_r);
	xyz_SedCompt.push_back(iy_r);
	xyz_SedCompt.push_back(iz_r);         // get x,y,z from source comparment in vector

	// Ignore if current sediment source cell is not the top one
	if (xyz_SedCompt[exchange_direction] != 0)
		return;

	// Return if there is snow
  // TODO: if (frostdepth<0. .AND. frostdepth>-9999.) then return
	double snow = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
				OpenWQ_wqconfig.TS_model->HypeHVB->get_erosion_inhibit_compartment(),
				ix_s,iy_s,iz_s);
	if (snow != 0.0)
		return;

	// ##################
	// CODE

	// Local variables
	
	// PARAMETERS IN
	double lusepar = OpenWQ_wqconfig.TS_model->HypeHVB->lusepar_entryArmaCube(ix_r, iy_r, iz_r);
	double soilpar = OpenWQ_wqconfig.TS_model->HypeHVB->soilpar_entryArmaCube(ix_r, iy_r, iz_r);
  double slopepar = OpenWQ_wqconfig.TS_model->HypeHVB->slopepar_entryArmaCube(ix_r, iy_r, iz_r);
  double precexppar = OpenWQ_wqconfig.TS_model->HypeHVB->precexppar_entryArmaCube(ix_r, iy_r, iz_r);
  double eroindexpar = OpenWQ_wqconfig.TS_model->HypeHVB->eroindexpar_entryArmaCube(ix_r, iy_r, iz_r);
  double slope = OpenWQ_wqconfig.TS_model->HypeHVB->slope_entryArmaCube(ix_r, iy_r, iz_r);
  double eroindex = OpenWQ_wqconfig.TS_model->HypeHVB->eroindex_entryArmaCube(ix_r, iy_r, iz_r);

	// DUMMY
	double erosionindex;
  double erodmonth;     // erodmonth
	double mobilisedsed;  // mobilised suspended sediment from rainfall (g/m2)

	// ###############################################
	// Calculate particles that are eroded by rain splash detachment and by overland flow (mobilised sediment)
	// ###############################################

  // Calculate mobilised sediments from erosion index, slope and rainfall
  erosionindex = eroindex / eroindexpar;
  mobilisedsed = pow((slope / 5.),slopepar)
                 * lusepar * soilpar * erosionindex 
                 * pow(wflux_s2r,precexppar); // !tonnes/km2 = g/m2
    
  // Eroded sediment calculated from mobilised sediment with a monthly factor
  // TODO (from hype): erodmonth = 1. + monthpar(m_erodmon,current_time%date%month)
  erodmonth = 1.0;
	(*OpenWQ_vars.d_sedmass_dt)(ix_r, iy_r, iz_r) = 1000. * mobilisedsed * erodmonth;  // kg/km2
	
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