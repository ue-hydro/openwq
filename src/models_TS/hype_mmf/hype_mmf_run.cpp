
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
#include<cmath>

// Calculate eroded particles from soil with a method based on the Morgan-Morgan-Finney erosion model. (from HYPE)
//
// Reference: R. P. C. Morgan, D. D. V. Morgan, and H. J. Finney, 1984. 
// A predictive model for the assessment of soil erosion risk, Journal 
// of Agricultural Engineering Research, vol. 30, pp. 245�253
// 

void OpenWQ_TS_model::mmf_hype_erosion_run(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    double wmass_source,
    std::string TS_type){
    
    // Local variables
    // IN  
    double cohesion = OpenWQ_wqconfig.TS_model->HypeMMF->cohesion_entryArmaCube(ix_r, iy_r, iz_r);
    double erodibility = OpenWQ_wqconfig.TS_model->HypeMMF->erodibility_entryArmaCube(ix_r, iy_r, iz_r);
    double sreroexp = OpenWQ_wqconfig.TS_model->HypeMMF->sreroexp_entryArmaCube(ix_r, iy_r, iz_r);
    double cropcover = OpenWQ_wqconfig.TS_model->HypeMMF->cropcover_entryArmaCube(ix_r, iy_r, iz_r);
    double groundcover = OpenWQ_wqconfig.TS_model->HypeMMF->groundcover_entryArmaCube(ix_r, iy_r, iz_r);
    double slope = OpenWQ_wqconfig.TS_model->HypeMMF->slope_entryArmaCube(ix_r, iy_r, iz_r);
    double snow = OpenWQ_wqconfig.TS_model->HypeMMF->snow_entryArmaCube(ix_r, iy_r, iz_r);
    double trans1 = OpenWQ_wqconfig.TS_model->HypeMMF->trans_1_entryArmaCube(ix_r, iy_r, iz_r);
    double trans2 = OpenWQ_wqconfig.TS_model->HypeMMF->trans_2_entryArmaCube(ix_r, iy_r, iz_r);

    // DUMMY
    double prec_erosion_threshold;
    double rainfall_energy;
    double intensity;
    double transportfactor;

    // OUT
    double mobilisedsed_EWF = 0.0; // mobilised suspended sediment from rainfall (g/m2)
    double mobilisedsed_LE = 0.0;  // mobilised suspended sediment from surface runoff (g/m2)
    double erodedsed;               // <eroded (transported) sediment (kg/km2)
    
    // #################################
    // Get relevant parameter vals at ix, iy, iz
    // #################################

    // if cohesion or erodibility is zero, then no erosion
    erodedsed = 0.0;
    if(cohesion == 0.0 || erodibility == 0.0 ){
      return;
    };       

    // Check for snow limiting erosion
    intensity = 1.0; 
    if(snow > 0.0){ 
      intensity = 0.0;
    };

    // Calculate particles that are eroded by rain splash detachment and by overland flow (mobilised sediment)
    
    // ###############################################
    // If erosion driven by EWF (e.g., precipitation)
    // ###############################################
    if (TS_type.compare("TS_type_EWF") && 
        wflux_s2r > 0.0){
      
        if(intensity > 0.0){

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

          OpenWQ_wqconfig.TS_model->HypeMMF->mobilisedsed_rain_potential = mobilisedsed_EWF;

        // ###############################################
        // If erosion driven by LE_model (e.g., runoff)
        // ###############################################
        } else if (TS_type.compare("TS_type_LE")
          && wflux_s2r > 0.0){

          // TODO 
          // wflux_s2r needs to be in mm
          // (((surfacerunoff * 365.) ** sreroexp) * (1. - common_groundcover) * (1./(0.5 * cohesion)) * SIN(basin(i)%slope / 100.)) / 365. !g/m2 
          mobilisedsed_LE = ( pow(wflux_s2r * 365, sreroexp) 
            * (1.0 - groundcover) * (1.0/(0.5f * cohesion)) * 
            sin(slope / 100.0)) / 365; // g/m2   


          // Transport capacity of fast flowing water may limit transport of sediment
          if(wflux_s2r  > 0.0){
            // transportfactor = MIN(1.,(flow / trans1)**trans2) 
            transportfactor = std::min(1.0d, pow(wflux_s2r / trans1, trans2));
          }else{
            transportfactor = 1.0;
          }

          // #####################################
          // Eroded sediment calculated from mobilised sediment, possibly limited by the transport capacity
          
          // Mobilization by rainfall if transport/flow exists
          mobilisedsed_EWF = OpenWQ_wqconfig.TS_model->HypeMMF->mobilisedsed_rain_potential(ix_r,iy_r,iz_r);

          (*OpenWQ_vars.d_sedmass_dt)(ix_r,iy_r,iz_r) 
            = 1000.0 * (
              mobilisedsed_EWF    // mobilization by rainfall
              + mobilisedsed_LE   // mobilization by runoff/flow
              ) * transportfactor;  // kg/km2

        }

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