
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


void OpenWQ_TS_model::mmf_hype_erosion(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, 
    double wmass_source){
    
    /*
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
    */

   /*
    // Local variables
    double rainfall_energy;
    double intensity;
    double cell_cropcover, cell_groundcover;
    double transportfactor;
    double mobilisedsed = 0.0f;         // mobilised suspended sediment from rainfall and surface runoff (g/m2)

    double erodedsed;     // <eroded (transported) sediment (kg/km2)

    erodedsed = 0.
    if(cohesion == 0.0 || erodibility == 0.0f ){
      return //no parameter values -> no erosion
    };       

    // Calculate current cropcover and groundcover, will limit erosion
    // XXXXXXXX = > get from input 
    cell_cropcover
    cell_groundcover
    cell_slope
    
    // Check for snow limiting erosion
    intensity = 1. // intenspar
    if(snow > 0.0f){ !snow
      intensity = 0.0f;
    };

    // Calculate particles that are eroded by rain splash detachment and by overland flow (mobilised sediment)
    if(prec > 0.0f){
      
      if(intensity > 0.0f){

        if(prec > 5.0){     // shorter timestep, other threshold?, holds for all over the world?, reference?
          
          rainfall_energy = 8.95 + 8.44 * log10(
            prec * (0.257 + sin(
               2 * 3.14 * ((pdayno - 70.0f)/365)) * 0.09) * 2.0f);

        }else{

          rainfall_energy = 0.0f;

        }

        rainfall_energy = prec * rainfall_energy;        // J/m2
        mobilisedsed = rainfall_energy * (1.0f - cell_cropcover) * erodibility;  //g/m2
      
      }
    }

    // if surfacerunoff exists
    if(surfacerunoff > 0.0f){   
      mobilisedsed = mobilisedsed + ((pow(surfacerunoff * 365),sreroexp) 
        * (1.0f - cell_groundcover) * (1.0f/(0.5f * cohesion)) * 
        sin(cell_slope / 100.0f)) / 365; // g/m2   
    }

    // Transport capacity of fast flowing water may limit transport of sediment
    if(flow  > 0.0f){
      transportfactor = min(1.0f, pow(flow / trans1), trans2));
    }else{
      transportfactor = 1.0f;
    }
      
    // Eroded sediment calculated from mobilised sediment, possibly limited by the transport capacity
    erodedsed = 1000.0f * mobilisedsed * transportfactor;  // kg/km2
*/ 
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