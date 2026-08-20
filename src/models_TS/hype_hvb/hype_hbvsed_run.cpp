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

/* ########################################
// Calculate eroded particles from soil by HBV-sed based model
#########################################*/
void OpenWQ_TS_model::hbvsed_hype_erosion_run(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r, double wmass_source,
    std::string TS_type){


   // Local variables
  int exchange_direction = OpenWQ_wqconfig.TS_model->HypeHVB->get_exchange_direction();
  int erodFlux_icmp = OpenWQ_wqconfig.TS_model->HypeHVB->get_eroding_transp_compartment();
  int erodInhibut_icmp = OpenWQ_wqconfig.TS_model->HypeHVB->get_erosion_inhibit_compartment();
  int sediment_icmp = OpenWQ_wqconfig.TS_model->HypeHVB->get_sediment_compartment();


	// ###################
	// CHECK IF RUN APPLICABLE

  // only applicable with precipitation input
  if (TS_type.compare("TS_type_EWF") == 0 && recipient == sediment_icmp){

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    std::vector<unsigned int> xyz_lowerComp;

    // Get number of cells in input_direction_index
    std::vector<unsigned int> xyz_SedCompt;
    xyz_SedCompt.push_back(ix_r);
    xyz_SedCompt.push_back(iy_r);
    xyz_SedCompt.push_back(iz_r);         // get x,y,z from source comparment in vector

    // Ignore if current sediment source cell is not the top one
    if (xyz_SedCompt[exchange_direction] != 0)
      return;

    // Return if there is snow (erosion inhibited by snow/frost cover).
    // erodInhibut_icmp == -1 means EROSION_INHIBIT_COMPARTMENT = "NONE": the
    // inhibition is disabled (e.g. river-routing hosts with no snow compartment).
    // TODO: if (frostdepth<0. .AND. frostdepth>-9999.) then return
    // The gate reads the RECIPIENT cell's HRU (top layer): EWF calls carry
    // (0,0,0) as source coordinates, so using ix_s would always test HRU 0.
    double snow = 0.0;
    if (erodInhibut_icmp >= 0){
      snow = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
            erodInhibut_icmp,
            ix_r, iy_r, 0);
    }
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

    // HYPE's formula expects precipitation as a DEPTH (mm); openWQ passes the
    // water as a VOLUME (m3). Convert using the cell area, which the host
    // coupling publishes as a dependency variable named "cellArea_m2". Look it
    // up BY NAME (not a fixed index) so any host works regardless of how many
    // dependency variables it registers (SUMMA and mizuRoute differ). If no
    // host publishes it, fall back to the raw volume so behaviour is unchanged.
    double cellArea_m2 = 0.0;
    const int nDep = (int) OpenWQ_hostModelconfig.get_num_HydroDepend();
    for (int di = 0; di < nDep; di++){
        if (OpenWQ_hostModelconfig.get_HydroDepend_name_at(di).compare("cellArea_m2") == 0){
            cellArea_m2 = OpenWQ_hostModelconfig.get_dependVar_at(di, ix_r, iy_r, 0);
            break;
        }
    }
    double prec_mm = (cellArea_m2 > 0.0) ? (wflux_s2r / cellArea_m2 * 1000.0) : wflux_s2r;

    mobilisedsed = pow((slope / 5.),slopepar)
                  * lusepar * soilpar * erosionindex
                  * pow(prec_mm,precexppar); // !tonnes/km2 = g/m2
    
    // Eroded sediment calculated from mobilised sediment with a monthly factor
    // From HYPE: erodmonth = 1. + monthpar(m_erodmon, current_time%date%month)
    {
      time_t simtime = OpenWQ_wqconfig.get_simtime();
      struct tm tm_sim;
      gmtime_r(&simtime, &tm_sim);
      int current_month = tm_sim.tm_mon; // 0-based (Jan=0, Dec=11)
      erodmonth = 1.0 + OpenWQ_wqconfig.TS_model->HypeHVB->monthpar[current_month];
    }
    // Convert the areal density (g/m2) to an ABSOLUTE mass (kg) for this cell:
    //   g/m2 * area_m2 = g ; / 1000 = kg.
    // If cell area is unavailable, keep the legacy 1000x scaling so results are
    // unchanged for couplings that do not publish cell area.
    double abs_scale = (cellArea_m2 > 0.0) ? (cellArea_m2 / 1000.0) : 1000.0;
    (*OpenWQ_vars.d_sedmass_mobilized_dt)(ix_r, iy_r, iz_r) = abs_scale * mobilisedsed * erodmonth;  // kg
  
  // #######################
  // Mobind with runoff if flow exists
  } else if (TS_type.compare("TS_type_LE") == 0
      && source == erodFlux_icmp
      && source == recipient
      && ix_r!=-1 
      && iy_r!=-1 
      && iz_r!=-1){

    /// 2)
	  // Mobilization with flow
	  // Assuming here that, since d_sedmass_mobilized_dt is the mobilized sediment,
	  // then it will all move with flow regardless of the flow/runoff intensity
    // [conservation fix] Do NOT fold the mobilized (eroded) sediment into d_sedmass_transport_dt.
    // d_sedmass_transport_dt carries ONLY transport (below); the solver adds d_sedmass_mobilized_dt
    // (the erosion source) separately, so the eroded mass is counted exactly once.

    // Advect sediment with the flow (residence-time-limited): move only the
    // fraction of this cell's suspended sediment that actually flows out this
    // step, frac = Qout*dt / (cell water volume). Moving ALL of it regardless of
    // flow (the previous scheme) left the routed transport undamped and made the
    // reach sediment oscillate. frac==1 recovers the full flush (e.g. SUMMA
    // runoff, where wmass_source == wflux_s2r).
    {
      double frac = (wmass_source > 0.0) ? (wflux_s2r / wmass_source) : 0.0;
      if (frac > 1.0) frac = 1.0;
      if (frac < 0.0) frac = 0.0;
      double moved = (*OpenWQ_vars.sedmass)(ix_s, iy_s, iz_s) * frac;
      // removing from source
      (*OpenWQ_vars.d_sedmass_transport_dt)(ix_s, iy_s, iz_s) -= moved;
      // adding to recipient
      (*OpenWQ_vars.d_sedmass_transport_dt)(ix_r, iy_r, iz_r) += moved;

      // Sorbed partner species (species-pair scheme, e.g. PO4-P => PP) are
      // particle-bound, so they ride the sediment: move the SAME fraction of
      // their chemass, through the existing transport derivative channel
      // (integrated by the chemistry solver). x% of sedmass moved => x% of
      // the sorbed species moved.
      {
        const unsigned int npairs = OpenWQ_wqconfig.SI_model->get_num_sorption_pairs();
        for (unsigned int pi = 0; pi < npairs; pi++){
            const int sorbi = OpenWQ_wqconfig.SI_model->get_sorbed_species_index_at(pi);
            if (sorbi < 0) continue;
            double smoved = std::fmax(OpenWQ_vars.live_mass(
                erodFlux_icmp, (unsigned int) sorbi, ix_s, iy_s, iz_s), 0.0) * frac;
            (*OpenWQ_vars.d_chemass_dt_transp_part)(erodFlux_icmp)(sorbi)(ix_s, iy_s, iz_s) -= smoved;
            (*OpenWQ_vars.d_chemass_dt_transp_part)(erodFlux_icmp)(sorbi)(ix_r, iy_r, iz_r) += smoved;
        }
      }
    }

  // #######################
  // Sediment leaving the domain with surface runoff (recipient == -1 = OUT)
  } else if (TS_type.compare("TS_type_LE") == 0
      && source == erodFlux_icmp
      && recipient == -1
      && wflux_s2r > 0.0
      && ix_s!=-1
      && iy_s!=-1
      && iz_s!=-1){

    // Surface runoff exiting the system carries its suspended sediment out.
    // Same residence-time-limited fraction as the reach-to-reach advection
    // above: only the sediment carried by the outflow leaves this step.
    {
      double frac = (wmass_source > 0.0) ? (wflux_s2r / wmass_source) : 0.0;
      if (frac > 1.0) frac = 1.0;
      if (frac < 0.0) frac = 0.0;
      (*OpenWQ_vars.d_sedmass_transport_dt)(ix_s, iy_s, iz_s) -= (*OpenWQ_vars.sedmass)(ix_s, iy_s, iz_s) * frac;

      // Sorbed partner species leave the domain with the sediment they are
      // bound to, at the same fraction (outlet sink; no recipient).
      {
        const unsigned int npairs = OpenWQ_wqconfig.SI_model->get_num_sorption_pairs();
        for (unsigned int pi = 0; pi < npairs; pi++){
            const int sorbi = OpenWQ_wqconfig.SI_model->get_sorbed_species_index_at(pi);
            if (sorbi < 0) continue;
            (*OpenWQ_vars.d_chemass_dt_transp_part)(erodFlux_icmp)(sorbi)(ix_s, iy_s, iz_s) -=
                std::fmax(OpenWQ_vars.live_mass(
                    erodFlux_icmp, (unsigned int) sorbi, ix_s, iy_s, iz_s), 0.0) * frac;
        }
      }
    }

  }
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