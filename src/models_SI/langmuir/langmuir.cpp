
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

#include "models_SI/headerfile_SI.hpp"

/* #################################################
// Langmuir Sorption Isotherm — SPECIES-PAIR scheme (dissolved => sorbed)
// q = (qmax * KL * C) / (1 + KL * C)
//
// Each sorbable species is PAIRED with an explicit sorbed partner species
// (e.g. PO4-P => PP), both regular chemass state variables. The isotherm
// shifts mass between the two, BGC-style, via its own derivative channel:
//     d_chemass_dt_sorpt(dissolved) -= flux
//     d_chemass_dt_sorpt(sorbed)    += flux
//
// q_current comes from the REAL sorbed pool (the partner species' chemass),
// making this a true kinetic dissolved<=>sorbed exchange. Equilibrium is the
// partition of the REAL total (dissolved + sorbed) mass:
//     total_mass = C * Vol + q * Msolid_g
// which leads to a quadratic in C:
//     KL * Vol * C^2 + (Vol + qmax * KL * Msolid_g - KL * total_mass) * C - total_mass = 0
// where Msolid_g = bulk_density * layer_thickness * cellArea / 1000
// (g of dissolved mass per mg/kg of sorbed conc; cellArea from the host
//  coupling's "cellArea_m2" dependency, fallback 1 m2 = legacy behaviour).
//
// Kinetic adsorption/desorption toward that equilibrium:
//     flux = (q_eq - q_current) * (1 - exp(-Kadsdes * dt)) * Msolid_g
//
// The sorbed partner is particle-bound: excluded from water advection
// (keep it out of BGC_GENERAL_MOBILE_SPECIES) and co-transported with the
// sediment by model_TS at the sediment transport fraction.
################################################# */

void OpenWQ_SI_model::langmuir(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_output& OpenWQ_output){

    // Get Langmuir parameters
    auto* LG = OpenWQ_wqconfig.SI_model->LANGMUIR;
    const unsigned int num_si_species = LG->num_species;

    if (num_si_species == 0) return;  // Nothing to do

    const double bulk_density = LG->bulk_density;    // kg/m3
    const double layer_thick = LG->layer_thickness;  // m
    const double rhoL = bulk_density * layer_thick;  // kg/m2

    // Get timestep for kinetic rate
    const double dt = OpenWQ_hostModelconfig.get_time_step();  // seconds
    if (dt <= 0.0) return;

    // Get number of compartments
    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Water volume minimum limit
    const double watervol_minlim = OpenWQ_hostModelconfig.get_watervol_minlim();

    // Locate the host coupling's cell-area dependency ("cellArea_m2") by
    // name, so the solid-phase mass basis is a real absolute mass. If the
    // host doesn't publish it, fall back to 1 m2 (legacy behaviour).
    int cellArea_dep = -1;
    {
        const int nDep = (int) OpenWQ_hostModelconfig.get_num_HydroDepend();
        for (int di = 0; di < nDep; di++){
            if (OpenWQ_hostModelconfig.get_HydroDepend_name_at(di).compare("cellArea_m2") == 0){
                cellArea_dep = di;
                break;
            }
        }
    }

    // Sediment-transport compartment index (if model_TS is active): there the
    // solid phase is the ACTUAL suspended sediment (sedmass), so sorption is
    // sediment-limited - no particles in suspension means nothing to sorb onto.
    int transp_icmp = -1;
    if (OpenWQ_wqconfig.is_TS_enabled) {
        const std::string& transp_name = OpenWQ_wqconfig.TS_model->ErodTranspCmpt;
        for (unsigned int ic = 0; ic < num_comps; ic++){
            if (transp_name.compare(OpenWQ_hostModelconfig.get_HydroComp_name_at(ic)) == 0){
                transp_icmp = (int) ic;
                break;
            }
        }
    }

    // Loop over each SI species pair (dissolved => sorbed)
    for (unsigned int si = 0; si < num_si_species; si++) {

        const unsigned int chemi = LG->species_index[si];
        const int sorbi = LG->sorbed_species_index[si];
        const double qmax = LG->qmax[si];
        const double KL = LG->KL[si];
        const double Kadsdes = LG->Kadsdes[si];

        // Skip if parameters are invalid or the pair is unresolved
        if (qmax <= 0.0 || KL <= 0.0 || Kadsdes <= 0.0 || sorbi < 0) continue;

        // Pre-compute kinetic factor: (1 - exp(-Kadsdes * dt))
        const double kinetic_factor = 1.0 - std::exp(-Kadsdes * dt);

        // Loop over all compartments and cells
        for (unsigned int icmp = 0; icmp < num_comps; icmp++) {

            const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
            const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
            const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

            // References to the sorption derivative cubes (the species pair)
            auto& d_chemass = (*OpenWQ_vars.d_chemass_dt_sorpt)(icmp)(chemi);
            auto& d_chemass_sorb = (*OpenWQ_vars.d_chemass_dt_sorpt)(icmp)(sorbi);

            for (unsigned int ix = 0; ix < nx; ix++) {
                for (unsigned int iy = 0; iy < ny; iy++) {
                    for (unsigned int iz = 0; iz < nz; iz++) {

                        // Get water volume [m3] in this cell
                        const double Vol = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
                            icmp, ix, iy, iz);

                        // Skip dry cells
                        if (Vol <= watervol_minlim) continue;

                        // Current dissolved and sorbed mass [g] in this cell,
                        // from the LIVE balance (state + claims already written
                        // this step, e.g. BGC reactions) so the sorption flux
                        // cannot overdraw a pool another process has claimed.
                        const double mass_dissolved = std::fmax(
                            OpenWQ_vars.live_mass(icmp, chemi, ix, iy, iz), 0.0);
                        const double mass_sorbed = std::fmax(
                            OpenWQ_vars.live_mass(icmp, (unsigned int) sorbi, ix, iy, iz), 0.0);
                        if (mass_dissolved <= 0.0 && mass_sorbed <= 0.0) continue;

                        // Solid-phase mass basis for this cell:
                        // Msolid_g converts q [mg/kg] to an absolute mass [g].
                        double Msolid_g;
                        if ((int) icmp == transp_icmp) {
                            // Transport compartment: the sorbent is the ACTUAL
                            // suspended sediment. q [mg/kg] * sedmass [kg] = mg ; /1000 = g
                            Msolid_g = (*OpenWQ_vars.sedmass)(ix, iy, iz) / 1000.0;
                            if (Msolid_g <= 0.0) {
                                // No particles in suspension: nothing to sorb
                                // onto; sorbed mass present desorbs at the
                                // kinetic rate (paired, conservative).
                                if (mass_sorbed > 0.0) {
                                    const double released = mass_sorbed * kinetic_factor;
                                    d_chemass(ix, iy, iz) += released;
                                    d_chemass_sorb(ix, iy, iz) -= released;
                                }
                                continue;
                            }
                        } else {
                            // Soil-matrix compartments:
                            //   q [mg/kg] * (rhoL [kg/m2] * area [m2]) = mg ; /1000 = g
                            double cellArea_m2 = 1.0;
                            if (cellArea_dep >= 0){
                                const double a = OpenWQ_hostModelconfig.get_dependVar_at(
                                    cellArea_dep, ix, iy, 0);
                                if (a > 0.0) cellArea_m2 = a;
                            }
                            Msolid_g = rhoL * cellArea_m2 / 1000.0;
                            if (Msolid_g <= 0.0) continue;
                        }

                        // Current sorbed concentration from the REAL pool [mg/kg]
                        const double q_current = (mass_sorbed > 0.0)
                            ? (mass_sorbed / Msolid_g) : 0.0;

                        // Total REAL mass in this cell (dissolved + sorbed) [g]
                        const double total_mass = mass_dissolved + mass_sorbed;

                        if (total_mass <= 0.0) continue;

                        // ############################
                        // Analytical solution for equilibrium C_eq via quadratic formula
                        //
                        // From mass conservation with Langmuir:
                        //   C * Vol + (qmax * KL * C / (1 + KL * C)) * Msolid_g = total_mass
                        //
                        // Multiply through by (1 + KL * C):
                        //   C * Vol * (1 + KL * C) + qmax * KL * C * Msolid_g = total_mass * (1 + KL * C)
                        //   C * Vol + KL * Vol * C^2 + qmax * KL * Msolid_g * C = total_mass + KL * total_mass * C
                        //
                        // Rearranging:
                        //   KL * Vol * C^2 + (Vol + qmax * KL * Msolid_g - KL * total_mass) * C - total_mass = 0
                        //
                        // a = KL * Vol
                        // b = Vol + qmax * KL * Msolid_g - KL * total_mass
                        // c = -total_mass
                        // ############################

                        const double a = KL * Vol;
                        const double b = Vol + qmax * KL * Msolid_g - KL * total_mass;
                        const double c = -total_mass;

                        const double discriminant = b * b - 4.0 * a * c;

                        double C_eq;
                        if (discriminant < 0.0 || a <= 0.0) {
                            // Fallback: assume all dissolved
                            C_eq = total_mass / Vol;
                        } else {
                            // Take positive root
                            C_eq = (-b + std::sqrt(discriminant)) / (2.0 * a);
                            if (C_eq < 0.0) {
                                C_eq = (-b - std::sqrt(discriminant)) / (2.0 * a);
                            }
                            if (C_eq < 0.0) {
                                C_eq = total_mass / Vol;  // Fallback
                            }
                        }

                        // Equilibrium sorbed concentration
                        const double q_eq = (qmax * KL * C_eq) / (1.0 + KL * C_eq);

                        // ############################
                        // Kinetic adsorption/desorption
                        // Mass flux from dissolved to sorbed phase
                        // positive = adsorption (remove from dissolved)
                        // negative = desorption (add to dissolved)
                        // ############################

                        const double delta_q = (q_eq - q_current) * kinetic_factor;
                        double flux_mass = delta_q * Msolid_g;  // [g]

                        // Safety: adsorption cannot remove more than the dissolved mass
                        if (flux_mass > 0.0 && flux_mass > mass_dissolved) {
                            flux_mass = mass_dissolved;
                        }

                        // Safety: desorption cannot remove more than the sorbed mass
                        if (flux_mass < 0.0 && -flux_mass > mass_sorbed) {
                            flux_mass = -mass_sorbed;
                        }

                        // Apply as a BGC-style shift between the species pair
                        // (positive flux_mass = adsorption: dissolved -> sorbed)
                        d_chemass(ix, iy, iz) -= flux_mass;
                        d_chemass_sorb(ix, iy, iz) += flux_mass;

                    } // iz
                } // iy
            } // ix
        } // icmp
    } // si species
}
