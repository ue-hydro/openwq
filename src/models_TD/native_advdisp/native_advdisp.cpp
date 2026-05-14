

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


#include "models_TD/headerfile_TD.hpp"

/* #################################################
// Mass transport
// Advection & Dispersion (Fickian dispersive flux approximation)
//
// The advective flux is: F_adv = (1 - exp(-Q*dt/V)) * chemass_source
//   which moves mass proportionally to the water flux (concentration-based).
//
// The dispersive flux uses a Fickian approximation between cell pairs with
// a harmonic-mean volume so the scheme stays stable and symmetric for any
// volume ratio:
//   V_h    = wmass_source * wmass_recipient / (wmass_source + wmass_recipient)
//   F_disp = D_eff * dt * (C_source - C_recipient) * V_h
// where:
//   D_eff = D_avg / L^2  [1/s]
//   D_avg = (Dx + Dy + Dz) / 3  [m2/s]  (averaged since we don't know connection orientation)
//   L = characteristic_length_m  [m]  (user-provided distance between cell centers)
//   C_source    = chemass_source    / wmass_source
//   C_recipient = chemass_recipient / wmass_recipient
//
// With V_h, the resulting concentration changes are
//   ΔC_source    = -D_eff * dt * (C_s - C_r) * V_r / (V_s + V_r)
//   ΔC_recipient = +D_eff * dt * (C_s - C_r) * V_s / (V_s + V_r)
// so the gradient simply relaxes by a factor of (1 - D_eff*dt) per step
// without overshoot, regardless of how asymmetric V_s and V_r are.
//
// Both fluxes are evaluated using the PRE-call state (start-of-call mass on
// both sides) and then summed before being applied. This avoids the previous
// failure mode where dispersion reacted to the artificial gradient that
// advection had just introduced and ended up moving mass back upstream,
// undoing advection.
//
// The dispersive flux is BIDIRECTIONAL: it can move mass in either direction
// depending on the concentration gradient. Sign convention: F_disp > 0 moves
// mass source -> recipient.
//
// OPTIMIZED: cached flags, pre-fetched refs, pre-computed D_eff
################################################# */
void OpenWQ_TD_model::AdvDisp(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    double wflux_s2r,
    double wmass_source){

    // Return if no flux: wflux_s2r == 0
    if(wflux_s2r == 0.0f){return;}

    // OPTIMIZED: use cached values instead of string comparisons
    const unsigned int numspec = OpenWQ_wqconfig.cached_num_mobile_species;
    const std::vector<unsigned int>& mobile_species = *OpenWQ_wqconfig.cached_mobile_species_ptr;

    // Pre-compute advective concentration factor using exponential decay
    // (analytical well-mixed CSTR solution): conc_factor = 1 - exp(-Q*dt/V)
    // For CFL << 1: ≈ CFL (same as linear upwind, difference < 5% for CFL < 0.1)
    // For CFL = 1:  0.632 (retains 36.8% of mass, prevents full cell flushing)
    // For CFL >> 1: → 1.0 (asymptotically flushes cell)
    // This eliminates alternating-zero oscillations when time step ≥ cell residence time
    const double conc_factor = 1.0 - std::exp(-wflux_s2r / wmass_source);

    // Get pre-computed effective dispersion rate [1/s]
    const double D_eff = OpenWQ_wqconfig.TD_model->NativeAdvDisp->D_eff;

    // Host model time step [s]. D_eff has units of 1/s, so multiplying by dt
    // gives a dimensionless mass fraction. Without dt the dispersive flux
    // silently used an implicit 1-second window regardless of the actual
    // routing time step, producing arbitrary magnitudes.
    const double dt = OpenWQ_hostModelconfig.get_time_step();

    // Get recipient water mass from host model for concentration calculation
    double wmass_recipient = 0.0;
    const bool has_recipient = (recipient != -1);
    if (has_recipient){
        wmass_recipient = OpenWQ_hostModelconfig.get_waterVol_hydromodel_at(
            recipient, ix_r, iy_r, iz_r);
    }

    // OPTIMIZED: pre-fetch field references to avoid repeated pointer dereferencing
    auto& chemass_source = (*OpenWQ_vars.chemass)(source);
    auto& d_transp_source = (*OpenWQ_vars.d_chemass_dt_transp)(source);

    // Loop for mobile chemical species
    for (unsigned int chemi=0;chemi<numspec;chemi++){

        const unsigned int ichem_mob = mobile_species[chemi];

        // ===========================
        // 0) Snapshot the START-OF-STEP mass in source and recipient.
        // Using chemass only (not chemass + d_transp) makes the scheme
        // order-independent: every AdvDisp call reads the same source and
        // recipient state regardless of which other reaches have been
        // processed earlier in the same timestep. This eliminates the
        // alternating ±large-flux oscillations that appeared whenever the
        // host (MPI/OpenMP/parallel routing) processed reaches out of
        // strict topological order. Trade-off: mass advances at most one
        // reach per timestep.
        // ===========================
        const double src_initial = std::fmax(
            chemass_source(ichem_mob)(ix_s,iy_s,iz_s), 0.0);

        double recip_initial = 0.0;
        if (has_recipient){
            recip_initial = std::fmax(
                (*OpenWQ_vars.chemass)(recipient)(ichem_mob)(ix_r,iy_r,iz_r), 0.0);
        }

        // ===========================
        // 1) ADVECTIVE FLUX (source -> recipient via water flux)
        // ===========================
        double chemass_flux_adv = conc_factor * src_initial;
        chemass_flux_adv = std::fmin(chemass_flux_adv, src_initial);

        // ===========================
        // 2) DISPERSIVE FLUX (Fickian, pre-call gradient, with dt)
        //
        // The mass exchanged uses the HARMONIC MEAN of source and recipient
        // volumes:
        //     V_h = V_s * V_r / (V_s + V_r)
        // This is the mathematically correct quantity for symmetric exchange
        // between two cells and keeps the scheme unconditionally stable for
        // any volume ratio. Using V_source alone (the textbook simplification
        // valid only when V_s == V_r) makes the concentration change at the
        // recipient scale as V_s/V_r, so when V_s >> V_r the recipient
        // overshoots the source's start-of-step concentration and dispersion
        // becomes an exponential amplifier instead of a smoother. mizuRoute
        // reach volumes can vary by orders of magnitude, so the symmetric
        // form is required for robust behaviour on any network.
        // ===========================
        double chemass_flux_disp = 0.0;
        if (has_recipient && D_eff > 0.0 && dt > 0.0 && wmass_recipient > 0.0){
            const double C_source    = src_initial   / wmass_source;
            const double C_recipient = recip_initial / wmass_recipient;
            const double V_harmonic  = (wmass_source * wmass_recipient)
                                     / (wmass_source + wmass_recipient);
            chemass_flux_disp = D_eff * dt * (C_source - C_recipient) * V_harmonic;
        }

        // ===========================
        // 3) Combined flux: cap at whichever side is donating mass so neither
        // side can be driven negative. The cap uses the LIVE balance
        // (chemass + d_transp) rather than start-of-step state, so that
        // successive calls in the same step (e.g. multiple upstream tribs
        // converging on the same recipient, all pulling backwards via
        // dispersion) cannot collectively drain more mass than the donating
        // side actually has. Using start-of-step for the cap permits N calls
        // to each pull recip_initial, creating N-1 copies of mass via the
        // end-of-step non-negativity clamp.
        // ===========================
        double chemass_flux_total = chemass_flux_adv + chemass_flux_disp;
        if (chemass_flux_total > 0.0){
            // Net transfer source -> recipient
            const double src_live = std::fmax(
                chemass_source(ichem_mob)(ix_s,iy_s,iz_s)
                + d_transp_source(ichem_mob)(ix_s,iy_s,iz_s), 0.0);
            chemass_flux_total = std::fmin(chemass_flux_total, src_live);
        } else if (chemass_flux_total < 0.0){
            // Net transfer recipient -> source (only possible with dispersion)
            const double recip_live = has_recipient ? std::fmax(
                (*OpenWQ_vars.chemass)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
                + (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r), 0.0)
                : 0.0;
            chemass_flux_total = std::fmax(chemass_flux_total, -recip_live);
        }

        // ===========================
        // 4) Apply
        // ===========================
        d_transp_source(ichem_mob)(ix_s,iy_s,iz_s) -= chemass_flux_total;

        if (!has_recipient){
            // OUT-flux: mass leaving the system. No dispersion in this branch
            // (no neighbour), so chemass_flux_total == chemass_flux_adv here.
            if (OpenWQ_vars.mass_balance.initialized &&
                ichem_mob < OpenWQ_vars.mass_balance.num_species) {
                OpenWQ_vars.mass_balance.cumulative_out_flux[ichem_mob] += chemass_flux_total;
            }
            continue;
        }
        (*OpenWQ_vars.d_chemass_dt_transp)(recipient)(ichem_mob)(ix_r,iy_r,iz_r)
            += chemass_flux_total;
    }

}
