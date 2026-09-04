// Copyright 2026, Diogo Costa (diogo.costa@uevora.pt)
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


#include "compute/headerfile_compute.hpp"
#include "output/headerfile_OUT.hpp"   // OpenWQ_output::ConsoleLog (ML_CLOSURES logging)
#include <cctype>   // std::toupper (ML_CLOSURES name resolution)


/* #################################################
// Forward Euler (explicit) solver
################################################# */

void OpenWQ_compute::Solve_with_ForwardEuler(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output,
    OpenWQ_CH_model& OpenWQ_CH_model){

    // OPTIMIZED: use cached flag instead of string comparison
    const unsigned int num_chem = OpenWQ_wqconfig.is_native_bgc_flex
        ? OpenWQ_wqconfig.CH_model->NativeFlex->num_chem
        : OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();
    const bool is_first_step = OpenWQ_hostModelconfig.is_first_interaction_step();

    // Hybrid physics-ML LAYER 2 (per-species derivative closures): resolve names
    // -> indices + build the conservation groups once (single-threaded, before
    // the parallel region). Empty -> ml_on false -> the original solver path.
    if (!OpenWQ_wqconfig.ml_closures_ready)
        Prepare_MLClosures(OpenWQ_hostModelconfig, OpenWQ_wqconfig, OpenWQ_output);
    const bool ml_on = !OpenWQ_wqconfig.ml_deriv_closures.empty();

    /* #####################################################
    // Compartment loop - parallelized over compartments and chemicals
    ##################################################### */

    #pragma omp parallel num_threads(num_threads)
    {
        // Local variables - thread private
        unsigned int nx, ny, nz;
        unsigned int ix, iy, iz;
        double dm_dt_chem, dm_dt_sorpt, dm_dt_trans, dm_dt_part, dm_ic, dm_ss, dm_ewf;

        #pragma omp for schedule(dynamic) collapse(2)
        for (unsigned int icmp = 0; icmp < num_comps; icmp++){
            for (unsigned int chemi = 0; chemi < num_chem; chemi++){

                // Get dimensions for this compartment (only once per icmp-chemi pair)
                nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
                ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
                nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

                // Pre-fetch pointers to reduce dereferencing overhead
                auto& chemass = (*OpenWQ_vars.chemass)(icmp)(chemi);
                auto& d_chemass_ic = (*OpenWQ_vars.d_chemass_ic)(icmp)(chemi);
                auto& d_chemass_ss = (*OpenWQ_vars.d_chemass_ss)(icmp)(chemi);
                auto& d_chemass_ss_out = (*OpenWQ_vars.d_chemass_ss_out)(icmp)(chemi);
                auto& d_chemass_ewf = (*OpenWQ_vars.d_chemass_ewf)(icmp)(chemi);
                auto& d_chemass_ewf_out = (*OpenWQ_vars.d_chemass_ewf_out)(icmp)(chemi);
                auto& d_chemass_dt_chem = (*OpenWQ_vars.d_chemass_dt_chem)(icmp)(chemi);
                auto& d_chemass_dt_chem_out = (*OpenWQ_vars.d_chemass_dt_chem_out)(icmp)(chemi);
                auto& d_chemass_dt_sorpt = (*OpenWQ_vars.d_chemass_dt_sorpt)(icmp)(chemi);
                auto& d_chemass_dt_sorpt_out = (*OpenWQ_vars.d_chemass_dt_sorpt_out)(icmp)(chemi);
                auto& d_chemass_dt_transp_diss = (*OpenWQ_vars.d_chemass_dt_transp_diss)(icmp)(chemi);
                auto& d_chemass_dt_transp_diss_out = (*OpenWQ_vars.d_chemass_dt_transp_diss_out)(icmp)(chemi);
                auto& d_chemass_dt_transp_part = (*OpenWQ_vars.d_chemass_dt_transp_part)(icmp)(chemi);
                auto& d_chemass_dt_transp_part_out = (*OpenWQ_vars.d_chemass_dt_transp_part_out)(icmp)(chemi);

                // LAYER 2 derivative-closure lookup for this (icmp, chemi). A
                // closure-driven species takes the ML loop below; every other
                // species keeps the original vectorized loop -> byte-identical.
                int _mlc_chem = -1, _mlc_sorpt = -1, _mlc_ss = -1;
                if (ml_on) {
                    _mlc_chem  = OpenWQ_wqconfig.ml_chem_cl[icmp][chemi];
                    _mlc_sorpt = OpenWQ_wqconfig.ml_sorpt_cl[icmp][chemi];
                    _mlc_ss    = OpenWQ_wqconfig.ml_ss_cl[icmp][chemi];
                }
                const bool _use_ml = (_mlc_chem >= 0 || _mlc_sorpt >= 0 || _mlc_ss >= 0);

                // Case-2 diagnostics: thread-local factor accumulators for the
                // up-to-3 active closures at this (icmp,chemi), merged into the
                // shared stats ONCE after the grid loop (no per-cell atomics).
                OpenWQ_wqconfig::MLClosureStats _acc_chem, _acc_sorpt, _acc_ss;
                if (_mlc_chem  >= 0) _acc_chem.max_correction  = OpenWQ_wqconfig.ml_deriv_closures[_mlc_chem ].net.max_correction;
                if (_mlc_sorpt >= 0) _acc_sorpt.max_correction = OpenWQ_wqconfig.ml_deriv_closures[_mlc_sorpt].net.max_correction;
                if (_mlc_ss    >= 0) _acc_ss.max_correction    = OpenWQ_wqconfig.ml_deriv_closures[_mlc_ss   ].net.max_correction;

                // X, Y, Z loops - process entire 3D grid for this compartment-chemical pair
                for (ix = 0; ix < nx; ix++){
                    for (iy = 0; iy < ny; iy++){
                      if (!_use_ml) {
                        // ---- original pure-physics path (unchanged) ----
                        // Vectorization hint for innermost loop
                        #pragma omp simd
                        for (iz = 0; iz < nz; iz++){

                            // ####################################
                            // 1. Single-time change at simulation start
                            dm_ic = is_first_step ? d_chemass_ic(ix, iy, iz) : 0.0;

                            // ####################################
                            // 2. SS (Sink & Sources)
                            dm_ss = d_chemass_ss(ix, iy, iz);
                            d_chemass_ss_out(ix, iy, iz) += dm_ss;

                            // ####################################
                            // 3. EWF (External Water Fluxes)
                            dm_ewf = d_chemass_ewf(ix, iy, iz);
                            d_chemass_ewf_out(ix, iy, iz) += dm_ewf;

                            // ####################################
                            // 4. Dynamic change (derivatives): chemistry and transport
                            dm_dt_chem = d_chemass_dt_chem(ix, iy, iz);
                            d_chemass_dt_chem_out(ix, iy, iz) += dm_dt_chem;

                            dm_dt_sorpt = d_chemass_dt_sorpt(ix, iy, iz);
                            d_chemass_dt_sorpt_out(ix, iy, iz) += dm_dt_sorpt;

                            dm_dt_trans = d_chemass_dt_transp_diss(ix, iy, iz);
                            d_chemass_dt_transp_diss_out(ix, iy, iz) += dm_dt_trans;

                            dm_dt_part = d_chemass_dt_transp_part(ix, iy, iz);
                            d_chemass_dt_transp_part_out(ix, iy, iz) += dm_dt_part;

                            // ####################################
                            // 5. Apply all changes to state variable
                            double new_mass = chemass(ix, iy, iz) + dm_ic + dm_ss + dm_ewf
                                              + dm_dt_chem + dm_dt_sorpt + dm_dt_trans + dm_dt_part;

                            // Ensure non-negativity and track mass lost to clamping
                            if (new_mass > 0.0) {
                                chemass(ix, iy, iz) = new_mass;
                            } else {
                                // Track mass lost to negative correction for mass balance
                                if (OpenWQ_vars.mass_balance.initialized &&
                                    chemi < OpenWQ_vars.mass_balance.num_species) {
                                    // Note: This is not thread-safe, but the error is small
                                    // For exact tracking, would need atomic or per-thread accumulation
                                    #pragma omp atomic
                                    OpenWQ_vars.mass_balance.cumulative_negative_correction[chemi] -= new_mass;
                                }
                                chemass(ix, iy, iz) = 0.0;
                            }
                        }
                      } else {
                        // ---- LAYER 2: bounded per-species derivative closures ----
                        // Same integration, but CHEM/SORPT/SS derivatives are
                        // scaled by their bounded factor first. CHEM/SORPT use a
                        // group-shared factor (conservative); the factor input is
                        // the driving (observed) species' state at this cell.
                        for (iz = 0; iz < nz; iz++){

                            dm_ic = is_first_step ? d_chemass_ic(ix, iy, iz) : 0.0;
                            dm_ss = d_chemass_ss(ix, iy, iz);
                            dm_ewf = d_chemass_ewf(ix, iy, iz);
                            dm_dt_chem = d_chemass_dt_chem(ix, iy, iz);
                            dm_dt_sorpt = d_chemass_dt_sorpt(ix, iy, iz);
                            dm_dt_trans = d_chemass_dt_transp_diss(ix, iy, iz);
                            dm_dt_part = d_chemass_dt_transp_part(ix, iy, iz);

                            if (_mlc_chem >= 0) {
                                const auto& _cl = OpenWQ_wqconfig.ml_deriv_closures[_mlc_chem];
                                const double _f = _cl.net.factor(arma::vec({
                                    (*OpenWQ_vars.chemass)(icmp)(_cl.chemi)(ix, iy, iz)}));
                                _acc_chem.add(_f, dm_dt_chem);   // w = unscaled tendency
                                dm_dt_chem *= _f;
                            }
                            if (_mlc_sorpt >= 0) {
                                const auto& _cl = OpenWQ_wqconfig.ml_deriv_closures[_mlc_sorpt];
                                const double _f = _cl.net.factor(arma::vec({
                                    (*OpenWQ_vars.chemass)(icmp)(_cl.chemi)(ix, iy, iz)}));
                                _acc_sorpt.add(_f, dm_dt_sorpt);
                                dm_dt_sorpt *= _f;
                            }
                            if (_mlc_ss >= 0) {
                                const auto& _cl = OpenWQ_wqconfig.ml_deriv_closures[_mlc_ss];
                                const double _f = _cl.net.factor(arma::vec({
                                    (*OpenWQ_vars.chemass)(icmp)(_cl.chemi)(ix, iy, iz)}));
                                _acc_ss.add(_f, dm_ss);
                                dm_ss *= _f;
                            }

                            d_chemass_ss_out(ix, iy, iz) += dm_ss;
                            d_chemass_ewf_out(ix, iy, iz) += dm_ewf;
                            d_chemass_dt_chem_out(ix, iy, iz) += dm_dt_chem;
                            d_chemass_dt_sorpt_out(ix, iy, iz) += dm_dt_sorpt;
                            d_chemass_dt_transp_diss_out(ix, iy, iz) += dm_dt_trans;
                            d_chemass_dt_transp_part_out(ix, iy, iz) += dm_dt_part;

                            double new_mass = chemass(ix, iy, iz) + dm_ic + dm_ss + dm_ewf
                                              + dm_dt_chem + dm_dt_sorpt + dm_dt_trans + dm_dt_part;
                            if (new_mass > 0.0) {
                                chemass(ix, iy, iz) = new_mass;
                            } else {
                                if (OpenWQ_vars.mass_balance.initialized &&
                                    chemi < OpenWQ_vars.mass_balance.num_species) {
                                    #pragma omp atomic
                                    OpenWQ_vars.mass_balance.cumulative_negative_correction[chemi] -= new_mass;
                                }
                                chemass(ix, iy, iz) = 0.0;
                            }
                        }
                      }
                    }
                }

                // Case-2 diagnostics: fold this (icmp,chemi)'s thread-local
                // factor stats into the shared per-closure accumulators.
                if (_use_ml) {
                    #pragma omp critical(mlclosure_stats)
                    {
                        if (_mlc_chem  >= 0) OpenWQ_wqconfig.ml_closure_stats[_mlc_chem ].merge(_acc_chem);
                        if (_mlc_sorpt >= 0) OpenWQ_wqconfig.ml_closure_stats[_mlc_sorpt].merge(_acc_sorpt);
                        if (_mlc_ss    >= 0) OpenWQ_wqconfig.ml_closure_stats[_mlc_ss   ].merge(_acc_ss);
                    }
                }
            }
        }
    }

}

void OpenWQ_compute::Solve_with_ForwardEuler_Sediment(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_json& OpenWQ_json,
    OpenWQ_output& OpenWQ_output){

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();
    const unsigned int num_threads = OpenWQ_wqconfig.get_num_threads_requested();

    // Find sediment compartment index
    unsigned int sed_icmp = 0;
    const std::string& erod_transp_cmpt = OpenWQ_wqconfig.TS_model->ErodTranspCmpt;
    
    for (unsigned int icmp = 0; icmp < num_comps; icmp++){
        if (erod_transp_cmpt.compare(OpenWQ_hostModelconfig.get_HydroComp_name_at(icmp)) == 0) {
            sed_icmp = icmp;
            break;
        }
    }

    // Get dimensions for sediment compartment
    const unsigned int nx = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(sed_icmp);
    const unsigned int ny = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(sed_icmp);
    const unsigned int nz = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(sed_icmp);

    // Pre-fetch array references
    auto& sedmass = *OpenWQ_vars.sedmass;
    auto& d_sedmass_transport_dt = *OpenWQ_vars.d_sedmass_transport_dt;
    auto& d_sedmass_mobilized_dt = *OpenWQ_vars.d_sedmass_mobilized_dt;
    auto& d_sedmass_transport_dt_out = *OpenWQ_vars.d_sedmass_transport_dt_out;
    auto& d_sedmass_mobilized_dt_out = *OpenWQ_vars.d_sedmass_mobilized_dt_out;

    /* #####################################################
    // Parallel 3D loop with better granularity
    ##################################################### */

    #pragma omp parallel num_threads(num_threads)
    {
        unsigned int ix, iy, iz;
        double dm_dt;

        // Parallelize over flattened 3D index for better load balancing
        #pragma omp for schedule(static)
        for (unsigned int idx = 0; idx < nx * ny * nz; idx++){
            ix = idx / (ny * nz);
            iy = (idx / nz) % ny;
            iz = idx % nz;

            dm_dt = d_sedmass_transport_dt(ix, iy, iz) + d_sedmass_mobilized_dt(ix, iy, iz);

            // Accumulate for debug outputs
            d_sedmass_transport_dt_out(ix, iy, iz) += d_sedmass_transport_dt(ix, iy, iz);
            d_sedmass_mobilized_dt_out(ix, iy, iz) += d_sedmass_mobilized_dt(ix, iy, iz);

            // Apply change and ensure non-negativity
            double new_sedmass = sedmass(ix, iy, iz) + dm_dt;
            sedmass(ix, iy, iz) = (new_sedmass > 0.0) ? new_sedmass : 0.0;
        }
    }
}


// #################################################
// Hybrid physics-ML LAYER 2 — prepare the per-species derivative closures.
// Resolve config (names -> indices), compute the CHEM/SORPT conservation groups
// (union-find over the reaction graph / sorption pairs), and build the
// per-(icmp,chemi) lookup arrays. Runs ONCE (guarded by ml_closures_ready).
// Conservation: a group shares one factor, so every between-species transfer
// stays balanced. SS is external forcing (per-species).
// #################################################
void OpenWQ_compute::Prepare_MLClosures(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output){

    OpenWQ_wqconfig.ml_closures_ready = true;    // mark done regardless of outcome

    if (OpenWQ_wqconfig.ml_closure_cfg.empty())
        return;                                  // no closures -> pure physics

    auto _upper = [](std::string s){
        for (char& c : s) c = (char)std::toupper((unsigned char)c);
        return s; };
    std::string _m;   // reusable message buffer (ConsoleLog takes std::string&)

    const unsigned int num_comps = OpenWQ_hostModelconfig.get_num_HydroComp();

    // Conservation groups require the native BGC-flex reaction graph.
    if (!OpenWQ_wqconfig.is_native_bgc_flex) {
        _m = "<OpenWQ> WARNING: ML_CLOSURES require the native BGC-flex "
             "chemistry (reaction graph) - closures ignored.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, _m, true, true);
        return;
    }

    auto& nf = *(OpenWQ_wqconfig.CH_model->NativeFlex);
    const unsigned int num_chem = nf.num_chem;
    const std::vector<std::string>& species = nf.chem_species_list;

    // ---- union-find over the species set ----
    std::vector<int> parent(num_chem);
    auto uf_reset = [&](){ for (unsigned i=0;i<num_chem;i++) parent[i]=(int)i; };
    auto uf_find  = [&](int x){ while(parent[x]!=x){ parent[x]=parent[parent[x]]; x=parent[x]; } return x; };
    auto uf_union = [&](int a,int b){ int ra=uf_find(a), rb=uf_find(b); if(ra!=rb) parent[ra]=rb; };

    // CHEM groups: an edge between the consumed and produced species of every
    // transformation (an empty side is stored as UINT_MAX -> no edge).
    uf_reset();
    for (const auto& t : nf.BGCexpressions_info){
        unsigned int c = std::get<3>(t);
        unsigned int p = std::get<4>(t);
        if (c < num_chem && p < num_chem) uf_union((int)c,(int)p);
    }
    OpenWQ_wqconfig.ml_chem_group.assign(num_chem, 0);
    for (unsigned i=0;i<num_chem;i++) OpenWQ_wqconfig.ml_chem_group[i]=uf_find((int)i);

    // SORPT groups: dissolved <-> sorbed pairs (only if SI is active).
    uf_reset();
    if (OpenWQ_wqconfig.SI_model){
        for (unsigned int pi=0; pi<OpenWQ_wqconfig.SI_model->get_num_sorption_pairs(); pi++){
            int d = OpenWQ_wqconfig.SI_model->get_dissolved_species_index_at(pi);
            int s = OpenWQ_wqconfig.SI_model->get_sorbed_species_index_at(pi);
            if (d>=0 && s>=0 && (unsigned)d<num_chem && (unsigned)s<num_chem) uf_union(d,s);
        }
    }
    OpenWQ_wqconfig.ml_sorpt_group.assign(num_chem, 0);
    for (unsigned i=0;i<num_chem;i++) OpenWQ_wqconfig.ml_sorpt_group[i]=uf_find((int)i);

    // ---- resolve config -> ml_deriv_closures ----
    OpenWQ_wqconfig.ml_deriv_closures.clear();
    for (const auto& cfg : OpenWQ_wqconfig.ml_closure_cfg){
        const std::string term = _upper(cfg.term);
        if (term!="CHEM" && term!="SORPT" && term!="SS"){
            _m = "<OpenWQ> WARNING: ML_CLOSURES unknown TERM '"+cfg.term+"' - skipped.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, _m, true, true);
            continue;
        }
        int sidx=-1; const std::string su=_upper(cfg.species);
        for (unsigned j=0;j<num_chem;j++){ if(_upper(species[j])==su){ sidx=(int)j; break; } }
        if (sidx<0){
            _m = "<OpenWQ> WARNING: ML_CLOSURES species '"+cfg.species+"' not found - skipped.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, _m, true, true);
            continue;
        }
        std::vector<unsigned int> comps;
        const std::string cu=_upper(cfg.compartment);
        if (cu=="ALL" || cu==""){ for(unsigned c=0;c<num_comps;c++) comps.push_back(c); }
        else { for(unsigned c=0;c<num_comps;c++){ if(_upper(OpenWQ_hostModelconfig.get_HydroComp_name_at((int)c))==cu) comps.push_back(c); } }
        if (comps.empty()){
            _m = "<OpenWQ> WARNING: ML_CLOSURES compartment '"+cfg.compartment+"' not found - skipped.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, _m, true, true);
            continue;
        }
        for (unsigned int ic : comps){
            OpenWQ_wqconfig::MLDerivClosure mdc;
            mdc.icmp=ic; mdc.chemi=(unsigned)sidx; mdc.term=term; mdc.net=cfg.net;
            OpenWQ_wqconfig.ml_deriv_closures.push_back(mdc);
        }
    }

    // ---- build per-(icmp,chemi) lookup arrays (-1 = no closure) ----
    OpenWQ_wqconfig.ml_chem_cl.assign(num_comps, std::vector<int>(num_chem,-1));
    OpenWQ_wqconfig.ml_sorpt_cl.assign(num_comps, std::vector<int>(num_chem,-1));
    OpenWQ_wqconfig.ml_ss_cl.assign(num_comps, std::vector<int>(num_chem,-1));
    for (unsigned idx=0; idx<OpenWQ_wqconfig.ml_deriv_closures.size(); idx++){
        const auto& mdc = OpenWQ_wqconfig.ml_deriv_closures[idx];
        if (mdc.term=="CHEM"){
            int g=OpenWQ_wqconfig.ml_chem_group[mdc.chemi];
            for(unsigned j=0;j<num_chem;j++) if(OpenWQ_wqconfig.ml_chem_group[j]==g) OpenWQ_wqconfig.ml_chem_cl[mdc.icmp][j]=(int)idx;
        } else if (mdc.term=="SORPT"){
            int g=OpenWQ_wqconfig.ml_sorpt_group[mdc.chemi];
            for(unsigned j=0;j<num_chem;j++) if(OpenWQ_wqconfig.ml_sorpt_group[j]==g) OpenWQ_wqconfig.ml_sorpt_cl[mdc.icmp][j]=(int)idx;
        } else { // SS
            OpenWQ_wqconfig.ml_ss_cl[mdc.icmp][mdc.chemi]=(int)idx;
        }
    }

    // ---- Case-2 diagnostics: one stats accumulator per closure ----
    // (labels + dials emitted verbatim to ml_closure_diagnostics.json so the
    // report is fully self-contained; accumulators start at zero.)
    OpenWQ_wqconfig.ml_closure_stats.assign(
        OpenWQ_wqconfig.ml_deriv_closures.size(),
        OpenWQ_wqconfig::MLClosureStats());
    for (unsigned idx=0; idx<OpenWQ_wqconfig.ml_deriv_closures.size(); idx++){
        const auto& mdc = OpenWQ_wqconfig.ml_deriv_closures[idx];
        auto& st = OpenWQ_wqconfig.ml_closure_stats[idx];
        st.species        = (mdc.chemi < num_chem) ? species[mdc.chemi] : std::string("?");
        st.compartment    = OpenWQ_hostModelconfig.get_HydroComp_name_at((int)mdc.icmp);
        st.term           = mdc.term;
        st.alpha          = mdc.net.alpha;
        st.max_correction = mdc.net.max_correction;
    }

    _m = "<OpenWQ> Hybrid ML: "
        + std::to_string(OpenWQ_wqconfig.ml_deriv_closures.size())
        + " derivative closure(s) active (per-species; CHEM/SORPT group-scaled -> mass-conserving).";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, _m, true, true);
}