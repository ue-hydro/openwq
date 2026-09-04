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


#include "headerfile_CC.hpp"
#include <fstream>


namespace {
// ################################################################
// Hybrid physics-ML LAYER 2 — write the closure DIAGNOSTICS accumulated by the
// solver ("how much the ML pulled") to <output_dir>/ml_closure_diagnostics.json.
// Overwritten on each print step; the final write holds the whole-run totals.
// Fully self-contained (labels, dials, factor stats, and net/gross mass moved
// in the model's output mass units) so the calibration results report reads it
// directly — no Python-side re-derivation. No-op when no closures are active.
// ################################################################
static void _write_ml_closure_diagnostics(OpenWQ_wqconfig& OpenWQ_wqconfig){

    if (OpenWQ_wqconfig.ml_closure_stats.empty()) return;

    const std::string path = OpenWQ_wqconfig.get_output_dir()
                             + "/ml_closure_diagnostics.json";
    std::ofstream f(path, std::ios::trunc);
    if (!f.is_open()) return;

    // native mass (g) -> output mass unit; take the mass part of e.g. "mg/l".
    const double num = OpenWQ_wqconfig.get_output_units_numerator();
    std::string umass = OpenWQ_wqconfig.get_output_units();
    { const std::size_t p = umass.find('/');
      if (p != std::string::npos) umass = umass.substr(0, p); }
    const std::string solver =
          OpenWQ_wqconfig.is_solver_sundials      ? "sundials"
        : OpenWQ_wqconfig.is_solver_forward_euler ? "forward_euler" : "unknown";

    auto jstr = [](const std::string& s){
        std::string o = "\"";
        for (char c : s){ if (c=='"' || c=='\\') o += '\\'; o += c; }
        return o + "\""; };

    f.precision(10);
    f << "{\n";
    f << "  \"schema\": \"openwq_ml_closure_diagnostics/1\",\n";
    f << "  \"solver\": "     << jstr(solver) << ",\n";
    f << "  \"mass_units\": " << jstr(umass)  << ",\n";
    f << "  \"note\": \"factor = tendency_ML / tendency_physics; 1.0 = pure "
         "physics; |factor-1| = the pull. CVode: factor sampled once per step "
         "at the step's final state.\",\n";
    f << "  \"closures\": [\n";

    const std::size_t N = OpenWQ_wqconfig.ml_closure_stats.size();
    for (std::size_t k = 0; k < N; k++){
        const auto& s = OpenWQ_wqconfig.ml_closure_stats[k];
        const double nd          = s.n > 0 ? (double)s.n : 1.0;
        const double factor_mean = s.n > 0 ? s.sum_factor / nd : 1.0;
        const double pull_mean   = s.n > 0 ? s.sum_absdev / nd : 0.0;
        const double pull_fw     = s.sum_abs_w > 0.0 ? s.sum_absdev_w / s.sum_abs_w : 0.0;
        const double net_fc      = s.sum_abs_w > 0.0 ? s.sum_dev_w   / s.sum_abs_w : 0.0;
        const double clampf      = s.n > 0 ? (double)s.n_clamp / nd : 0.0;

        f << "    {\n";
        f << "      \"species\": "        << jstr(s.species)     << ",\n";
        f << "      \"compartment\": "    << jstr(s.compartment) << ",\n";
        f << "      \"term\": "           << jstr(s.term)        << ",\n";
        f << "      \"alpha\": "          << s.alpha             << ",\n";
        f << "      \"max_correction\": " << s.max_correction    << ",\n";
        f << "      \"n_samples\": "      << s.n                 << ",\n";
        f << "      \"factor_mean\": "    << factor_mean         << ",\n";
        f << "      \"factor_min\": "     << (s.n > 0 ? s.min_factor : 1.0) << ",\n";
        f << "      \"factor_max\": "     << (s.n > 0 ? s.max_factor : 1.0) << ",\n";
        f << "      \"pull_mean\": "      << pull_mean           << ",\n";
        f << "      \"pull_max\": "       << s.max_absdev        << ",\n";
        f << "      \"pull_flux_weighted\": " << pull_fw         << ",\n";
        f << "      \"net_flux_change\": "    << net_fc          << ",\n";
        f << "      \"clamp_fraction\": "     << clampf          << ",\n";
        f << "      \"physics_flux\": "       << s.sum_abs_w   * num << ",\n";
        f << "      \"net_mass_moved\": "     << s.sum_dev_w   * num << ",\n";
        f << "      \"gross_mass_moved\": "   << s.sum_absdev_w * num << "\n";
        f << "    }" << (k + 1 < N ? "," : "") << "\n";
    }

    f << "  ]\n";
    f << "}\n";
}
} // namespace


// ################################################################
// Calls all functions required inside time loop
// But AFTER space loop has been finalized
// ################################################################
void OpenWQ_couplercalls::RunTimeLoopEnd(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
    OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
    OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
    OpenWQ_utils& OpenWQ_utils,                    // utility methods/functions
    OpenWQ_readjson& OpenWQ_readjson,            // read json files
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
    OpenWQ_TD_model& OpenWQ_TD_model,      // transport modules
    OpenWQ_LE_model& OpenWQ_LE_model,           // LE_model model
    OpenWQ_CH_model& OpenWQ_CH_model,                   // biochemistry modules
    OpenWQ_SI_model& OpenWQ_SI_model,
    OpenWQ_TS_model& OpenWQ_TS_model,
    OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
    OpenWQ_compute& OpenWQ_compute,
    OpenWQ_output& OpenWQ_output,
    time_t simtime){


    /* ########################################
    Solver 
    (needs to be inside the time loop, but outside the space loop
    Only place where the state-variables are changed
    ######################################## */ 

    OpenWQ_compute.SOLVER_driver(
        OpenWQ_hostModelconfig,
        OpenWQ_wqconfig,
        OpenWQ_vars, 
        OpenWQ_json,
        OpenWQ_output,
        OpenWQ_CH_model,
        OpenWQ_TS_model,
        OpenWQ_utils);

    // ################################
    // Apply deferred concentration-based ICs
    // (for cells that had zero water volume at init)
    // Must run AFTER solver and BEFORE output so that
    // transport doesn't redistribute IC mass on the first timestep
    // ################################
    OpenWQ_initiate.applyDeferredIC_Conc(
        OpenWQ_vars,
        OpenWQ_hostModelconfig,
        OpenWQ_wqconfig,
        OpenWQ_output);

    // ########################################
    // Output Results
    // ###########################################

    // Only print if time to print -> Needs to be adapted to host model time conventions
    // Note that OpenWQ_wqconfig.timetep_out converted to seconds in OpenWQ_readjson    

    // Print/Save results
    // Return if not time to print
    if (simtime >= OpenWQ_wqconfig.get_nexttime_out()){

        // Get compartment volumes for calculation of concentration (if requested)
        // neet to convert from mm (CRHM native units) to m3 (OpenWQ native units)

        // Call main printing function
        OpenWQ_output.writeResults(
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_hostModelconfig,
            OpenWQ_wqconfig,
            OpenWQ_compute,
            simtime);  // needs to be in seconds since 00:00 hours, Jan 1, 1970 UTC

        // Hybrid physics-ML LAYER 2: refresh the closure-diagnostics JSON
        // (cumulative; the last write at run end holds the whole-run totals).
        _write_ml_closure_diagnostics(OpenWQ_wqconfig);

    }

}