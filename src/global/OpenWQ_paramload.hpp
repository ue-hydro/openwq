// Copyright 2026, Diogo Costa
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

#pragma once

#include "global/OpenWQ_param.hpp"
#include "global/OpenWQ_ML.hpp"
#include "global/OpenWQ_hostModelconfig.hpp"
#include "readjson/headerfile_nlohmann.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

// #################################################
// OpenWQ_load_param
//
// Build an OpenWQ_param from a JSON parameter entry. Spatial cubes are sized
// from the host model's compartment dimensions, which are available by
// config-read time (add_HydroComp -> initmemory run before the modules read
// their parameters).
//
//   number   -> GLOBAL scalar (historical, byte-identical)
//   object   -> SPATIAL (per-cell). Supported forms:
//       {"UNIFORM": v}                          every cell = v
//                                               (equivalence test: a spatial
//                                                field that must reproduce the
//                                                scalar-v result exactly)
//       {"DEFAULT": d, "CELLS": [[c,x,y,z,v]]}  d everywhere, listed cells = v
//                                               (an actual spatially-varying
//                                                parameter)
//
// Any other JSON shape -> GLOBAL 0.0 (defensive; caller's validity checks then
// skip it, matching how an invalid scalar parameter is handled today).
//
// This is the loader that ACTIVATES the spatial path wired into the modules:
// returning a spatial OpenWQ_param makes OpenWQ_param::is_spatial() true, so the
// module uses its per-cell branch. The ML regionalization layer (Layer 1) will
// later produce these same spatial cubes from a learned mapping.
// #################################################
// #################################################
// OpenWQ_ML_from_json — build a small MLP (OpenWQ_ML) from its exported weights.
// Weights JSON shape:
//   {"layers":[{"W":[[...]],"b":[...],"activation":"tanh"}, ...],
//    "in_transform":"log1p",                    (optional input transform, before std.)
//    "in_mean":[...],  "in_std":[...],          (optional input standardization)
//    "out_scale":[...],"out_offset":[...]}      (optional output de-normalization)
// #################################################
inline OpenWQ_ML OpenWQ_ML_from_json(const nlohmann::json& jw)
{
    OpenWQ_ML net;
    auto to_vec = [](const nlohmann::json& a){
        arma::vec v(a.size());
        for (unsigned int i=0;i<a.size();i++) v(i)=a[i].get<double>();
        return v;
    };
    if (jw.contains("in_transform") && jw["in_transform"].is_string())
        net.in_transform = jw["in_transform"].get<std::string>();
    if (jw.contains("in_mean") && jw.contains("in_std")) {
        net.in_mean = to_vec(jw["in_mean"]);
        net.in_std  = to_vec(jw["in_std"]);
    }
    if (jw.contains("out_scale"))  net.out_scale  = to_vec(jw["out_scale"]);
    if (jw.contains("out_offset")) net.out_offset = to_vec(jw["out_offset"]);
    if (jw.contains("layers")) {
        for (const auto& lj : jw["layers"]) {
            OpenWQ_ML::Layer L;
            const auto& Wj = lj["W"];
            const unsigned int nout = Wj.size();
            const unsigned int nin  = nout ? Wj[0].size() : 0;
            L.W.set_size(nout, nin);
            for (unsigned int i=0;i<nout;i++)
                for (unsigned int k=0;k<nin;k++)
                    L.W(i,k) = Wj[i][k].get<double>();
            L.b = to_vec(lj["b"]);
            L.activation = lj.value("activation", std::string("linear"));
            net.layers.push_back(std::move(L));
        }
    }
    return net;
}

// #################################################
// OpenWQ_load_closure — LAYER 2 (learned flux closure) from its config block.
//
// In the MAIN openWQ config the block is FLAT (openWQ's loader upper-cases keys
// AND values and can't handle nested objects/booleans), so the trained network
// lives in a SEPARATE file referenced by a *FILEPATH key (whose value is left
// untouched by the normalizer):
//   {"ALPHA": 0.3, "MAX_CORRECTION": 0.5, "WEIGHTS_FILEPATH": "openwq_in/g.json"}
// (`alpha == 0` -> factor() == 1.0 = exact physics; there is no on/off boolean.)
//
// Inline "weights" are still accepted for direct (non-normalized) callers such
// as unit tests. Keys are read case-insensitively (upper- from a real config,
// lower- from a hand-built json).
// #################################################
inline OpenWQ_ML_closure OpenWQ_load_closure(const nlohmann::json& j)
{
    OpenWQ_ML_closure cl;

    auto num = [&](const char* lo, const char* up, double def)->double{
        if (j.contains(up)) return j[up].get<double>();
        if (j.contains(lo)) return j[lo].get<double>();
        return def;
    };
    cl.alpha          = num("alpha", "ALPHA", 0.0);
    cl.max_correction = num("max_correction", "MAX_CORRECTION", 1.0);
    cl.enabled        = (cl.alpha != 0.0);   // alpha=0 => pure physics

    // Weights: a separate file (real configs) or inline (direct callers).
    std::string wpath;
    if (j.contains("WEIGHTS_FILEPATH")) wpath = j["WEIGHTS_FILEPATH"].get<std::string>();
    else if (j.contains("weights_filepath")) wpath = j["weights_filepath"].get<std::string>();
    if (!wpath.empty()) {
        std::ifstream wf(wpath);
        if (wf) { nlohmann::json jw; wf >> jw; cl.net = OpenWQ_ML_from_json(jw); }
    } else if (j.contains("weights")) {
        cl.net = OpenWQ_ML_from_json(j["weights"]);
    } else if (j.contains("WEIGHTS")) {
        cl.net = OpenWQ_ML_from_json(j["WEIGHTS"]);
    }
    return cl;
}

// #################################################
// OpenWQ_load_param_runtime — LAYER 1, Mode B (runtime NN inference).
// Evaluate a trained MLP per cell AT LOAD TIME to fill the spatial cubes, so no
// pre-baked value map is needed. The config carries the net + per-cell attribute
// vectors (already mapped to internal indices):
//   {"weights": {...}, "default": d,
//    "attributes": [[icmp, ix, iy, iz, a1, a2, ...], ...]}
// This is the same trained network as Mode A; Mode B just moves the evaluation
// from Python into openWQ (config stays transparent + retrainable).
// #################################################
inline OpenWQ_param OpenWQ_load_param_runtime(
    const nlohmann::json& j,
    OpenWQ_hostModelconfig& host)
{
    // Weights: a separate file (real config; its lowercase keys survive openWQ's
    // config normalizer, which the nested inline form does not) or inline (for
    // direct / non-normalized callers such as unit tests).
    OpenWQ_ML net;
    {
        std::string wpath;
        if (j.contains("WEIGHTS_FILEPATH")) wpath = j["WEIGHTS_FILEPATH"].get<std::string>();
        else if (j.contains("weights_filepath")) wpath = j["weights_filepath"].get<std::string>();
        if (!wpath.empty()) {
            std::ifstream wf(wpath);
            if (wf) { nlohmann::json jw; wf >> jw; net = OpenWQ_ML_from_json(jw); }
        } else if (j.contains("weights")) net = OpenWQ_ML_from_json(j["weights"]);
        else if (j.contains("WEIGHTS")) net = OpenWQ_ML_from_json(j["WEIGHTS"]);
    }

    const unsigned int ncmp = host.get_num_HydroComp();
    const double def = j.contains("DEFAULT") ? j["DEFAULT"].get<double>()
                     : (j.contains("default") ? j["default"].get<double>() : 0.0);

    std::vector<arma::Cube<double>> cubes;
    cubes.reserve(ncmp);
    for (unsigned int c=0;c<ncmp;c++){
        const unsigned int nx = host.get_HydroComp_num_cells_x_at(c);
        const unsigned int ny = host.get_HydroComp_num_cells_y_at(c);
        const unsigned int nz = host.get_HydroComp_num_cells_z_at(c);
        cubes.emplace_back(nx, ny, nz);
        cubes.back().fill(def);
    }

    // Per-cell attribute rows: a separate file (real config) or inline.
    nlohmann::json attrs;
    {
        std::string apath;
        if (j.contains("ATTRIBUTES_FILEPATH")) apath = j["ATTRIBUTES_FILEPATH"].get<std::string>();
        else if (j.contains("attributes_filepath")) apath = j["attributes_filepath"].get<std::string>();
        if (!apath.empty()) { std::ifstream af(apath); if (af) af >> attrs; }
        else if (j.contains("attributes")) attrs = j["attributes"];
        else if (j.contains("ATTRIBUTES")) attrs = j["ATTRIBUTES"];
    }
    unsigned int n_eval = 0, n_skip = 0;
    double vmin = 0.0, vmax = 0.0;
    if (!attrs.is_null()) {
        for (const auto& e : attrs) {
            if (!e.is_array() || e.size() < 5) continue;  // icmp,ix,iy,iz + >=1 attr
            // [icmp, ix, iy, iz, a1, a2, ...] — ix/iy/iz ONE-BASED (openWQ cell
            // convention == xyz_elements), icmp 0-based compartment index.
            const long long c1 = e[0].get<long long>();
            const long long x1 = e[1].get<long long>();
            const long long y1 = e[2].get<long long>();
            const long long z1 = e[3].get<long long>();
            // -1 = wildcard ("every compartment" / "every layer"): a unit
            // (SUMMA HRU column, mizuRoute reach) gets the value in all the
            // compartments and vertical elements it spans.
            if (c1 < -1 || x1 < 1 || y1 < -1 || y1 == 0 || z1 < -1 || z1 == 0) { n_skip++; continue; }
            arma::vec xin(e.size()-4);
            for (unsigned int k=4;k<e.size();k++) xin(k-4)=e[k].get<double>();
            const double v = net.forward_scalar(xin);
            bool hit = false;
            for (unsigned int c = 0; c < cubes.size(); c++) {
                if (c1 >= 0 && c != (unsigned int)c1) continue;
                const unsigned int x = (unsigned int)(x1 - 1);
                if (x >= cubes[c].n_rows) continue;
                const unsigned int y0 = (y1 < 0) ? 0 : (unsigned int)(y1 - 1);
                const unsigned int y9 = (y1 < 0) ? cubes[c].n_cols : y0 + 1;
                const unsigned int z0 = (z1 < 0) ? 0 : (unsigned int)(z1 - 1);
                const unsigned int z9 = (z1 < 0) ? cubes[c].n_slices : z0 + 1;
                for (unsigned int y = y0; y < y9 && y < cubes[c].n_cols; y++)
                    for (unsigned int z = z0; z < z9 && z < cubes[c].n_slices; z++)
                        { cubes[c](x,y,z) = v; hit = true; }
            }
            if (hit) {
                if (n_eval == 0) { vmin = v; vmax = v; }
                else { vmin = std::min(vmin, v); vmax = std::max(vmax, v); }
                n_eval++;
            } else {
                n_skip++;
            }
        }
    }
    // Activation trace (stdout -> the run's model_output.log): proves the
    // network was found, evaluated and applied, and to how many cells.
    std::cout << "<OpenWQ> ML_RUNTIME (Layer 1B): net "
              << (net.empty() ? "EMPTY (weights not loaded -> DEFAULT everywhere)"
                              : std::to_string(net.layers.size()) + " layer(s)")
              << ", attributes rows evaluated = " << n_eval
              << (n_skip ? (", out-of-range rows skipped = " + std::to_string(n_skip)) : std::string(""))
              << ", default = " << def
              << (n_eval ? (", field min/max = " + std::to_string(vmin) + "/" + std::to_string(vmax))
                         : std::string(""))
              << std::endl;

    OpenWQ_param p;
    p.set_spatial(std::move(cubes));
    return p;
}

inline OpenWQ_param OpenWQ_load_param(
    const nlohmann::json& jval,
    OpenWQ_hostModelconfig& host)
{
    // ---- GLOBAL scalar (default / historical) ----
    if (jval.is_number())
        return OpenWQ_param(jval.get<double>());

    // ---- SPATIAL field ----
    if (jval.is_object()) {

        // LAYER 1, Mode B: runtime NN inference (evaluate the net per cell now)
        if (jval.contains("ML_RUNTIME"))
            return OpenWQ_load_param_runtime(jval["ML_RUNTIME"], host);

        const unsigned int ncmp = host.get_num_HydroComp();

        // Baseline value for every cell
        double def = 0.0;
        if (jval.contains("UNIFORM"))      def = jval["UNIFORM"].get<double>();
        else if (jval.contains("DEFAULT")) def = jval["DEFAULT"].get<double>();

        // One cube per compartment, sized to that compartment's cells
        std::vector<arma::Cube<double>> cubes;
        cubes.reserve(ncmp);
        for (unsigned int c = 0; c < ncmp; c++) {
            const unsigned int nx = host.get_HydroComp_num_cells_x_at(c);
            const unsigned int ny = host.get_HydroComp_num_cells_y_at(c);
            const unsigned int nz = host.get_HydroComp_num_cells_z_at(c);
            cubes.emplace_back(nx, ny, nz);
            cubes.back().fill(def);
        }

        // Per-cell overrides: [compartment, ix, iy, iz, value]
        // CELLS rows are [icmp, ix, iy, iz, value] with ix/iy/iz ONE-BASED
        // (-1 in icmp / iy / iz = wildcard: every compartment / column / layer) —
        // the openWQ cell convention used by every JSON input and by the
        // xyz_elements written to the HDF5 output (ix+1), which is where the
        // calibration tools take them from. icmp is the 0-based index into the
        // host model's compartment list (produced by the tools, not by hand).
        unsigned int n_applied = 0, n_skipped = 0;
        if (jval.contains("CELLS")) {
            for (const auto& e : jval["CELLS"]) {
                if (!e.is_array() || e.size() < 5) { n_skipped++; continue; }
                const long long c1 = e[0].get<long long>();
                const long long x1 = e[1].get<long long>();
                const long long y1 = e[2].get<long long>();
                const long long z1 = e[3].get<long long>();
                // icmp -1 = every compartment; iy/iz -1 = every column/layer
                // (a SUMMA HRU column spans several compartments and soil
                // layers; a mizuRoute reach is one cell either way).
                if (c1 < -1 || x1 < 1 || y1 < -1 || y1 == 0 || z1 < -1 || z1 == 0) { n_skipped++; continue; }
                const double       v = e[4].get<double>();
                unsigned int n_hit = 0;
                for (unsigned int c = 0; c < cubes.size(); c++) {
                    if (c1 >= 0 && c != (unsigned int)c1) continue;
                    const unsigned int x = (unsigned int)(x1 - 1);
                    if (x >= cubes[c].n_rows) continue;
                    const unsigned int y0 = (y1 < 0) ? 0 : (unsigned int)(y1 - 1);
                    const unsigned int y9 = (y1 < 0) ? cubes[c].n_cols : y0 + 1;
                    const unsigned int z0 = (z1 < 0) ? 0 : (unsigned int)(z1 - 1);
                    const unsigned int z9 = (z1 < 0) ? cubes[c].n_slices : z0 + 1;
                    for (unsigned int y = y0; y < y9 && y < cubes[c].n_cols; y++)
                        for (unsigned int z = z0; z < z9 && z < cubes[c].n_slices; z++)
                            { cubes[c](x, y, z) = v; n_hit++; }
                }
                if (n_hit) n_applied += n_hit; else n_skipped++;
            }
        }

        // Trace (stdout -> model_output.log): a map that applies to 0 cells is a
        // silent no-op otherwise — say how many cells it reached.
        std::cout << "<OpenWQ> SPATIAL PARAM: " << n_applied << " cell(s) set from CELLS"
                  << (n_skipped ? (", " + std::to_string(n_skipped) + " row(s) skipped (out of range / invalid)") : std::string(""))
                  << ", default = " << def << std::endl;
        OpenWQ_param p;
        p.set_spatial(std::move(cubes));
        return p;
    }

    // ---- defensive fallback ----
    return OpenWQ_param(0.0);
}
