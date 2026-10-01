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

#include <armadillo>
#include <vector>

// #################################################
// OpenWQ_param
//
// A single model parameter that is EITHER a global scalar OR a spatial field.
//
//   * GLOBAL mode (default): one scalar value, identical for every cell. This
//     is the historical openWQ behaviour and is what every parameter falls back
//     to when no spatial map / ML regionalization is configured. Wiring a module
//     to OpenWQ_param therefore never changes its numerical result until a
//     spatial field is actually supplied.
//
//   * SPATIAL mode: one arma::Cube<double> per hydrological compartment, giving
//     a distinct value per (icmp, ix, iy, iz). The cubes are populated either
//     from an input map (Phase 0 - spatial parameters) or, later, by the ML
//     regionalization layer (OpenWQ_ML, Layer 1). The module code is identical
//     in both cases - it just reads at(icmp,ix,iy,iz).
//
// This is the shared foundation both hybrid physics-ML layers build on:
//   Layer 1 (parameter learning) writes the SPATIAL cubes;
//   with ML switched off the parameter stays a GLOBAL scalar.
// #################################################
class OpenWQ_param {

    public:

        // Storage mode
        enum class Mode { GLOBAL, SPATIAL };

        OpenWQ_param() = default;
        explicit OpenWQ_param(double scalar_value)
            : mode_(Mode::GLOBAL), scalar_(scalar_value) {}

        // #############################
        // Configuration
        // #############################

        // Set as a global scalar (default / backward-compatible)
        void set_global(double scalar_value){
            mode_   = Mode::GLOBAL;
            scalar_ = scalar_value;
            cubes_.clear();
        }

        // Set as a spatial field: one cube per compartment (moved in)
        void set_spatial(std::vector<arma::Cube<double>>&& cubes){
            mode_  = Mode::SPATIAL;
            cubes_ = std::move(cubes);
        }

        // #############################
        // Queries
        // #############################

        bool   is_spatial() const { return mode_ == Mode::SPATIAL; }
        double scalar()     const { return scalar_; }

        // Per-cell value.
        //   GLOBAL  -> returns the scalar (negligible cost, so a module wired to
        //              OpenWQ_param is byte-identical to the scalar version until
        //              a spatial field is supplied).
        //   SPATIAL -> cube lookup for this compartment/cell.
        inline double at(unsigned int icmp,
                         unsigned int ix,
                         unsigned int iy,
                         unsigned int iz) const {
            if (mode_ == Mode::GLOBAL) return scalar_;
            return cubes_[icmp](ix, iy, iz);
        }

    private:

        Mode   mode_   = Mode::GLOBAL;
        double scalar_ = 0.0;
        std::vector<arma::Cube<double>> cubes_;   // one per compartment (SPATIAL)
};
