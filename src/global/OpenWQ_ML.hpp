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
#include <string>
#include <cmath>

// #################################################
// OpenWQ_ML
//
// A small feed-forward neural network (MLP) for the hybrid physics-ML layers.
// Hand-rolled with Armadillo (matrix-vector products) - NO ONNX / Torch / other
// runtime dependency, since the networks used for regionalization and flux
// closures are tiny.
//
//   * LAYER 1, Mode B (parameter regionalization at runtime): openWQ evaluates
//     the trained mapping theta = g(attributes) itself, per cell, at init, and
//     fills the OpenWQ_param spatial cubes - so no pre-baked map file is needed.
//     (Mode A bakes the SAME network's output into {"DEFAULT","CELLS"} maps in
//     Python; Mode B just moves the evaluation into openWQ.)
//
//   * LAYER 2 (learned flux closures): the same class evaluates a bounded
//     correction g(state, attributes) applied as flux * (1 + alpha * g).
//
// Training happens in Python (autodiff); this class does INFERENCE only. The
// trained weights are exported to a small JSON (see OpenWQ_paramload.hpp for
// the json -> OpenWQ_ML builder).
// #################################################
class OpenWQ_ML {

    public:

        struct Layer {
            arma::mat   W;            // weight matrix [n_out, n_in]
            arma::vec   b;            // bias vector   [n_out]
            std::string activation;   // "tanh" | "relu" | "sigmoid" | "linear"
        };

        std::vector<Layer> layers;

        // Optional per-feature input standardization: (x - mean) / std, applied
        // before the first layer (matches how the network was trained). Empty
        // vectors = no standardization.
        arma::vec in_mean;
        arma::vec in_std;

        // Optional element-wise input transform applied BEFORE the
        // standardization: "" / "none" = identity; "log1p" = log(1 + max(x,0)).
        // Lets a raw state (e.g. a mass in g spanning many decades across
        // compartments) drive a tanh network without saturating it.
        std::string in_transform;

        // Optional output affine de-normalization: y * scale + offset (applied
        // after the last layer). Empty vectors = no de-normalization.
        arma::vec out_scale;
        arma::vec out_offset;

        bool empty() const { return layers.empty(); }

        // #############################
        // Forward pass: attribute vector for ONE cell -> output vector.
        // #############################
        arma::vec forward(const arma::vec& x_in) const {
            arma::vec h = x_in;

            if (in_transform == "log1p")
                h.transform([](double v){ return std::log1p(v > 0.0 ? v : 0.0); });

            if (!in_mean.is_empty())
                h = (h - in_mean) / in_std;

            for (const auto& L : layers) {
                h = L.W * h + L.b;
                if      (L.activation == "tanh")    h = arma::tanh(h);
                else if (L.activation == "relu")    h = arma::clamp(h, 0.0, arma::datum::inf);
                else if (L.activation == "sigmoid") h = 1.0 / (1.0 + arma::exp(-h));
                // "linear" (or anything else) -> identity
            }

            if (!out_scale.is_empty())
                h = h % out_scale + out_offset;    // element-wise

            return h;
        }

        // Convenience for single-output networks -> scalar.
        double forward_scalar(const arma::vec& x_in) const {
            arma::vec y = forward(x_in);
            return y.is_empty() ? 0.0 : y(0);
        }
};


// #################################################
// OpenWQ_ML_closure — LAYER 2 (learned flux closure).
//
// A BOUNDED multiplicative correction to a physics flux:
//     y = Phi_physics(state, theta) * factor(inputs)
//     factor = 1 + alpha * g
// where
//   * g   = the network output bounded to [-1, 1] (the network is trained to
//           end in tanh; we also clamp defensively),
//   * alpha = the user's dial. alpha = 0 recovers the governing equation
//           EXACTLY (pure physics), so the layer is a true opt-in,
//   * the total correction (alpha * g) is clamped to +/- max_correction, so the
//     learned term can NUDGE the flux but never flip its sign or replace
//     mass/energy conservation.
//
// The network is trained in Python against a DIFFERENTIABLE EMULATOR of the
// module (autodiff); this class does inference only.
// #################################################
class OpenWQ_ML_closure {

    public:

        OpenWQ_ML net;
        double    alpha = 0.0;            // the dial (0 = pure physics)
        double    max_correction = 1.0;  // |alpha*g| clamped to this
        bool      enabled = false;

        // Multiplicative correction factor for a flux. `inputs` is the feature
        // vector (state variables + static attributes) the closure was trained
        // on. Returns 1.0 (exact physics) when disabled / alpha == 0 / no net.
        double factor(const arma::vec& inputs) const {
            if (!enabled || alpha == 0.0 || net.empty())
                return 1.0;                          // exact physics
            double g = net.forward_scalar(inputs);
            if (g >  1.0) g =  1.0;                  // bound the network output
            if (g < -1.0) g = -1.0;
            double corr = alpha * g;                 // scaled by the dial
            if (corr >  max_correction) corr =  max_correction;
            if (corr < -max_correction) corr = -max_correction;
            return 1.0 + corr;                       // 1 at alpha=0 -> physics
        }
};
