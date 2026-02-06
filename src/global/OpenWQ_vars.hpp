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
#include <memory>
#include <vector>
#include <algorithm> 


/* #################################################
 Main openWQ data structures
################################################# */
class OpenWQ_vars
{
    public:

        size_t num_HydroComp;
        size_t num_EWF;

        std::unique_ptr<arma::field<arma::field<arma::Cube<double>>>> 
            chemass,                    // state-variable (chemistry)
            d_chemass_dt_chem,          // derivative (chemistry)
            d_chemass_dt_transp,        // derivative (transport)
            d_chemass_ss,               // derivative (isolated) SS
            d_chemass_ewf,              // derivative (isolated) EWF
            d_chemass_ic,               // derivative (at start) IC
            d_chemass_dt_chem_out,      // cumulative derivative for output in debug model
            d_chemass_dt_transp_out,    // cumulative derivative for output in debug model
            d_chemass_ss_out,           // cumulative derivative for output in debug model
            d_chemass_ewf_out,          // cumulative derivative for output in debug model
            ewf_conc,
            d_chemass;                   // concentration of external water fluxes

        std::unique_ptr<arma::Cube<double>>
            sedmass,                        // sediment mass in flow-transport compartment (e.g, runoff)
            d_sedmass_dt,                   // derivative of sediments in flow-transport compartment (e.g. runoff)
            d_sedmass_mobilized_dt;         // sediments mobilized in sediment comparment (e.g., soil)

        // ############################################
        // Sorbed mass for sorption isotherm tracking
        // ############################################
        std::unique_ptr<arma::field<arma::field<arma::Cube<double>>>>
            sorbed_mass;                    // sorbed mass on solid phase [g]

        // ############################################
        // Mass Balance Tracking
        // ############################################
        struct MassBalanceData {
            // Per-species totals (summed across all compartments and cells)
            std::vector<double> total_dissolved_mass;   // Current dissolved mass [g]
            std::vector<double> total_sorbed_mass;      // Current sorbed mass [g]
            std::vector<double> total_mass;             // Total = dissolved + sorbed [g]

            // Cumulative fluxes since simulation start
            std::vector<double> cumulative_sources;     // Total mass added (IC + SS + EWF) [g]
            std::vector<double> cumulative_sinks;       // Total mass removed (transport out, sinks) [g]
            std::vector<double> cumulative_out_flux;    // Mass leaving system (recipient=-1) [g]
            std::vector<double> cumulative_negative_correction; // Mass lost to negative clamping [g]

            // Initial mass at simulation start
            std::vector<double> initial_mass;           // Mass at t=0 [g]

            // Mass balance error
            std::vector<double> mass_balance_error;     // Expected - Actual [g]
            std::vector<double> mass_balance_error_pct; // Error as percentage [%]

            // Flags
            bool initialized = false;
            unsigned int num_species = 0;

            // Initialize vectors for given number of species
            void initialize(unsigned int n_species) {
                num_species = n_species;
                total_dissolved_mass.resize(n_species, 0.0);
                total_sorbed_mass.resize(n_species, 0.0);
                total_mass.resize(n_species, 0.0);
                cumulative_sources.resize(n_species, 0.0);
                cumulative_sinks.resize(n_species, 0.0);
                cumulative_out_flux.resize(n_species, 0.0);
                cumulative_negative_correction.resize(n_species, 0.0);
                initial_mass.resize(n_species, 0.0);
                mass_balance_error.resize(n_species, 0.0);
                mass_balance_error_pct.resize(n_species, 0.0);
                initialized = true;
            }

            // Reset cumulative values (but keep initial mass)
            void reset_cumulative() {
                std::fill(cumulative_sources.begin(), cumulative_sources.end(), 0.0);
                std::fill(cumulative_sinks.begin(), cumulative_sinks.end(), 0.0);
                std::fill(cumulative_out_flux.begin(), cumulative_out_flux.end(), 0.0);
                std::fill(cumulative_negative_correction.begin(), cumulative_negative_correction.end(), 0.0);
            }
        };

        MassBalanceData mass_balance;

        // Constructor
        OpenWQ_vars(size_t num_HydroComp, size_t num_EWF);


};