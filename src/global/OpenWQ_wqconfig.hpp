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
#include <ctime>
#include <memory>
#include "exprtk.hpp"
#include "OpenWQ_param.hpp"    // scalar-or-spatial model parameter (Phase 0 / ML foundation)
#include "OpenWQ_ML.hpp"       // OpenWQ_ML_closure (Layer 2 learned flux closure)
#include <string>
#include <unordered_map>
#include <sys/stat.h>
#include "PhreeqcRM.h"
#ifdef _OPENMP
#include <omp.h>
#endif



// Forward declaration
class OpenWQ_hostModelconfig;

/* #################################################
 General information for openWQ about the project
################################################# */
class OpenWQ_wqconfig
{

    // #################################################
    // Compiling and re-structuring of input data for quicker access during runtime
    // Source and source forcing (SinkSource_FORC)
    // AND
    // External fluxes (ExtFlux_FORC_jsonAscii)
    private:

        // Master file location
        std::string OpenWQ_masterjson;


        size_t num_coldata_jsonAscii = 20;
        // num_coldata_jsonAscii is, the moment, equal to 20
        // 0 - chemical
        // 1 - compartment id (from HydroComp) / external flux id (from HydroExtFlux)
        // 2 - source(=0) or sink(=1)
        // 3 - YYYY
        // 4 - MM
        // 5 - DD
        // 6 - HH
        // 7 - MIN
        // 8 - SEC
        // 9 - ix
        // 10 - iy
        // 11 - iz
        // 12 - value (already converted to mg/l (concentration) or g(mass))
        // 13 - load scheme (0) not applicable, (1) discrete or (2) continuous
        // 14, 15 ,16, 17, 18, 19 - flag to deal with "ALL" entries in YYYY, MM, DD, HH, MIN, SEC
            // if there are no "all"s, then it's to use one time only and 
            // and it is set to -1, which after use becomes -2 for not use again
            // otherwise, it gets updated everytime the load is added
            // and it provides the time increment for the next load
        size_t num_coldata_h5 = 4;
        // 1 - ix
        // 2 - iy
        // 3 - iz
        // 4 - value (already converted to mg/l (concentration) or g(mass))


        // JSON classes or output class
        std::string LogFile_name = "Log_OpenWQ.txt";       // name of log file
        std::string LogFile_name_fullpath;

        // Thread Settings
        int num_threads_system;         // number of threads in the system
        int num_threads_requested;      // number of threads requested by user
        int num_threads_default = 4;    // if requested num threads is invalid, defaults to this value 

        // Time settings
        double time_previous;       // previous time (in seconds) for calculation of dynamic
                                    // timesteps that critical for calculation of transformation rates
        double timestep_out;             // time step (in seconds)
        std::string timestep_out_unit;  // time step unit
        double nexttime_out = 0.0f;     // iteractive next printing time (in seconds)
        time_t simtime = 0;             // current simulation time (seconds since 00:00 hours, Jan 1, 1970 UTC)


        // Sink and Source AND External fluxes
        std::string h5EWF_interpMethod;     // interpolation method for h5 EWF 
        int allSS_flag = -1;                // number to replace in SinkSource_FORC to denote "all"
        bool tstep1_flag = true;            // flag to note that it's the first time step, so need to exclude loads prior to that

        // output format
        unsigned long output_type;      // 1) CSV = 0, 2) HDF5 = 1
        // Output folder
        std::string output_dir;
        void check_mkdir(); // private method used by set_output_dir and set_output_type
        std::tuple<
                std::string,                // output units as provided by the user
                double,                     // numerator multiplier (determined by Convert_Units)
                double,                     // denominator multiplier (determined by Convert_Units)
                bool                        // flag if requested concentration (needs to be divided by water volume)
                > output_units;             // Tuple with info about output units

    public:

        //###############################################
        // Hybrid physics-ML for SOURCE/SINK loading (extwatflux_ss)
        //###############################################
        // LAYER 1: optional per-cell multiplier on the load (global 1.0 by
        // default -> byte-identical; a regionalized {DEFAULT,CELLS} map makes it
        // spatial). (The Layer-2 SS closure lives with the solver derivative
        // closures, applied on dm_ss.)
        OpenWQ_param       ss_scale = OpenWQ_param(1.0);

        //###############################################
        // Hybrid physics-ML LAYER 2 — per-species DERIVATIVE closures
        //###############################################
        // A bounded, learned correction applied in the SOLVER to a species' net
        // tendency for ONE process term: CHEM (biogeochemistry), SORPT
        // (sorption) or SS (source/sink loads). It replaces the old per-reaction
        // /per-module closures with a single solver hook. Conservation is
        // guaranteed for CHEM/SORPT by scaling the whole COUPLING GROUP (species
        // linked by reactions / a dissolved-sorbed pair) by one factor, so every
        // between-species transfer stays balanced. SS is external forcing.
        // Empty ml_deriv_closures -> factors default 1.0 -> byte-identical.

        // Raw config as read from OPENWQ_INPUT > ML_CLOSURES (names, resolved to
        // indices at first solve by OpenWQ_compute::Prepare_MLClosures).
        struct MLClosureConfig {
            std::string compartment;   // compartment name, or "ALL"
            std::string species;       // target (observed) species name
            std::string term;          // "CHEM" | "SORPT" | "SS"
            OpenWQ_ML_closure net;     // the bounded factor(state)
        };
        std::vector<MLClosureConfig> ml_closure_cfg;

        // Resolved closures (indices into the model's compartment/species space).
        struct MLDerivClosure {
            unsigned int icmp;         // target compartment index
            unsigned int chemi;        // driving (observed) species index
            std::string  term;         // "CHEM" | "SORPT" | "SS"
            OpenWQ_ML_closure net;
        };
        std::vector<MLDerivClosure> ml_deriv_closures;

        // Global conservation groups (connected components): [chemi] -> group id.
        std::vector<int> ml_chem_group;     // reaction-coupled species
        std::vector<int> ml_sorpt_group;    // dissolved<->sorbed pairs
        // Per [icmp][chemi] -> index into ml_deriv_closures (or -1 = none). A
        // group member points at the closure driving its group.
        std::vector<std::vector<int>> ml_chem_cl;
        std::vector<std::vector<int>> ml_sorpt_cl;
        std::vector<std::vector<int>> ml_ss_cl;
        bool ml_closures_ready = false;     // Prepare_MLClosures ran once

        //###############################################
        // Hybrid physics-ML LAYER 2 — closure DIAGNOSTICS ("how much the ML
        // pulled"). Accumulated in the solver during the run and written to
        // <output_dir>/ml_closure_diagnostics.json on print steps, so the
        // calibration results report can quantify the correction the closure
        // applied WITHOUT re-deriving anything on the Python side.
        //   * factor = tendency_ML / tendency_physics (dimensionless);
        //     factor = 1 -> pure physics, |factor-1| = the "pull".
        //   * w = the term's (unscaled) tendency, accumulated exactly as openWQ
        //     accumulates its own `_out` flux, so sum_dev_w / sum_absdev_w are
        //     the NET / GROSS mass the closure moved in the SAME units as the
        //     model's reported flux (× output_units_numerator -> display mass).
        //   * flux-weighted sums weight the pull by |w| so near-zero-flux cells
        //     don't dominate.
        // Solver-agnostic: Forward Euler accumulates exactly per cell·step;
        // CVode evaluates the factor once per step at the step's final state.
        // Empty when no closures are active -> nothing written, no overhead.
        //###############################################
        struct MLClosureStats {
            // labels (filled by Prepare_MLClosures; emitted for the report)
            std::string species, compartment, term;
            double alpha = 0.0;             // the closure dial (0 = physics)
            double max_correction = 0.0;    // |alpha*g| clamp bound
            // accumulators
            long long n = 0;                // factor evaluations (cell·step)
            long long n_clamp = 0;          // times |factor-1| hit max_correction
            double sum_factor = 0.0;        // Σ factor
            double sum_absdev = 0.0;        // Σ |factor-1|            (mean pull)
            double max_absdev = 0.0;        // peak |factor-1|         (peak pull)
            double min_factor = 0.0;        // running min factor (seeded at n==1)
            double max_factor = 0.0;        // running max factor
            double sum_abs_w = 0.0;         // Σ |w|          (gross physics flux)
            double sum_absdev_w = 0.0;      // Σ |w·(factor-1)|   (GROSS mass moved)
            double sum_dev_w = 0.0;         // Σ  w·(factor-1)    (NET   mass moved)
            // Add one factor sample for a term whose (unscaled) tendency is `w`.
            // Thread-UNSAFE: callers accumulate into a thread-local and merge()
            // under a critical section (Forward Euler), or call serially (CVode).
            inline void add(double factor, double w){
                const double dev  = factor - 1.0;
                const double adev = dev < 0.0 ? -dev : dev;
                const double aw   = w   < 0.0 ? -w   : w;
                ++n;
                sum_factor   += factor;
                sum_absdev   += adev;
                sum_abs_w    += aw;
                sum_absdev_w += aw * adev;      // = |w·(factor-1)|
                sum_dev_w    += w  * dev;
                if (adev > max_absdev) max_absdev = adev;
                if (n == 1) { min_factor = factor; max_factor = factor; }
                else { if (factor < min_factor) min_factor = factor;
                       if (factor > max_factor) max_factor = factor; }
                if (adev >= max_correction - 1e-12) ++n_clamp;
            }
            // Fold a thread-local accumulator into this (shared) one.
            inline void merge(const MLClosureStats& a){
                if (a.n == 0) return;
                const bool first = (n == 0);
                n += a.n; n_clamp += a.n_clamp;
                sum_factor += a.sum_factor; sum_absdev += a.sum_absdev;
                sum_abs_w += a.sum_abs_w; sum_absdev_w += a.sum_absdev_w;
                sum_dev_w += a.sum_dev_w;
                if (a.max_absdev > max_absdev) max_absdev = a.max_absdev;
                if (first) { min_factor = a.min_factor; max_factor = a.max_factor; }
                else { if (a.min_factor < min_factor) min_factor = a.min_factor;
                       if (a.max_factor > max_factor) max_factor = a.max_factor; }
            }
        };
        std::vector<MLClosureStats> ml_closure_stats;   // 1:1 with ml_deriv_closures

        //###############################################
        // Methods
        //###############################################

        /***********************************************
        * OpenWQ_masterjson
        ************************************************/
        std::string get_OpenWQ_masterjson();
        void set_OpenWQ_masterjson(std::string OpenWQ_masterjson);

        /***********************************************
        * LogFile // TODO: Check if output class is a better place for this
        ************************************************/
        std::string get_LogFile_name();
        std::string get_LogFile_name_fullpath();
        void create_LogFile_name_fullpath(std::string output_format, std::string output_folder);

        /***********************************************
        * Methods for Threads
        ************************************************/
        int get_num_threads_default();
        void set_num_threads_system(int num_threads_system);
        int get_num_threads_system();
        void set_num_threads_requested(int num_threads_requested);
        int get_num_threads_requested();
        void set_threads_requested_to_system();
        void set_threads_requested_to_default();
        bool is_num_threads_requested_valid();
        std::string get_num_threads_warning();
        std::string get_num_threads_info();

        /***********************************************
        * Time Methods
        ************************************************/
        time_t get_simtime() { return simtime; }
        void set_simtime(time_t simtime_in) { simtime = simtime_in; }
        double get_time_previous();
        void set_time_previous(double time_previous);
        double get_nexttime_out();
        void set_nexttime_out(double nexttime_out);
        // update new nexttime_out with timstesp_out
        void update_nexttime_out();
        void set_timestep_out(double timestep_out);
        double get_timestep_out();
        void set_timestep_out_unit(std::string timestep_out_unit);
        std::string get_timestep_out_unit();
        void convert_units_timestep_out(std::vector<double> unit_multiplers);

        void set_h5EWF_interpMethod(std::string h5EWF_interpMethod);
        bool is_h5EWF_interpMethod(std::string h5EWF_interpMethod);
        std::string h5EWF_interp_warning_msg();

        void set_tstep1_flag_false();
        bool is_tstep1();

        int get_allSS_flag();

        /***********************************************
        * Output Format Methods //TODO: Check if output class is a better place for this
        ************************************************/
        bool is_output_type_csv();
        bool is_output_type_hdf5();
        void set_output_type_csv();
        void set_output_type_hdf5();
        std::string get_output_dir();
        void set_output_dir(std::string output_dir);
        
        // Output Units
        std::string get_output_units();
        double get_output_units_numerator();
        double get_output_units_denominator();
        bool is_concentration_requested();
        void set_output_units(std::string output_units_name);
        void set_output_units_numerator(double numerator);
        void set_output_units_denominator(double denominator);
        void set_output_units_concentration(bool concentration_requested);

        // Since the unix time epoch is 1970, which is used as a reference for timegm,
        // the seconds become negative for years below 1970, 
        // which will mess up time management.
        // Thus, the number of seconds since 00:00 1 Jan 1970 GMT, 
        // which is 2,208,988,800, is added 
        // (which is saved in OpenWQ_vars.secSinceUnixTimeEpoch).
        const unsigned long long secFrom1900toUnixTimeEpoch1970 = 2208988800;

        std::unordered_map<std::string, hid_t> files;

        // [REDESIGN] Open handles of the per-file time-extensible datasets,
        // kept open for the whole run (keyed by filename, same as `files`) so
        // each timestep just extends+writes — the active chunk stays in the
        // chunk cache and is flushed once per full chunk, instead of a
        // read-modify-write of the whole chunk on every step. Closed in the
        // destructor before the files.
        std::unordered_map<std::string, hid_t> conc_dsets;  // /concentrations
        std::unordered_map<std::string, hid_t> time_dsets;  // /timestamps

        // TODO: Below needs to be moved to private

         // ##########################
        // 3) Solver
        // General info 
        // Native Forward Euler and not native SUNDIALS
        std::string SOLVER_module;  // Get module name

        /***********************************************
        * Sink and Source AND External fluxes
        ************************************************/
        std::unique_ptr<            
            arma::Mat<double>
            > SinkSource_FORC;              // SS
        std::unique_ptr<
            arma::Mat<double>
            > ExtFlux_FORC_jsonAscii;       // External fluxes (JSON and ASCII)

        // OPTIMIZED (perf): flat row-major staging buffers used to build the
        // FORC matrices in a SINGLE allocation instead of growing them one row
        // at a time with arma::insert_rows (which reallocates+copies the whole
        // matrix on every call -> O(n^2) for large sink-source files).
        // Filled by AppendRow_SS_EWF_FORC_jsonAscii(), drained by
        // Finalize_FORC_jsonAscii() at the end of Set_EWFandSS_driver().
        std::vector<double> SinkSource_FORC_buffer;
        std::vector<double> ExtFlux_FORC_jsonAscii_buffer;
        unsigned int SinkSource_FORC_buffer_ncol = 0;        // columns per staged row
        unsigned int ExtFlux_FORC_jsonAscii_buffer_ncol = 0; // columns per staged row

        // OPTIMIZED (perf): per-row cached "next apply" time (time_t) so the
        // per-timestep SS/EWF apply loop can skip not-yet-due rows without
        // recomputing timegm() for every row at every timestep. Built on the
        // first timestep in CheckApply_EWFandSS_jsonAscii().
        std::vector<time_t> SinkSource_FORC_nextTime;
        std::vector<time_t> ExtFlux_FORC_jsonAscii_nextTime;
        std::unique_ptr<
            std::vector<       
            std::vector<time_t>
            >> ExtFlux_FORC_HDF5vec_time;   // External fluxes HDF5 vector (timestamps as time_t)
        std::unique_ptr<                    // EWF compartment id
            std::vector<unsigned int>
            > ExtFlux_FORC_HDF5vec_ewfCompID;
        std::unique_ptr<                    // per-request: BGC chem id for each
            std::vector<                    // *dense* h5 chem slot loaded (the
            std::vector<int>                // mother model may export only a
            >> ExtFlux_FORC_HDF5vec_chemID; // subset of species, so dense != BGC)
        std::unique_ptr<            
            arma::Cube<double>
            > ExtFlux_FORC_data_tStep;      // External fluxes HDF5 vector (one timestep)
        std::unique_ptr<
            std::vector<                    // JSON-h5-EWF request (blocks)   
            std::vector<                    // Chemical species
            std::vector<                    // Time steps
            arma::Cube<double>
            >>>> ExtFlux_FORC_HDF5vec_data;   // External fluxes HDF5 vector (data)

            
        bool debug_mode = false;        // set to true if debug mode is requested

        // chemicals, compartments and cells/elements to export
        std::vector<int> chem2print;
        std::vector<int> compt2print;
        std::vector<int> fluxconc2print;   // flux-concentration exports selected in OUTPUT
        // Cells to print per selected flux export (key = flux-export index).
        // Each mat is (num_cells x 3) of 0-based (ix,iy,iz) in the SOURCE
        // compartment (the cells the flux leaves from).
        std::unordered_map<int, arma::mat> fluxconc_cells2print;
        std::vector<bool> cells2print_bool;
        std::vector<arma::mat> cells2print_vec;
        // No water concentration (as a marker/flag)
        // TODO: Getter Setter
        int noWaterConc = -9999; // setting a default value

        // Flag for printing coordinates once

        bool print_oneStep = true;

        // Error message flags
        bool readSet_print_errmsg = true;
        bool BGC_Transform_print_errmsg = true;
        bool invalid_bgc_entry_errmsg = true;

        // ########################################
        // PERFORMANCE: Cached flags to avoid repeated string comparisons at runtime
        // Call cache_runtime_flags() after all modules are configured
        // ########################################
        bool is_native_bgc_flex = false;    // CH_model->BGC_module == "NATIVE_BGC_FLEX"
        bool is_TD_enabled = false;         // TD_model->TD_module != "NONE"
        bool is_LE_enabled = false;         // LE_model->LE_module != "NONE"
        bool is_TS_enabled = false;         // TS_model->TS_module != "NONE"
        bool is_TD_advdisp = false;         // TD_model->TD_module == "OPENWQ_NATIVE_TD_ADVDISP"
        bool is_TD_adv = false;             // TD_model->TD_module == "NATIVE_TD_ADV"
        unsigned int cached_num_mobile_species = 0;
        const std::vector<unsigned int>* cached_mobile_species_ptr = nullptr;
        unsigned int cached_num_chem = 0;
        const std::vector<std::string>* cached_chem_species_list_ptr = nullptr;

        // Pre-built BGC lookup: cycling framework name -> vector of indices into BGCexpressions_info
        std::unordered_map<std::string, std::vector<unsigned int>> bgc_cycle_to_transf_indices;

        // Cached solver flags
        bool is_solver_sundials = false;
        bool is_solver_forward_euler = false;

        void cache_runtime_flags();
        void build_bgc_lookup();
        void build_thread_local_expressions(
            OpenWQ_hostModelconfig& hostModelconfig);
        
        // ########################################
        // MODULES
        // ########################################

        // 1) Transport / Erosion (TE)
        class TD_model_; 
        TD_model_* TD_model;

        // 2) Biogeochemistry
        class CH_model_; 
        CH_model_* CH_model;

        // 4) Lateral Exchange (LE_model)
        class LE_model_; 
        LE_model_* LE_model;

        // 4) Transport Sediments (ST)
        class TS_model_; 
        TS_model_* TS_model;
        
        // 5) Sorption isotherm (SI_model)
        class SI_model_; 
        SI_model_* SI_model;

        OpenWQ_wqconfig();
        ~OpenWQ_wqconfig();

};


// ##############################
// Define INNER MODEL CLASSES
// ##############################

class OpenWQ_wqconfig::TD_model_{

    public:
    std::string TD_module;

    // constructor/destructor
    TD_model_(); 
    ~TD_model_();

    // ######################
    // classes for each model

    class NativeAdv_;
    NativeAdv_* NativeAdv;

    class NativeAdvDisp_;
    NativeAdvDisp_* NativeAdvDisp;

};

class OpenWQ_wqconfig::CH_model_{

    public:
    std::string BGC_module;

    // constructor/destructor
    CH_model_(); 
    ~CH_model_();

    // ######################
    // classes for each model

    class NativeFlex_;
    NativeFlex_* NativeFlex;

    class PHREEQC_;
    PHREEQC_* PHREEQC;

    // add more classes for other models here
    // e.g., PHREEQC

};

class OpenWQ_wqconfig::LE_model_{

    public:
    std::string LE_module;

    // constructor/destructor
    LE_model_(); 
    ~LE_model_();

    // ######################
    // classes for each model

    class BoundMix_;
    BoundMix_* BoundMix;

    // add more classes for other models here


};

class OpenWQ_wqconfig::TS_model_{

    public: 
    std::string TS_module;
    std::string ErodTranspCmpt;
    std::string SedCmpt;

    // constructor/destructor
    TS_model_(); 
    ~TS_model_();

    // ######################
    // classes for each model

    class HypeHVB_;
    HypeHVB_* HypeHVB;

    class HypeMMF_;
    HypeMMF_* HypeMMF;
    
};

class OpenWQ_wqconfig::SI_model_{

    public:
    std::string SI_module;

    // constructor/destructor
    SI_model_(); 
    ~SI_model_();

     // ######################
    // classes for each model

    class LANGMUIR_;
    LANGMUIR_* LANGMUIR;

    class FREUNDLICH_;
    FREUNDLICH_* FREUNDLICH;

    // ######################
    // Species-pair accessors (sorption scheme: dissolved => sorbed partner,
    // e.g. PO4-P => PP). Used by model_TS to co-transport the sorbed partner
    // species with the sediment. Defined inline at the end of this header
    // (after the inner classes are complete).
    unsigned int get_num_sorption_pairs();
    int get_sorbed_species_index_at(unsigned int pair_i);    // global chemi; -1 if unresolved
    int get_dissolved_species_index_at(unsigned int pair_i); // global chemi of the dissolved form

};

// #######################################
// Create INNER-MOST class for each model
// Inside are Model-specific variables  
// e.g., Module LE_model -> BoundMix
// #######################################

class OpenWQ_wqconfig::TD_model_::NativeAdv_{

    public:

};

class OpenWQ_wqconfig::TD_model_::NativeAdvDisp_{

    public:

    // Dispersion coefficients [m2/s] (from JSON config)
    double dispersion_x = 0.0;  // dispersion coefficient in x direction
    double dispersion_y = 0.0;  // dispersion coefficient in y direction
    double dispersion_z = 0.0;  // dispersion coefficient in z direction

    // Characteristic length [m] between cell centers
    // Used to convert D [m2/s] to an effective dispersion rate D_eff = D_avg / L^2 [1/s]
    // where D_avg = (Dx + Dy + Dz) / 3
    double characteristic_length_m = 1.0;

    // Pre-computed effective dispersion rate [1/s] = D_avg / L^2
    // Computed once after JSON parsing to avoid repeated division at runtime
    double D_eff = 0.0;

    // LAYER 1: D_eff as a scalar-or-spatial parameter. Defaults to the scalar
    // D_eff above (global -> byte-identical); becomes a per-cell field when the
    // TD config provides a regionalized "D_EFF_1/S" map (ML Layer 1).
    OpenWQ_param D_eff_param;

};

// Module LE_model -> BoundMix
class OpenWQ_wqconfig::LE_model_::BoundMix_{

    public:

    // information vector-tuple for the boundary
    std::vector
        <std::tuple<
            unsigned int,   // input_direction_index
            unsigned int,   // input_upper_compartment_index
            unsigned int,   // input_lower_compartment_index
            double          // data_format
            >> info_vector;

        

    // LAYER 1: per-entry exchange coefficient k as a scalar-or-spatial param
    // (index-aligned with info_vector). Global by default (K_VAL a number ->
    // .at() returns the scalar -> byte-identical); becomes a per-cell field
    // when K_VAL is a regionalized {DEFAULT,CELLS} map (ML Layer 1).
    std::vector<OpenWQ_param> k_param;

    // METHODS
    int get_exchange_direction(unsigned int entry_i);
    int get_upper_compartment(unsigned int entry_i);
    int get_lower_compartment(unsigned int entry_i);
    int get_k_value(unsigned int entry_i);

};

// Module CH_model -> NativeFlex
class OpenWQ_wqconfig::CH_model_::NativeFlex_{

    public:

    unsigned int num_chem;                  //Number of chemical species  

    std::vector
        <std::string> chem_species_list;    // Chemical species list
    std::vector
        <unsigned int> mobile_species;      // index of mobile chem species

    // OpenWQ native module: NATIVE_BGC_FLEX
    
    // BGC kinetic formulas (tuple with all the info needed)
    // It includes also the formulas parsed and ready to be used
    // for each BGC cyle provided by the user
    typedef exprtk::expression<double> expression_t;
    std::vector<
        std::tuple<
            std::string,                // Biogeochemical cycle name
            std::string,                // Transformation name
            std::string,                // kinetic equation provided
            unsigned int,               // index of consumed species       
            unsigned int,               // index of produced species
            std::vector<unsigned int>   // index of chemical in transformation equation (needs to be here for loop reset)
        >> BGCexpressions_info;
    
    std::vector<
        exprtk::expression<double>                      // Expression (exprtk) parsed
        >BGCexpressions_eq;            // BGC kinetic formulas for all biogeochemical cycles
    
    std::vector<double> chemass_InTransfEq; // chemical mass involved in transformation (needs to be here for loop reset)

    // Store modified expression strings (with variable substitutions and unit multipliers)
    // for re-compilation in per-thread expression copies
    std::vector<std::string> BGCexpressions_modif_strings;

    // ########################################
    // PARALLEL: Per-thread copies of expressions and their bound data vectors
    // Each thread gets its own chemass vector, dependVar_scalar vector,
    // and compiled expressions that reference the thread-local data.
    // Thread 0 uses the original (shared) data above.
    // ########################################
    int num_omp_threads = 1;
    // Per-thread chemass vectors: [thread_id][chem_index]
    std::vector<std::vector<double>> thread_chemass_InTransfEq;
    // Per-thread dependency scalar vectors: [thread_id][dep_index]
    std::vector<std::vector<double>> thread_dependVar_scalar;
    // Per-thread compiled expressions: [thread_id][transf_index]
    std::vector<std::vector<exprtk::expression<double>>> thread_BGCexpressions_eq;
    // Flag indicating thread-local copies are ready
    bool thread_local_ready = false;

    // ########################################
    // ML / SPATIAL BGC parameters (exprtk refactor)
    // A BGC parameter is either GLOBAL (string-substituted as a literal into the
    // expression - historical, byte-identical) or SPATIAL (per-cell). Spatial
    // params are left in the expression as elements of a bound "openWQ_BGCparam"
    // vector and refreshed per cell, exactly like chemass_InTransfEq. All of the
    // members below stay empty / max_BGCparam_size==0 when NO parameter is
    // spatial, so existing configs compile to identical expression strings, the
    // openWQ_BGCparam vector is never bound, and the run is byte-identical.
    // ########################################
    // Per-expression ordered list of SPATIAL params; the list index equals the
    // k used inside that expression as openWQ_BGCparam[k].
    std::vector<std::vector<OpenWQ_param>> BGCexpr_spatial_params;
    // Max number of spatial params across all expressions (vector sizing)
    unsigned int max_BGCparam_size = 0;
    // Serial path: the current cell's spatial-parameter values (bound once)
    std::vector<double> BGCparam_InTransfEq;
    // Parallel path: per-thread copies [thread_id][k]
    std::vector<std::vector<double>> thread_BGCparam_InTransfEq;

};
// Module CH_model -> PHREEQC

class OpenWQ_wqconfig::CH_model_::PHREEQC_{

    public:

    unsigned int num_chem = 0;              //Number of chemical species

    std::vector
        <std::string> chem_species_list;    // Chemical species list
    std::vector
        <unsigned int> mobile_species;      // index of mobile chem species

    // Gram formula weights (GFW) for each component [g/mol]
    // Used for unit conversion between OpenWQ (mg/L) and PHREEQC (mol/kgw)
    // Index matches chem_species_list order
    std::vector<double> gfw;                // Gram formula weights [g/mol]

    // Unit conversion factors (computed once during setup)
    // OpenWQ uses mg/L internally, PHREEQC uses mol/kgw
    // To convert mg/L to mol/kgw: conc_mol = conc_mg_L / (1000 * gfw)
    // To convert mol/kgw to mg/L: conc_mg_L = conc_mol * 1000 * gfw
    // Note: Assumes water density ~1 kg/L, so mg/L ≈ mg/kgw

    // PHREEQC
    std::unique_ptr<PhreeqcRM> phreeqcrm;
};

// Module TS_model -> HypeHVB
class OpenWQ_wqconfig::TS_model_::HypeHVB_{

    public: 

    // General configuration (compartments, direction)
    std::tuple<
        unsigned int,   // direction_index
        unsigned int,   // ErodTranspCmpt index
        unsigned int,   // Sdiment comparment index
        unsigned int,   // input_compartment_inhibitErosion_index
        std::string>    // data_format
            info_vector;

    // #################
    // PARAMETERS
    // #################
    
    // lusepar       ! soil erosion factor (land use dependence) 
    arma::Cube<double> lusepar_entryArmaCube;  

    // soilpar       ! soil erosion factor (soil dependence)
    arma::Cube<double> soilpar_entryArmaCube;  

    // slopepar      ! slope erosion factor (exponent)         
    arma::Cube<double> slopepar_entryArmaCube;  

    // precexppar    ! erosion precipitation dependence factor (exponent)        
    arma::Cube<double> precexppar_entryArmaCube;  

    // eroindexpar   ! model parameter for scaling of erosion index          
    arma::Cube<double> eroindexpar_entryArmaCube;  

     // slope       ! basin slope
    arma::Cube<double> slope_entryArmaCube;  

    // erosion index
    arma::Cube<double> eroindex_entryArmaCube;

    // Monthly erosion factor (12 values, one per month, Jan=index 0)
    // In HYPE: erodmonth = 1.0 + monthpar(m_erodmon, current_month)
    std::vector<double> monthpar = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


    // METHODS
    int get_exchange_direction();
    int get_eroding_transp_compartment();
    int get_erosion_inhibit_compartment();
    int get_sediment_compartment();
    std::string get_data_format();

};

// Module TS_model -> HypeMMF
class OpenWQ_wqconfig::TS_model_::HypeMMF_{

    public: 

    // General configuration (compartments, direction)
    std::tuple<
        unsigned int,   // direction_index
        unsigned int,   // ErodTranspCmpt index
        unsigned int,   // Sdiment comparment index
        unsigned int,   // input_compartment_inhibitErosion_index
        std::string>    // data_format
            info_vector;

    // Pseudo day-of-year (1-365), dynamically updated from simtime at runtime
    unsigned pdayno = 1;

    // #################
    // PARAMETERS
    // #################

    // Soil coheson (kPa)  
    arma::Cube<double> cohesion_entryArmaCube;  

    // Soil erodibility (g/J)          
    arma::Cube<double> erodibility_entryArmaCube;  

    // Surface runoff erosion exponent           
    arma::Cube<double> sreroexp_entryArmaCube;  

     // Crop cover           
    arma::Cube<double> cropcover_entryArmaCube;  

    // Ground cover           
    arma::Cube<double> groundcover_entryArmaCube;  

    // Slope           
    arma::Cube<double> slope_entryArmaCube; 

    // transport factor 1      
    arma::Cube<double> trans_1_entryArmaCube;  

    // transport factor 2     
    arma::Cube<double> trans_2_entryArmaCube;  

    // #################
    // INTERMEDIATE VARIABLE
    // #################

    // Potential mobilization of sediments by rainfall
    // but it can only be mobilized if there is runoff
    arma::Cube<double> mobilisedsed_rain_potential;


    // METHODS
    int get_exchange_direction();
    int get_eroding_transp_compartment();
    int get_erosion_inhibit_compartment();
    int get_sediment_compartment();
    std::string get_data_format();

};

// Module SI_model -> Langmuir Isotherm
// q = (qmax * KL * C) / (1 + KL * C)
// where q = sorbed concentration [M/M_soil], C = dissolved concentration [M/L^3]
// Kinetic: sorption_rate = Kadsdes * (q_eq - q_current) * bulk_density * layer_thickness
class OpenWQ_wqconfig::SI_model_::LANGMUIR_{

    public:

    // Number of species with Langmuir sorption configured
    unsigned int num_species = 0;

    // Per-species parameters (vectors indexed by local species index)
    std::vector<unsigned int> species_index;     // Global chemical species index (dissolved form)
    std::vector<std::string> species_name;       // Species name (for logging)
    std::vector<int> sorbed_species_index;       // Global index of the SORBED partner species
                                                 // (species-pair scheme, e.g. PO4-P => PP)
    std::vector<std::string> sorbed_species_name;// Sorbed partner name (for logging)
    // Per-species parameters as OpenWQ_param: GLOBAL scalar by default
    // (backward-compatible), or a per-cell SPATIAL field (map / ML).
    std::vector<OpenWQ_param> qmax;              // Maximum adsorption capacity [mg/kg_soil]
    std::vector<OpenWQ_param> KL;                // Langmuir equilibrium constant [L/mg]
    std::vector<OpenWQ_param> Kadsdes;           // Kinetic adsorption/desorption rate [1/s]

    // Soil/medium properties (uniform or per-compartment in future)
    double bulk_density = 1500.0;                // Bulk density [kg/m^3]
    double layer_thickness = 1.0;                // Representative layer thickness [m]


};

// Module SI_model -> Freundlich Isotherm
// q = Kfr * C^(1/Nfr)
// where q = sorbed concentration [M/M_soil], C = dissolved concentration [M/L^3]
// Uses Newton-Raphson for equilibrium + kinetic adsorption/desorption
class OpenWQ_wqconfig::SI_model_::FREUNDLICH_{

    public:

    // Number of species with Freundlich sorption configured
    unsigned int num_species = 0;

    // Per-species parameters (vectors indexed by local species index)
    std::vector<unsigned int> species_index;     // Global chemical species index (dissolved form)
    std::vector<std::string> species_name;       // Species name (for logging)
    std::vector<int> sorbed_species_index;       // Global index of the SORBED partner species
                                                 // (species-pair scheme, e.g. PO4-P => PP)
    std::vector<std::string> sorbed_species_name;// Sorbed partner name (for logging)
    // Per-species parameters. Each is an OpenWQ_param: a GLOBAL scalar by
    // default (backward-compatible, identical for every cell) that can instead
    // hold a per-cell SPATIAL field (from an input map or ML regionalization).
    std::vector<OpenWQ_param> Kfr;               // Freundlich coefficient [mg/kg / (mg/L)^(1/Nfr)]
    std::vector<OpenWQ_param> Nfr;               // Freundlich exponent [-]
    std::vector<OpenWQ_param> Kadsdes;           // Kinetic adsorption/desorption rate [1/s]

    // Soil/medium properties
    double bulk_density = 1500.0;                // Bulk density [kg/m^3]
    double layer_thickness = 1.0;                // Representative layer thickness [m]


};

// #######################################
// SI_model_ species-pair accessors (inline; defined here because they need
// the complete LANGMUIR_/FREUNDLICH_ types). They expose the sorption pairs
// of whichever isotherm module is active, so callers (e.g. model_TS moving
// the sorbed partner species with the sediment) don't need to know which
// isotherm is configured.
// #######################################

inline unsigned int OpenWQ_wqconfig::SI_model_::get_num_sorption_pairs(){
    if (SI_module.compare("LANGMUIR") == 0)   return LANGMUIR->num_species;
    if (SI_module.compare("FREUNDLICH") == 0) return FREUNDLICH->num_species;
    return 0;
}

inline int OpenWQ_wqconfig::SI_model_::get_sorbed_species_index_at(unsigned int pair_i){
    if (SI_module.compare("LANGMUIR") == 0)
        return (pair_i < LANGMUIR->sorbed_species_index.size())
            ? LANGMUIR->sorbed_species_index[pair_i] : -1;
    if (SI_module.compare("FREUNDLICH") == 0)
        return (pair_i < FREUNDLICH->sorbed_species_index.size())
            ? FREUNDLICH->sorbed_species_index[pair_i] : -1;
    return -1;
}

inline int OpenWQ_wqconfig::SI_model_::get_dissolved_species_index_at(unsigned int pair_i){
    if (SI_module.compare("LANGMUIR") == 0)
        return (pair_i < LANGMUIR->species_index.size())
            ? (int) LANGMUIR->species_index[pair_i] : -1;
    if (SI_module.compare("FREUNDLICH") == 0)
        return (pair_i < FREUNDLICH->species_index.size())
            ? (int) FREUNDLICH->species_index[pair_i] : -1;
    return -1;
}
