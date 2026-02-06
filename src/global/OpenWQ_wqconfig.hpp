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
        std::unique_ptr<
            std::vector<       
            std::vector<time_t>
            >> ExtFlux_FORC_HDF5vec_time;   // External fluxes HDF5 vector (timestamps as time_t)
        std::unique_ptr<                    // EWF compartment id
            std::vector<unsigned int>
            > ExtFlux_FORC_HDF5vec_ewfCompID;
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

};

// #######################################
// Create INNER-MOST class for each model
// Inside are Model-specific variables  
// e.g., Module LE_model -> BoundMix
// #######################################

class OpenWQ_wqconfig::TD_model_::NativeAdv_{

    public:

    // add variables here if needed
    // current not passing any variables

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
    std::vector<unsigned int> species_index;     // Global chemical species index
    std::vector<std::string> species_name;       // Species name (for logging)
    std::vector<double> qmax;                    // Maximum adsorption capacity [mg/kg_soil]
    std::vector<double> KL;                      // Langmuir equilibrium constant [L/mg]
    std::vector<double> Kadsdes;                 // Kinetic adsorption/desorption rate [1/s]

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
    std::vector<unsigned int> species_index;     // Global chemical species index
    std::vector<std::string> species_name;       // Species name (for logging)
    std::vector<double> Kfr;                     // Freundlich coefficient [mg/kg / (mg/L)^(1/Nfr)]
    std::vector<double> Nfr;                     // Freundlich exponent [-]
    std::vector<double> Kadsdes;                 // Kinetic adsorption/desorption rate [1/s]

    // Soil/medium properties
    double bulk_density = 1500.0;                // Bulk density [kg/m^3]
    double layer_thickness = 1.0;                // Representative layer thickness [m]

};
