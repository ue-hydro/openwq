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

#include "OpenWQ_wqconfig.hpp"
#include "OpenWQ_hostModelConfig.hpp"


// Constructor
OpenWQ_wqconfig::OpenWQ_wqconfig() {
      
    // Sink and source forcing
    SinkSource_FORC = 
        std::unique_ptr<
        arma::Mat<double>>
        (new  arma::mat(0,this->num_coldata_jsonAscii));

    // External fluxes forcing (JSON or ASCII datatypes)
    ExtFlux_FORC_jsonAscii = 
        std::unique_ptr<
        arma::Mat<double>>
        (new  arma::mat(0,this->num_coldata_jsonAscii));
    
    // External fluxes forcing (HDF5) 
    // Storing timestamps as time_t
    ExtFlux_FORC_HDF5vec_time =
        std::unique_ptr<        // EWF-h5 json block/request
        std::vector<            
        std::vector<time_t>>>
        (new  std::vector<std::vector<time_t>>);
    // Storing external_flux id
    ExtFlux_FORC_HDF5vec_ewfCompID =
        std::unique_ptr<
        std::vector<unsigned int>>
        (new  std::vector<unsigned int>);
    // Saving 1 timestep
    ExtFlux_FORC_data_tStep = 
        std::unique_ptr<
        arma::Cube<double>>
        (new  arma::cube);
    // Storing all timesteps
    ExtFlux_FORC_HDF5vec_data = 
        std::unique_ptr<
        std::vector<           // EWF-h5 json block/request
        std::vector<           // ChemSpecies 
        std::vector<           // timestamps
        arma::Cube<double>>>>>
        (new  std::vector<std::vector<std::vector<arma::cube>>>);

    // Modules
    TD_model = new TD_model_();
    TS_model = new TS_model_();
    CH_model = new CH_model_();
    LE_model = new LE_model_();
    SI_model = new SI_model_();

}

OpenWQ_wqconfig::TD_model_::TD_model_() {
    NativeAdv = new NativeAdv_();
    NativeAdvDisp = new NativeAdvDisp_();
}
OpenWQ_wqconfig::TD_model_::~TD_model_() {
    delete NativeAdvDisp;
}

OpenWQ_wqconfig::LE_model_::LE_model_() {
    BoundMix = new BoundMix_();
}
OpenWQ_wqconfig::LE_model_::~LE_model_() {
    delete BoundMix;
}

OpenWQ_wqconfig::CH_model_::CH_model_() {
    NativeFlex = new NativeFlex_();
    PHREEQC = new PHREEQC_();
}
OpenWQ_wqconfig::CH_model_::~CH_model_() {
    delete NativeFlex;
    delete PHREEQC;
}

OpenWQ_wqconfig::TS_model_::TS_model_() {
    HypeHVB = new HypeHVB_();
    HypeMMF = new HypeMMF_();
}
OpenWQ_wqconfig::TS_model_::~TS_model_() {   
    delete HypeHVB;
    delete HypeMMF;
}

OpenWQ_wqconfig::SI_model_::SI_model_() {
    LANGMUIR = new LANGMUIR_();
    FREUNDLICH = new FREUNDLICH_();
}
OpenWQ_wqconfig::SI_model_::~SI_model_() {   
    delete LANGMUIR;
    delete FREUNDLICH;
}

// Destructor
OpenWQ_wqconfig::~OpenWQ_wqconfig() {

    for (auto x: files) {
        H5Fclose(x.second);
    }

    delete TD_model;
    delete TS_model;
    delete CH_model;
    delete LE_model;
    delete SI_model;

}

//###############################################
// Methods 
//###############################################


/***********************************************
 * OpenWQ_masterjson
************************************************/
std::string OpenWQ_wqconfig::get_OpenWQ_masterjson() 
{
    return this->OpenWQ_masterjson;
}
void OpenWQ_wqconfig::set_OpenWQ_masterjson(std::string OpenWQ_masterjson) 
{
    this->OpenWQ_masterjson = OpenWQ_masterjson;  
}

/***********************************************
 * LogFile
************************************************/
std::string OpenWQ_wqconfig::get_LogFile_name() 
{
    return this->LogFile_name;
}

std::string OpenWQ_wqconfig::get_LogFile_name_fullpath() 
{
    return this->LogFile_name_fullpath;
}
void OpenWQ_wqconfig::create_LogFile_name_fullpath(
    std::string output_format, std::string output_folder) 
{
    this->LogFile_name_fullpath = this->LogFile_name;
    this->LogFile_name_fullpath.insert(0,"/");
    this->LogFile_name_fullpath.insert(0, output_format);
    this->LogFile_name_fullpath.insert(0,"/");
    this->LogFile_name_fullpath.insert(0, output_folder);

    std::ofstream outfile(this->LogFile_name_fullpath);
}

/***********************************************
 * Methods for Threads
************************************************/
int OpenWQ_wqconfig::get_num_threads_default() {
    return this->num_threads_default;
}

void OpenWQ_wqconfig::set_num_threads_system(int num_threads_system) {
    this->num_threads_system = num_threads_system;
}

void OpenWQ_wqconfig::set_num_threads_requested(int num_threads_requested) {
    this->num_threads_requested = num_threads_requested;
}
int OpenWQ_wqconfig::get_num_threads_requested() {
    return this->num_threads_requested;
}

void OpenWQ_wqconfig::set_threads_requested_to_system() {
  this->num_threads_requested = this->num_threads_system;
}

void OpenWQ_wqconfig::set_threads_requested_to_default() {
  this->num_threads_requested = this->num_threads_default;
}

bool OpenWQ_wqconfig::is_num_threads_requested_valid() {
  return this->num_threads_system >= this->num_threads_requested;
}

std::string OpenWQ_wqconfig::get_num_threads_warning() {
  return "<OpenWQ> WARNING: Number of threads in the system (" 
      + std::to_string(this->num_threads_system)
      + ") is lower than requested ("
      + std::to_string(this->num_threads_requested)
      + "). All system threads available have been engaged.";
}

std::string OpenWQ_wqconfig::get_num_threads_info() {
  return "<OpenWQ> INFO: Number of threads requested and used: " 
    + std::to_string(this->num_threads_requested) + ".";;
}

/***********************************************
 * Time Methods
************************************************/
double OpenWQ_wqconfig::get_time_previous() {
    return this->time_previous;
}
void OpenWQ_wqconfig::set_time_previous(double time_previous) {
    this->time_previous = time_previous;  
}

double OpenWQ_wqconfig::get_nexttime_out() {
    return this->nexttime_out;
}
void OpenWQ_wqconfig::set_nexttime_out(double nexttime_out) {
    this->nexttime_out = nexttime_out;  
}
void OpenWQ_wqconfig::update_nexttime_out() {
    this->nexttime_out += this->timestep_out;
}
void OpenWQ_wqconfig::set_timestep_out(double timestep_out) {
    this->timestep_out = timestep_out;
}
double OpenWQ_wqconfig::get_timestep_out() {
    return this->timestep_out;
}
void OpenWQ_wqconfig::set_timestep_out_unit(std::string timestep_out_unit) {
    this->timestep_out_unit = timestep_out_unit;
}
std::string OpenWQ_wqconfig::get_timestep_out_unit() {
    return this->timestep_out_unit;
}
// TODO: Reduntant but necessary for now
void OpenWQ_wqconfig::convert_units_timestep_out(std::vector<double> unit_multiplers){    
    this->timestep_out = this->timestep_out * unit_multiplers[0] / unit_multiplers[1];
}

/***********************************************
* Sink and Source AND External fluxes Methods
************************************************/

void OpenWQ_wqconfig::set_h5EWF_interpMethod(std::string h5EWF_interpMethod) {
    this->h5EWF_interpMethod = h5EWF_interpMethod;
}

bool OpenWQ_wqconfig::is_h5EWF_interpMethod(std::string h5EWF_interpMethod) {
    return this->h5EWF_interpMethod.compare(h5EWF_interpMethod) == 0;
}

std::string OpenWQ_wqconfig::h5EWF_interp_warning_msg() {
    return  "<OpenWQ> WARNING: EWF load using HDF5 file 'INTERPOLATION' unkown:"
            + this->h5EWF_interpMethod
            + ". 'INTERPOLATION' defaulted to 'STEP'";
}

void OpenWQ_wqconfig::set_tstep1_flag_false() 
{
    this->tstep1_flag = false;
}

bool OpenWQ_wqconfig::is_tstep1() 
{
    return this->tstep1_flag;
}

int OpenWQ_wqconfig::get_allSS_flag() 
{
    return this->allSS_flag;
}

/***********************************************
* Output Format Methods //TODO: Check if output class is a better place for this
************************************************/
bool OpenWQ_wqconfig::is_output_type_csv() 
{
    return this->output_type == 0;
}
bool OpenWQ_wqconfig::is_output_type_hdf5() 
{
    return this->output_type == 1;
}

void OpenWQ_wqconfig::set_output_type_csv()
{
    this->output_type = 0;
    this->output_dir.append("/CSV");
    check_mkdir();
}

void OpenWQ_wqconfig::set_output_type_hdf5()
{
    this->output_type = 1;
    this->output_dir.append("/HDF5");
    check_mkdir();
}

std::string OpenWQ_wqconfig::get_output_dir()
{
    return this->output_dir;
}

void OpenWQ_wqconfig::set_output_dir(std::string output_dir)
{
    this->output_dir = output_dir;
    check_mkdir();
}

// Check if directory exists and create it
void OpenWQ_wqconfig::check_mkdir()
{
    struct stat st = {0};
    
    // convert to *char
    const char *cstr = this->output_dir.c_str();

    // mkdir
    if (stat(cstr, &st) == -1) {
        mkdir(cstr, 0700);
    }
}

// Output Units
std::string OpenWQ_wqconfig::get_output_units() 
{return std::get<0>(this->output_units);}

double OpenWQ_wqconfig::get_output_units_numerator() 
{return std::get<1>(this->output_units);}

double OpenWQ_wqconfig::get_output_units_denominator() 
{return std::get<2>(this->output_units);}

bool OpenWQ_wqconfig::is_concentration_requested()
{return std::get<3>(this->output_units);}

void OpenWQ_wqconfig::set_output_units(std::string output_units)
{std::get<0>(this->output_units) = output_units;}

void OpenWQ_wqconfig::set_output_units_numerator(double output_units_numerator)
{std::get<1>(this->output_units) = output_units_numerator;}

void OpenWQ_wqconfig::set_output_units_denominator(double output_units_denominator)
{std::get<2>(this->output_units) = output_units_denominator;}

void OpenWQ_wqconfig::set_output_units_concentration(bool conentration_requested)
{std::get<3>(this->output_units) = conentration_requested;}

/***********************************************
 * Module-specific methods 
************************************************/

// #################
// LE_model -> Boundmix
int OpenWQ_wqconfig::LE_model_::BoundMix_::get_exchange_direction(unsigned entry_i)
{return std::get<0>(this->info_vector[entry_i]);}

int OpenWQ_wqconfig::LE_model_::BoundMix_::get_upper_compartment(unsigned entry_i)
{return std::get<1>(this->info_vector[entry_i]);}

int OpenWQ_wqconfig::LE_model_::BoundMix_::get_lower_compartment(unsigned entry_i)
{return std::get<2>(this->info_vector[entry_i]);}

int OpenWQ_wqconfig::LE_model_::BoundMix_::get_k_value(unsigned entry_i)
{return std::get<3>(this->info_vector[entry_i]);}

// #################
// TS_mdoel -> Hype MMF Erosion 
int OpenWQ_wqconfig::TS_model_::HypeMMF_::get_exchange_direction()
{return std::get<0>(this->info_vector);}

int OpenWQ_wqconfig::TS_model_::HypeMMF_::get_eroding_transp_compartment()
{return std::get<1>(this->info_vector);}

int OpenWQ_wqconfig::TS_model_::HypeMMF_::get_erosion_inhibit_compartment()
{return std::get<2>(this->info_vector);}

int OpenWQ_wqconfig::TS_model_::HypeMMF_::get_sediment_compartment()
{return std::get<3>(this->info_vector);}

std::string OpenWQ_wqconfig::TS_model_::HypeMMF_::get_data_format()
{return std::get<4>(this->info_vector);}

// #################
// TS_mdoel -> Hype HypeHVB Erosion 
int OpenWQ_wqconfig::TS_model_::HypeHVB_::get_exchange_direction()
{return std::get<0>(this->info_vector);}

int OpenWQ_wqconfig::TS_model_::HypeHVB_::get_eroding_transp_compartment()
{return std::get<1>(this->info_vector);}

int OpenWQ_wqconfig::TS_model_::HypeHVB_::get_erosion_inhibit_compartment()
{return std::get<2>(this->info_vector);}

int OpenWQ_wqconfig::TS_model_::HypeHVB_::get_sediment_compartment()
{return std::get<3>(this->info_vector);}

std::string OpenWQ_wqconfig::TS_model_::HypeHVB_::get_data_format()
{return std::get<4>(this->info_vector);}

/***********************************************
 * Performance: cache runtime flags to avoid
 * repeated string comparisons in hot paths
************************************************/
void OpenWQ_wqconfig::cache_runtime_flags()
{
    is_native_bgc_flex = (CH_model->BGC_module.compare("NATIVE_BGC_FLEX") == 0);
    is_TD_enabled = (TD_model->TD_module.compare("NONE") != 0);
    is_LE_enabled = (LE_model->LE_module.compare("NONE") != 0);
    is_TS_enabled = (TS_model->TS_module.compare("NONE") != 0);
    is_TD_advdisp = (TD_model->TD_module.compare("OPENWQ_NATIVE_TD_ADVDISP") == 0);
    is_TD_adv = (TD_model->TD_module.compare("NATIVE_TD_ADV") == 0);
    is_solver_sundials = (SOLVER_module.compare("SUNDIALS") == 0);
    is_solver_forward_euler = (SOLVER_module.compare("FORWARD_EULER") == 0);

    if (is_native_bgc_flex) {
        cached_num_mobile_species = CH_model->NativeFlex->mobile_species.size();
        cached_mobile_species_ptr = &(CH_model->NativeFlex->mobile_species);
        cached_num_chem = CH_model->NativeFlex->num_chem;
        cached_chem_species_list_ptr = &(CH_model->NativeFlex->chem_species_list);
    } else {
        cached_num_mobile_species = CH_model->PHREEQC->mobile_species.size();
        cached_mobile_species_ptr = &(CH_model->PHREEQC->mobile_species);
        cached_num_chem = CH_model->PHREEQC->num_chem;
        cached_chem_species_list_ptr = &(CH_model->PHREEQC->chem_species_list);
    }
}

/***********************************************
 * Performance: build BGC cycle name -> transformation
 * index lookup to avoid O(N) scan per compartment
************************************************/
void OpenWQ_wqconfig::build_bgc_lookup()
{
    bgc_cycle_to_transf_indices.clear();
    if (!is_native_bgc_flex) return;

    for (unsigned int j = 0; j < CH_model->NativeFlex->BGCexpressions_info.size(); j++) {
        const std::string& cycle_name = std::get<0>(CH_model->NativeFlex->BGCexpressions_info[j]);
        bgc_cycle_to_transf_indices[cycle_name].push_back(j);
    }
}

/***********************************************
 * Performance: Build per-thread copies of exprtk
 * expressions and their bound data vectors to
 * enable OpenMP parallelization of the BGC
 * spatial loops. Each thread gets its own
 * chemass_InTransfEq vector, dependVar_scalar
 * vector, and compiled expression objects.
 * Must be called AFTER bgc_flex_setBGCexpressions().
************************************************/
void OpenWQ_wqconfig::build_thread_local_expressions(
    OpenWQ_hostModelconfig& hostModelconfig)
{
    if (!is_native_bgc_flex) return;

    auto* nativeFlex = CH_model->NativeFlex;
    const unsigned int num_expressions = nativeFlex->BGCexpressions_eq.size();
    if (num_expressions == 0) return;

    // Determine number of threads
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    nativeFlex->num_omp_threads = nthreads;

    if (nthreads <= 1) {
        nativeFlex->thread_local_ready = false;
        return;
    }

    const unsigned int num_depend = hostModelconfig.get_num_HydroDepend();

    // Allocate per-thread storage
    nativeFlex->thread_chemass_InTransfEq.resize(nthreads);
    nativeFlex->thread_dependVar_scalar.resize(nthreads);
    nativeFlex->thread_BGCexpressions_eq.resize(nthreads);

    // For each thread, we need to:
    // 1. Create a local chemass vector (same size as the max needed)
    // 2. Create a local dependVar_scalar vector
    // 3. Re-compile each expression with a symbol table pointing to thread-local data

    // Find the maximum chemass vector size needed across all transformations
    unsigned int max_chemass_size = 0;
    for (unsigned int j = 0; j < num_expressions; j++) {
        const auto& bgc_info = nativeFlex->BGCexpressions_info[j];
        unsigned int sz = std::get<5>(bgc_info).size();
        if (sz > max_chemass_size) max_chemass_size = sz;
    }

    // Check that modified expression strings were stored during setBGCexpressions.
    // These contain the fully-substituted expressions (chemical names replaced by
    // vector indices, unit multipliers applied) needed to re-compile per-thread copies.
    if (nativeFlex->BGCexpressions_modif_strings.size() != num_expressions) {
        nativeFlex->thread_local_ready = false;
        return;
    }

    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    for (int t = 0; t < nthreads; t++) {
        // Allocate thread-local data vectors
        nativeFlex->thread_chemass_InTransfEq[t].resize(max_chemass_size, 0.0);
        nativeFlex->thread_dependVar_scalar[t].resize(num_depend, 0.0);

        // Compile each expression for this thread
        nativeFlex->thread_BGCexpressions_eq[t].resize(num_expressions);

        for (unsigned int j = 0; j < num_expressions; j++) {
            // Create symbol table bound to thread-local data
            symbol_table_t sym_table;

            // Resize thread-local chemass to match this expression's needs
            // (all threads share the same max_chemass_size allocation)
            sym_table.add_vector(
                "openWQ_BGCnative_chemass_InTransfEq",
                nativeFlex->thread_chemass_InTransfEq[t]);

            // Add dependency variables bound to thread-local scalars
            for (unsigned int depi = 0; depi < num_depend; depi++) {
                sym_table.add_variable(
                    hostModelconfig.get_HydroDepend_name_at(depi),
                    nativeFlex->thread_dependVar_scalar[t][depi]);
            }

            // Create and compile expression
            expression_t expr;
            expr.register_symbol_table(sym_table);

            parser_t parser;
            parser.compile(
                nativeFlex->BGCexpressions_modif_strings[j],
                expr);

            nativeFlex->thread_BGCexpressions_eq[t][j] = expr;
        }
    }

    nativeFlex->thread_local_ready = true;
}


