
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

#include "models_CH/headerfile_CH.hpp"
#include "global/OpenWQ_paramload.hpp"   // OpenWQ_load_param: scalar->GLOBAL / object->SPATIAL
#include <algorithm>
#include <cctype>
#include <cstdio>

/* #################################################
// Whole-identifier search inside a kinetics expression.
// A plain std::string::find would match a short symbol inside a longer
// one (parameter "T" inside "InTransfEq", "K" inside "K_NIT", species "N"
// inside "NH4-N"), silently corrupting the expression. A symbol only
// matches when it is not glued to identifier characters on either side.
// '-' is deliberately NOT an identifier character: it is the minus
// operator, and species names that contain it (NO3-N) are matched as a
// whole because longer names are substituted first.
################################################# */
static inline bool bgc_is_ident_char(char c){
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '.';
}
// Numeric literal substituted into an expression. std::to_string keeps only
// six decimal places, so a value such as 2.5e-7 became "0.000000" and a
// calibrated 0.000123456 became "0.000123": small parameters were switched off
// or quantized. %.17g keeps every significant digit of a double.
static std::string bgc_number_literal(double v){
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

static size_t bgc_find_whole_symbol(
    const std::string& text,
    const std::string& name,
    size_t from = 0){
    size_t p = text.find(name, from);
    while (p != std::string::npos){
        const bool left_ok  = (p == 0) || !bgc_is_ident_char(text[p - 1]);
        const size_t end    = p + name.size();
        const bool right_ok = (end >= text.size()) || !bgc_is_ident_char(text[end]);
        if (left_ok && right_ok) return p;
        p = text.find(name, p + 1);
    }
    return std::string::npos;
}


/* #################################################
// Parse biogeochemical expressions
################################################# */
void OpenWQ_CH_model::bgc_flex_setBGCexpressions(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output){

    // Local variables for expression evaluator: exprtk
    std::string 
        chemname,                   // chemical name
        consumed_spec,
        produced_spec,
        expression_string,          // expression string
        kinetics_units,             // kinetics units
        expression_string_modif;    // expression string replaced by variable in code
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    
    // Other local variables
    typedef std::tuple<
        std::string,                    // Biogeochemical cycle name
        std::string,                    // Transformation name
        std::string,                    // kinetic equation provided
        unsigned int,                   // index of consumed species       
        unsigned int,                   // index of produced species
        std::vector<unsigned int>       // index of chemical in transformation equation (needs to be here for loop reset)
        > BGCTransfTuple_info;          // Tuple with info and expression for BGC cyling
    std::vector<double> unit_multiplers;// multiplers (numerator and denominator)
    unsigned int num_BGCcycles, 
        num_transf,
        index_cons,index_prod,           // indexed for consumed and produced chemical
        index_new_chemass_InTransfEq;    // interactive index to build index_new_chemass_InTransfEq
    std::vector<std::string> BGCcycles_namelist;
    std::string BGCcycles_name, Transf_name;
    double param_val; // prameter value
    std::string msg_string;             // error/warning message string


    // Number of BCG cycles defined
    num_BGCcycles = OpenWQ_json.BGC_module["CYCLING_FRAMEWORKS"].size();

    // Reset SPATIAL-parameter storage (exprtk refactor). These stay empty /
    // max==0 for all-GLOBAL configs, keeping the expression strings and the
    // symbol tables byte-identical to the historical (literal-substituted) code.
    OpenWQ_wqconfig.CH_model->NativeFlex->BGCexpr_spatial_params.clear();
    OpenWQ_wqconfig.CH_model->NativeFlex->BGCparam_InTransfEq.clear();
    OpenWQ_wqconfig.CH_model->NativeFlex->BGCparam_InTransfEq.reserve(64);
    OpenWQ_wqconfig.CH_model->NativeFlex->max_BGCparam_size = 0;

    // Get all biogeochemical cycling names
    for (auto it: OpenWQ_json.BGC_module["CYCLING_FRAMEWORKS"].items())
    {
        BGCcycles_namelist.push_back(it.key());
    }

    /* ########################################
    // Loop over biogeochemical cycling frameworks
    ######################################## */
    for (unsigned int bgci=0;bgci<num_BGCcycles;bgci++){

        /* ########################################
        // Loop over transformations in biogeochemical cycle bgci
        ######################################## */

        // Get BGC cycle name
        BGCcycles_name = BGCcycles_namelist[bgci];

        // Get number of transformations inside BGCcycles_name
        num_transf = OpenWQ_json.BGC_module
            ["CYCLING_FRAMEWORKS"]
            [BGCcycles_name]
            ["LIST_TRANSFORMATIONS"].size();
        
        // The species-mass vector that every compiled expression is bound to
        // (exprtk add_vector keeps a view of the vector's CURRENT buffer) is
        // allocated ONCE here, to the largest size any expression can need
        // (all species), and is never cleared, grown or resized afterwards.
        // Growing it per transformation reallocated the buffer, so every
        // expression compiled before the last reallocation kept a dangling
        // view and read its species masses as 0 (first-order kinetics inert).
        if (OpenWQ_wqconfig.CH_model->NativeFlex->chemass_InTransfEq.size()
                < OpenWQ_wqconfig.CH_model->NativeFlex->num_chem)
            OpenWQ_wqconfig.CH_model->NativeFlex->chemass_InTransfEq.assign(
                OpenWQ_wqconfig.CH_model->NativeFlex->num_chem, 0.0);

        for (unsigned int transi=0;transi<num_transf;transi++){

            // Get Transformation name
            Transf_name = OpenWQ_json.BGC_module
                ["CYCLING_FRAMEWORKS"]
                [BGCcycles_name]
                ["LIST_TRANSFORMATIONS"]
                [std::to_string(transi+1)];

            std::vector<unsigned int> index_transf; // index of chemical in transformation equation (needs to be here for loop reset)

            // Get transformation transi info
            consumed_spec =  OpenWQ_json.BGC_module
                ["CYCLING_FRAMEWORKS"]
                [BGCcycles_name]
                [std::to_string(transi+1)]
                ["CONSUMED"];
            produced_spec =  OpenWQ_json.BGC_module
                ["CYCLING_FRAMEWORKS"]
                [BGCcycles_name]
                [std::to_string(transi+1)]
                ["PRODUCED"];
            expression_string = OpenWQ_json.BGC_module
                ["CYCLING_FRAMEWORKS"]
                [BGCcycles_name]
                [std::to_string(transi+1)]
                ["KINETICS"].at(0);
            kinetics_units = OpenWQ_json.BGC_module
                ["CYCLING_FRAMEWORKS"]
                [BGCcycles_name]
                [std::to_string(transi+1)]
                ["KINETICS"].at(1);
            std::vector<std::string> parameter_names = OpenWQ_json.BGC_module
                ["CYCLING_FRAMEWORKS"]
                [BGCcycles_name]
                [std::to_string(transi+1)]
                ["PARAMETER_NAMES"];

            expression_string_modif = expression_string;

            /* ########################################
            // Adjust expression to change time units to sec
            ######################################## */
            
            // Calculate unit multiplers
            std::vector<std::string> units;          // units (numerator and denominator)
            OpenWQ_units.Calc_Unit_Multipliers(
                OpenWQ_wqconfig,
                OpenWQ_output,
                unit_multiplers,    // multiplers (numerator and denominator)
                kinetics_units,     // input units
                units,
                true);              // direction of the conversion: 
                                    // to native (true) or 
                                    // from native to desired output units (false)
            
            // Adjuct expression
            expression_string_modif = 
                "(" + expression_string_modif 
                + ")*" + bgc_number_literal(unit_multiplers[0])
                + "/" + bgc_number_literal(unit_multiplers[1]);

            /* ########################################
            // Find species indexes: consumed, produced and in the expression
            ######################################## */
            // these indexes need to be reset here
            index_cons = -1;
            index_prod = -1; // indexed for consumed and produced chemical
            index_new_chemass_InTransfEq = 0;

            for(unsigned int chemi=0;chemi<(OpenWQ_wqconfig.CH_model->NativeFlex->num_chem);chemi++){

                // Get chemical species name
                chemname = (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi];

                // Consumedchemass_consumed, chemass_produced;ty())
                if(consumed_spec.compare(chemname) == 0 && !consumed_spec.empty()){
                    index_cons = chemi; // index
                }

                // Produced
                if(produced_spec.compare(chemname) == 0 && !produced_spec.empty()){
                    index_prod = chemi; // index
                }
            }

            // ########################################
            // Substitute species and parameters in the expression.
            //   species       -> element of the bound species-mass vector
            //                    (one element per occurrence, in the order the
            //                    occurrences are substituted; index_transf
            //                    records which species each element holds)
            //   GLOBAL param  -> its literal value (std::to_string)
            //   SPATIAL param -> element k of the bound openWQ_BGCparam vector,
            //                    refreshed per cell
            // A parameter is SPATIAL when its OpenWQ_param carries a per-cell
            // field (JSON object); a JSON number is GLOBAL.
            //
            // All symbols (species + parameters) are substituted longest name
            // first, whole identifiers only, every occurrence. This is what
            // makes the substitution safe for one-letter parameters (T, K),
            // parameters that prefix each other (K and K_NIT), species that
            // contain other species (N inside NH4-N) and parameters that
            // appear more than once in the same expression.
            // ########################################
            std::vector<OpenWQ_param> expr_spatial_params; // this expression's spatial params (index == k)

            // symbol list: (name, species index) or (name, -(param index)-1)
            std::vector<std::pair<std::string,int>> symbols;
            for(unsigned int chemi=0;chemi<(OpenWQ_wqconfig.CH_model->NativeFlex->num_chem);chemi++)
                symbols.emplace_back(
                    (OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list)[chemi],
                    static_cast<int>(chemi));
            for (unsigned int i=0;i<parameter_names.size();i++)
                symbols.emplace_back(parameter_names[i], -static_cast<int>(i)-1);
            std::stable_sort(symbols.begin(), symbols.end(),
                [](const std::pair<std::string,int>& a, const std::pair<std::string,int>& b){
                    return a.first.size() > b.first.size();
                });

            for (const auto& sym : symbols){
                const std::string& name = sym.first;
                if (name.empty()) continue;

                if (sym.second >= 0){
                    // Species: replace every whole-word occurrence
                    const unsigned int chemi = static_cast<unsigned int>(sym.second);
                    size_t pos = bgc_find_whole_symbol(expression_string_modif, name);
                    while (pos != std::string::npos){
                        const std::string repl =
                            "openWQ_BGCnative_chemass_InTransfEq["
                            + std::to_string(index_new_chemass_InTransfEq) + "]";
                        index_transf.push_back(chemi);
                        expression_string_modif.replace(pos, name.size(), repl);
                        index_new_chemass_InTransfEq ++;
                        pos = bgc_find_whole_symbol(expression_string_modif, name, pos + repl.size());
                    }
                    continue;
                }

                // Parameter: read its JSON entry and build an OpenWQ_param
                // (a number -> GLOBAL scalar; an object such as {"UNIFORM":v}
                // or {"DEFAULT":d,"CELLS":[...]} -> SPATIAL).
                // A name listed in PARAMETER_NAMES is looked up first in the
                // reaction's PARAMETER_VALUES and, when absent there, in the
                // module-level GLOBAL_PARAMETERS block. A global parameter is
                // one value shared by several reactions (e.g. the rate constant
                // of a reaction written as two transformations, one per species
                // it consumes), so calibrating it keeps them consistent.
                const unsigned int i = static_cast<unsigned int>(-sym.second - 1);
                const json& rxn_json = OpenWQ_json.BGC_module
                    ["CYCLING_FRAMEWORKS"]
                    [BGCcycles_name]
                    [std::to_string(transi+1)];
                json param_jval;
                if (rxn_json.contains("PARAMETER_VALUES")
                        && rxn_json["PARAMETER_VALUES"].contains(parameter_names[i])){
                    param_jval = rxn_json["PARAMETER_VALUES"][parameter_names[i]];
                } else if (OpenWQ_json.BGC_module.contains("GLOBAL_PARAMETERS")
                        && OpenWQ_json.BGC_module["GLOBAL_PARAMETERS"].contains(parameter_names[i])){
                    param_jval = OpenWQ_json.BGC_module["GLOBAL_PARAMETERS"][parameter_names[i]];
                } else {
                    msg_string = "<OpenWQ> FATAL: parameter '" + parameter_names[i]
                        + "' of CYCLING_FRAMEWORKS > " + BGCcycles_name + " > "
                        + std::to_string(transi+1) + " (" + Transf_name + ") has no "
                        "value in PARAMETER_VALUES nor in GLOBAL_PARAMETERS.";
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    exit(EXIT_FAILURE);
                }

                OpenWQ_param param_i = OpenWQ_load_param(
                    param_jval, OpenWQ_hostModelconfig);
                param_val = param_i.scalar(); // scalar value (GLOBAL substitution + logging)

                std::string repl;
                if (!param_i.is_spatial()){
                    repl = "(" + bgc_number_literal(param_val) + ")";
                } else {
                    repl = "openWQ_BGCparam[" + std::to_string(expr_spatial_params.size()) + "]";
                }

                size_t pos = bgc_find_whole_symbol(expression_string_modif, name);
                if (pos == std::string::npos){
                    // Declared but not used in this expression
                    msg_string = "<OpenWQ> Parameter ignored in CYCLING_FRAMEWORKS > "
                        + BGCcycles_name + " > "
                        + std::to_string(transi+1)
                        + ". Parameter " + parameter_names[i]
                        + " with value " + std::to_string(param_val);
                    OpenWQ_output.ConsoleLog(
                        OpenWQ_wqconfig,    // for Log file name
                        msg_string,         // message
                        true,               // print in console
                        true);              // print in log file
                    continue;
                }
                if (param_i.is_spatial()) expr_spatial_params.push_back(param_i);
                while (pos != std::string::npos){
                    expression_string_modif.replace(pos, name.size(), repl);
                    pos = bgc_find_whole_symbol(expression_string_modif, name, pos + repl.size());
                }
            }

            // Add variables to symbol_table
            symbol_table_t symbol_table;
            
            // Bind the (fixed-size, never reallocated) species-mass vector; the
            // expression only indexes the first index_transf.size() entries.
            symbol_table.add_vector("openWQ_BGCnative_chemass_InTransfEq",OpenWQ_wqconfig.CH_model->NativeFlex->chemass_InTransfEq);

            // Bind the SPATIAL-parameter vector too (only when this expression
            // references it). BGCparam_InTransfEq was reserve()d up front, so
            // growing it to hold this expression's params does not reallocate
            // and the pointer exprtk stores stays valid. When there are no
            // spatial params openWQ_BGCparam is never bound - the symbol table
            // is byte-identical to the historical one.
            if (!expr_spatial_params.empty()){
                auto* nf_ = OpenWQ_wqconfig.CH_model->NativeFlex;
                if (nf_->BGCparam_InTransfEq.size() < expr_spatial_params.size())
                    nf_->BGCparam_InTransfEq.resize(expr_spatial_params.size(), 0.0);
                symbol_table.add_vector("openWQ_BGCparam", nf_->BGCparam_InTransfEq);
            }

            // Add variable dependencies to table of symbols (in case they are used).
            // exprtk binds a variable by reference, so it must be bound to the
            // scalar that the transform updates for every cell. It used to be
            // bound to a local copy, which went out of scope here: every
            // expression evaluated on the serial path (one thread, or a
            // compartment with a single cell such as a lumped aquifer) read a
            // dangling value instead of the cell's temperature or moisture.
            for (unsigned int depi=0;depi<OpenWQ_hostModelconfig.get_num_HydroDepend();depi++){

                symbol_table.add_variable(
                    OpenWQ_hostModelconfig.get_HydroDepend_name_at(depi),    // Dependency Var name
                    OpenWQ_hostModelconfig.get_dependVar_scalar_ref(depi)    // Variable data
                );

            }

            // Water volume of the cell [m3] (see Vw_m3 in OpenWQ_wqconfig.hpp)
            symbol_table.add_variable(
                "Vw_m3", OpenWQ_wqconfig.CH_model->NativeFlex->cell_watervol_m3);
            
            // Create Object
            expression_t expression;
            expression.register_symbol_table(symbol_table);

            // Parse expression and compile 
            parser_t parser;
            if (!parser.compile(expression_string_modif,expression)) {
                // A kinetics expression that does not compile (typically a
                // symbol that is neither a species, a declared parameter nor
                // a dependency) would silently evaluate to NaN -> 0 mass
                // transferred, i.e. the reaction is quietly switched off.
                // That is a configuration error: stop and say which one.
                msg_string = "<OpenWQ> FATAL: the KINETICS expression of "
                    "CYCLING_FRAMEWORKS > " + BGCcycles_name + " > "
                    + std::to_string(transi+1) + " (" + Transf_name + ") does not "
                    "compile: '" + expression_string + "'. exprtk: "
                    + parser.error() + ". Every symbol must be a chemical "
                    "species, a name listed in PARAMETER_NAMES (with a value "
                    "in PARAMETER_VALUES), or a declared dependency.";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                exit(EXIT_FAILURE);
            }

            // Save expressions in openWQ_BGCnative_BGCexpressions_info and openWQ_BGCnative_BGCexpressions_eq
            OpenWQ_wqconfig.CH_model->NativeFlex->BGCexpressions_info.push_back(
                BGCTransfTuple_info(
                    BGCcycles_name,
                    Transf_name,
                    expression_string,
                    index_cons,
                    index_prod,
                    index_transf));

            OpenWQ_wqconfig.CH_model->NativeFlex->BGCexpressions_eq.push_back(expression);

            // PARALLEL: Store modified expression string for re-compilation
            // in per-thread copies (needed for OpenMP parallelization)
            OpenWQ_wqconfig.CH_model->NativeFlex->BGCexpressions_modif_strings.push_back(
                expression_string_modif);

            // Store this expression's SPATIAL-parameter list (index-aligned with
            // BGCexpressions_eq / _modif_strings) and track the max count so the
            // per-thread param vectors can be sized once at thread-local init.
            OpenWQ_wqconfig.CH_model->NativeFlex->BGCexpr_spatial_params.push_back(
                expr_spatial_params);
            if (expr_spatial_params.size() >
                    OpenWQ_wqconfig.CH_model->NativeFlex->max_BGCparam_size)
                OpenWQ_wqconfig.CH_model->NativeFlex->max_BGCparam_size =
                    expr_spatial_params.size();

        }
    }
}
