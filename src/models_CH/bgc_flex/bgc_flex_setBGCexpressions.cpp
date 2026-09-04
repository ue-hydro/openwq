
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
    int index_i;                         // iteractive index (can be zero because it is determined in a .find())
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
                + ")*" + std::to_string(unit_multiplers[0])
                + "/" + std::to_string(unit_multiplers[1]);

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

                // In expression (replace chemical species name by index)
                // while loop to replace all occurances
                index_i = expression_string_modif.find(chemname);
                while(index_i!=-1){
                    if (index_i!=-1 && !expression_string_modif.empty()){
                        index_transf.push_back(chemi); // index
                        expression_string_modif.replace(
                            index_i,    // start position
                            chemname.size(),    // length
                            "openWQ_BGCnative_chemass_InTransfEq["+std::to_string(index_new_chemass_InTransfEq)+"]"); // string to replace by
                        index_new_chemass_InTransfEq ++;    
                    }
                    index_i = expression_string_modif.find(chemname);
                }
            }

            // ########################################
            // Replace parameters in the expression.
            //   GLOBAL param  -> substitute its literal value (historical path,
            //                    byte-identical via std::to_string).
            //   SPATIAL param -> leave it in the expression as element k of the
            //                    bound openWQ_BGCparam vector, refreshed per cell.
            // A parameter is SPATIAL when its OpenWQ_param carries a per-cell
            // field. Today values are JSON scalars (all GLOBAL) until the loader
            // / ML layer populates a spatial field, so this stays byte-identical.
            // ########################################
            std::vector<OpenWQ_param> expr_spatial_params; // this expression's spatial params (index == k)
            for (unsigned int i=0;i<parameter_names.size();i++){
                index_i = expression_string_modif.find(parameter_names[i]);

                // Read the parameter's JSON entry and build an OpenWQ_param:
                // a number -> GLOBAL scalar (historical); an object (e.g.
                // {"UNIFORM":v} or {"DEFAULT":d,"CELLS":[...]}) -> SPATIAL.
                json param_jval = OpenWQ_json.BGC_module
                    ["CYCLING_FRAMEWORKS"]
                    [BGCcycles_name]
                    [std::to_string(transi+1)]
                    ["PARAMETER_VALUES"]
                    [parameter_names[i]];

                OpenWQ_param param_i = OpenWQ_load_param(
                    param_jval, OpenWQ_hostModelconfig);
                param_val = param_i.scalar(); // scalar value (GLOBAL substitution + logging)

                // Try replacing
                try{
                    if (!param_i.is_spatial()){
                        // GLOBAL: literal substitution (unchanged behaviour)
                        expression_string_modif.replace(
                            index_i,
                            parameter_names[i].size(),
                            std::to_string(param_val));
                    } else {
                        // SPATIAL: reference the bound vector element k
                        const unsigned int k = expr_spatial_params.size();
                        expression_string_modif.replace(
                            index_i,
                            parameter_names[i].size(),
                            "openWQ_BGCparam[" + std::to_string(k) + "]");
                        expr_spatial_params.push_back(param_i);
                    }
                }
                catch(...){

                    // Create Message
                    msg_string = "<OpenWQ> Parameter ignored in CYCLING_FRAMEWORKS > "
                        + BGCcycles_name + " > "
                        + std::to_string(transi+1)
                        + ". Parameter " + parameter_names[i]
                        + " with value " + std::to_string(param_val);

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(
                        OpenWQ_wqconfig,    // for Log file name
                        msg_string,         // message
                        true,               // print in console
                        true);              // print in log file

                };
            }

            // Add variables to symbol_table
            symbol_table_t symbol_table;
            
            OpenWQ_wqconfig.CH_model->NativeFlex->chemass_InTransfEq.clear();
            for (unsigned int i=0;i<index_transf.size();i++){
                OpenWQ_wqconfig.CH_model->NativeFlex->chemass_InTransfEq.push_back(0); // creating the vector
            }

            // Add vectors to table of symbols
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

            // Add variable dependencies to table of symbols (in case they are used)
            for (unsigned int depi=0;depi<OpenWQ_hostModelconfig.get_num_HydroDepend();depi++){
                
                double var = OpenWQ_hostModelconfig.get_dependVar_scalar_at(depi);

                symbol_table.add_variable(
                    OpenWQ_hostModelconfig.get_HydroDepend_name_at(depi),    // Dependency Var name
                    var    // Variable data
                );

            }
            
            // Create Object
            expression_t expression;
            expression.register_symbol_table(symbol_table);

            // Parse expression and compile 
            parser_t parser;
            parser.compile(expression_string_modif,expression);

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
