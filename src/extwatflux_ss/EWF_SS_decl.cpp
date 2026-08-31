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

#include "headerfile_EWF_SS.hpp"
#include <cstdlib>   // std::exit, EXIT_FAILURE (EWF-HDF5 SPATIAL_MODE validation)
#include <cmath>     // std::llround (DISTRIBUTED reach-id matching)


/* #################################################
 // Prepare SS and EWF input data for use at running time: main driver
 ################################################# */
void OpenWQ_extwatflux_ss::Set_EWFandSS_driver(
    json &EWF_SS_json,                                  // SS or EWF json    
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output,
    std::string inputType) {                             // flag for SS or EWF
    
    // Local variables
    bool foundflag = false;                     // iteractive boolean to identify if comp or chem was found
    unsigned int num_srcfiles;                  // number of sink-source files 
                                                // (saved as sub-structure of SinkSource)
    unsigned int num_srchem;                    // number of chemical loads per file
    std::string DataFormat;                     // from JSON file (JSON or ASCII)
    std::string main_keyName;                   // interactive json-key name
    std::string msg_string;                     // error/warning message string


    // Get number of sub-structures of SS/EWF data
    num_srcfiles = EWF_SS_json.size();

    /* ########################################
    // Loop over file (saved as sub-structure of SS/EWF data)
    ######################################## */

    for (unsigned int ssf=0;ssf<num_srcfiles;ssf++){
        
        // Get number of loads in each sub-structure
        // (corresponding to different SinkSource/ExtWatFlux json files)
        num_srchem = EWF_SS_json[std::to_string(ssf+1)].size();
        
        /* ########################################
        // Loop over loads per sub-structure
        ######################################## */

        foundflag = false;      // reset to false at every 

        for (unsigned int ssi=0;ssi<num_srchem;ssi++){

            /* ########
            // Get data format
            // try-catch block because we may have other sub-structures
            // Example: Metadata 
            ###########*/
            try{
                DataFormat = EWF_SS_json // units
                    [std::to_string(ssf+1)]
                    [std::to_string(ssi+1)]
                    ["DATA_FORMAT"];
            }catch(...){continue;}

            /* ########
            // Call appropriate function depending on data format
            ###########*/
            if (DataFormat.compare("JSON")==0 || 
                DataFormat.compare("ASCII")==0) {
                
                // if JSON or ASCII format
                Set_EWFandSS_jsonAscii(
                    OpenWQ_hostModelconfig,
                    OpenWQ_wqconfig,
                    OpenWQ_utils,
                    OpenWQ_units,
                    OpenWQ_output,
                    ssf, ssi,
                    DataFormat,
                    EWF_SS_json[std::to_string(ssf+1)][std::to_string(ssi+1)],  // relevant sub-json
                    inputType, // ss or ewf
                    foundflag);

            }else if(DataFormat.compare("HDF5")==0){
                
                // h5 format only supported for ewf
                if(inputType.compare("ewf")==0){

                    // if H5 format
                    Set_EWF_h5(
                        OpenWQ_hostModelconfig,
                        OpenWQ_wqconfig,
                        OpenWQ_utils,
                        OpenWQ_units,
                        OpenWQ_output,
                        EWF_SS_json[std::to_string(ssf+1)][std::to_string(ssi+1)],  // relevant sub-json
                        inputType, // ss or ewf
                        foundflag);

                }else if(inputType.compare("ss")==0){

                    // Create Message (Warning Message)
                    msg_string = 
                        "<OpenWQ> WARNING: HDF5 input only supported for'" 
                        " EWF forcing. For SS forcing, please use JSON or "
                        "ASCII format (entry skipped).";

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true); 

                }

            }else{

                // Create Message (Warning Message)
                msg_string =
                    "<OpenWQ> WARNING: Unkown data format='"
                    + DataFormat
                    + "' in SS or EWF json files > "
                    + " > DATA_FORMAT (only supports JSON, ASCII and HD5F) "
                    + "(entry skipped)";

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

            }
        }
    }

    // OPTIMIZED (perf): all JSON/ASCII rows for this inputType have now been
    // staged in a flat buffer; build the FORC matrix in a single allocation
    // (see AppendRow_SS_EWF_FORC_jsonAscii / Finalize_FORC_jsonAscii).
    Finalize_FORC_jsonAscii(OpenWQ_wqconfig, inputType);
}


/* #################################################
 // Prepare SS and EWF input data for use at running time: 
 // Case: JSON and ASCII format
 ################################################# */
void OpenWQ_extwatflux_ss::Set_EWFandSS_jsonAscii(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    unsigned int ssf, unsigned int ssi,   // file-structure and substructure indexes
    std::string DataFormat,         // (JSON or ASCII)
    json EWF_SS_json_sub,           // relevant sub-json
    std::string inputType,
    bool foundflag){                // ss or ewf

    // Local variables
    std::string main_keyName;
    std::string Element_name;
    std::string Chemical_name;                  // chemical name
    std::vector<std::string> elm_list;          // model compartment list
    std::string err_text;                       // iteractive string for text to pass to error messages
    unsigned long cmpi_ssi;                     // model index for compartment Compartment_name_name
    unsigned long chem_ssi;                     // model index for compartment Compartment_name_name
    unsigned long sinksource_ssi;               // = 0 (source), = 1 (sink)
    unsigned int num_rowdata = 0;               // number of rows of data in JSON (YYYY, MM, DD, HH,...)
    std::string Type;                           // from JSON filec (only used in SS)
    std::string ascii_FilePath;                 // additional information for ASCII data input
    std::string ascii_delimiter;                // additional information for ASCII data input
    unsigned int lineCount;                     // row count for ASCII file
    std::string rowEntryASCII;                  // row data
    std::string headerRowASCII;                 // row with header data
    std::vector<std::string> ASCIIRowdata;      // all row data
    std::vector<std::string> ASCIIRowElemEntry; // vector with header keys
    std::vector<std::string> headerKeys;        // vector with header keys
    std::string headerWord;                     // interactive header words
    std::string elemName;                       // temporary element name
    json EWF_SS_json_sub_rowi;                  // json row of json-substructure EWF_SS_json_sub 
    int YYYY_json;                              // Year in JSON-sink_source (interactive)
    int MM_json;                                // Month in JSON-sink_source (interactive)
    int DD_json;                                // Day in JSON-sink_source (interactive)
    int HH_json;                                // Hour in JSON-sink_source (interactive)
    int MIN_json;                               // Minutes in JSON-sink_source (interactive)
    int SEC_json;                               // Seconds in JSON-sink_source (interactive)
    int ix_json;                                // iteractive ix info for sink-source row data 
    int iy_json;                                // iteractive iy info for sink-source row data 
    int iz_json;                                // iteractive iz info for sink-source row data
    double ss_data_json;                        // data (sink or source) from row data
    std::string ss_units_json;                  // units of row data
    std::string ss_units_json_mass_base;             // units of row data (mass)
    std::vector<double> unit_multiplers;        // multiplers (numerator and denominator)
    arma::vec row_data_col;                     // new row data (initially as col data)
    std::string msg_string;                     // error/warning message string
    std::string errorMsgIdentifier;             // error message section identifier
    bool validEntryFlag;                        // valid entry flag to skip problematic row data
    std::string loadScheme_str;                 // Load scheme string: (1) discrete or (2) continuous
    double loadScheme_id = 9999;                // Load scheme id number: (1) discrete or (2) continuous
    std::string contDt_str;                     // time units of continuous load
    bool cell_id_lookup_used = false;           // Flag to track if cell_id mapping was used for ix/iy/iz


    // Get element list (compartment for SS, and External fluxes for EWF)
    if (inputType.compare("ss")==0)     elm_list = OpenWQ_hostModelconfig.get_HydroComp_names();
    if (inputType.compare("ewf")==0)    elm_list = OpenWQ_hostModelconfig.get_HydroExtFlux_names();

    // Element names: Compartment or External Flux
    if (inputType.compare("ss")==0)       main_keyName = "COMPARTMENT_NAME"; // SS
    else if (inputType.compare("ewf")==0) main_keyName = "EXTERNAL_INPUTFLUX_NAME";    // EWF

    // Get mainKeyName
    // "COMPARTMENT_NAME" in the case of "ss"
    // "EXTERNAL_INPUTFLUX_NAME" in the case of "ewf"
    // Needs try-catch because there may be other irrelevant entries e.g. COMMENTS
    try{
        Element_name = EWF_SS_json_sub[main_keyName];
    }catch(...){
        return;
    }

    /* ########
    // Get chemical name and compartment/external-flux name 
    // if SS, then get also SS_type (source or sink)
    ###########*/
    errorMsgIdentifier = inputType + " json block";
    Chemical_name = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "CHEMICAL_NAME",
        errorMsgIdentifier,
        true);
            
    // Type (sink or source) (only used in SS)
    if (inputType.compare("ss")==0) {
        errorMsgIdentifier = inputType + " json block";
        Type = OpenWQ_utils.RequestJsonKeyVal_str(
            OpenWQ_wqconfig, OpenWQ_output,
            EWF_SS_json_sub, "TYPE",
            errorMsgIdentifier,
            true);
    }else if (inputType.compare("ewf")==0){
        Type = "SOURCE";}

    /* ########
    // Check if the requests are valid
    // chemical name, compartment name and SS_type (source or sink)
    ###########*/

    // Get chemical index
    err_text.assign("Chemical name");
    if ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0) {
        foundflag = getModIndex(
            OpenWQ_wqconfig,
            OpenWQ_output,
            OpenWQ_wqconfig.CH_model->NativeFlex->chem_species_list,
            Chemical_name,
            err_text,
            chem_ssi);
    } else {
        foundflag = getModIndex(
            OpenWQ_wqconfig,
            OpenWQ_output,
            OpenWQ_wqconfig.CH_model->PHREEQC->chem_species_list,
            Chemical_name,
            err_text,
            chem_ssi);
    }

    // Get Units
    errorMsgIdentifier = inputType + " json block";
    ss_units_json_mass_base= OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "UNITS",
        errorMsgIdentifier,
        true);
    std::transform(ss_units_json_mass_base.begin(), ss_units_json_mass_base.end(),
                  ss_units_json_mass_base.begin(), ::toupper);

    if (foundflag == false) return; // skip if chem not found

    // Get compartment/water flux index
    if (inputType.compare("ss")==0) err_text.assign("Compartment name");
    if (inputType.compare("ewf")==0) err_text.assign("External Water Flux name");

    foundflag = getModIndex(
        OpenWQ_wqconfig,
        OpenWQ_output,
        elm_list,
        Element_name,
        err_text,
        cmpi_ssi);

    if (foundflag == false) return; // skip if comp/ext-flux not found

    // Set type flag (sink or source)
    // if EWF, then Type has been defined above as "SOURCE"
    std::transform(Type.begin(), Type.end(), Type.begin(), ::toupper);
    if (Type.compare("SOURCE") == 0){ sinksource_ssi = 0;
    }else if (Type.compare("SINK") == 0){ sinksource_ssi = 1;
    }else{return;} // skip if Type is unknown
    
    /* ########
    // Check if DATA_FORMAT=ASCII
    // if yes, then get file structure info
    ###########*/

    // Get additional info for ASCII
    if (DataFormat.compare("ASCII")==0){

        // Check if substructure DATA exists
        errorMsgIdentifier = inputType + " json block with DataFormat=" + DataFormat;
        OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            EWF_SS_json_sub, "DATA",
            errorMsgIdentifier,
            true);

        try{

            // Error msg identifier in case json key not found
            errorMsgIdentifier = inputType + " json block 'DATA' with DataFormat=" + DataFormat;

            // file path
            ascii_FilePath = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                EWF_SS_json_sub["DATA"], "FILEPATH",
                errorMsgIdentifier,
                true);

            // delimiter
            ascii_delimiter = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                EWF_SS_json_sub["DATA"], "DELIMITER",
                errorMsgIdentifier,
                true);

            // Open ASCII file and auto-detect header row by scanning for "YYYY"
            {
                std::ifstream asciiFile (ascii_FilePath, std::ios::in);

                if (asciiFile.is_open()){

                    bool headerFound = false;
                    lineCount = 1;

                    while(std::getline(asciiFile, rowEntryASCII)){

                        // Convert line to uppercase for header detection
                        std::string lineUpper = OpenWQ_utils.ConvertStringToUpperCase(rowEntryASCII);

                        // Auto-detect: header row is the first line containing "YYYY"
                        if (!headerFound && lineUpper.find("YYYY") != std::string::npos){
                            headerFound = true;
                            headerRowASCII = lineUpper;
                            headerKeys = OpenWQ_utils.StringSplit(headerRowASCII, ascii_delimiter);
                        }
                        // Everything after the header row is data
                        else if (headerFound){
                            ASCIIRowdata.push_back(rowEntryASCII);
                        }

                        ++lineCount;
                    }
                    asciiFile.close();

                    // Check that we found a header
                    if (!headerFound){
                        msg_string = "<OpenWQ> WARNING: SS/EWF ASCII file has no header row "
                            "containing 'YYYY' (entry skipped): File=" + ascii_FilePath;
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        return;
                    }

                    num_rowdata = ASCIIRowdata.size();

                }else{
                    // If there is an issue with the ASCII input data
                    // through a warning message and skip entry
                    msg_string =
                        "<OpenWQ> WARNING: SS/EWF '"
                        " load/sink/conc ASCII file cannot be found (entry skipped): File="
                        + ascii_FilePath
                        + " in JSON SS file "
                        + std::to_string(ssf+1)
                        + ", Sub_structure=" + std::to_string(ssi+1);
                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    return;
                }
            } // end scope block for ASCII file reading

        }catch(...){
            // If there is an issue with the ASCII input data
            // through a warning message and skip entry
            msg_string = 
                "<OpenWQ> WARNING: SS/EWF '" 
                " load/sink/conc with ASCII format has an issue with json-keys or data structure (entry skipped): File=" 
                + std::to_string(ssf+1)
                + ", Sub_structure=" + std::to_string(ssi+1);
            // Print it (Console and/or Log file)
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            return;
        }
        
    // If JSON
    }else if (DataFormat.compare("JSON")==0){
        // Get number of rows of data in JSON (YYYY, MM, DD, HH,...)
        errorMsgIdentifier = inputType + " json block 'DATA' with DataFormat=" + DataFormat;
        // check if field exists, return if not with error message
        num_rowdata = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            EWF_SS_json_sub, "DATA",
            errorMsgIdentifier,
            true).size();
    }

    /* ########################################
    // Loop over row data in sink-source file
    ######################################## */

    // OPTIMIZED (perf): resolve the "DATA" block once before the loop instead of
    // doing EWF_SS_json_sub["DATA"] on every row. Combined with
    // RequestJsonKeyVal_json now taking its argument by const reference, this
    // removes a full deep-copy of the entire DATA block on every row iteration
    // (previously O(rows^2) — the dominant sink-source ingestion cost).
    json EWF_SS_DATA_block;
    if (DataFormat.compare("JSON")==0){
        EWF_SS_DATA_block = EWF_SS_json_sub["DATA"];
    }

    for (unsigned int di=0;di<num_rowdata;di++){

        // Reset the size to zero (the object will have no elements)
        row_data_col.reset();

        // Get row-json di from EWF_SS_json_sub ["DATA"]
        // Only needed for JSON format (ASCII reads from file directly)
        if (DataFormat.compare("JSON")==0){
            errorMsgIdentifier = inputType + " json block 'DATA', row "
                                + std::to_string(di) + " with DataFormat="
                                + DataFormat;
            EWF_SS_json_sub_rowi = OpenWQ_utils.RequestJsonKeyVal_json(
                OpenWQ_wqconfig, OpenWQ_output,
                EWF_SS_DATA_block, std::to_string(di+1),
                errorMsgIdentifier,
                true);
        }

        // If DataFormat=ASCII, then get row data
        // and convert-to-upper-case and split it by element entry
        if (DataFormat.compare("ASCII")==0){
            ASCIIRowElemEntry =
                OpenWQ_utils.StringSplit(
                    OpenWQ_utils.ConvertStringToUpperCase(ASCIIRowdata[di]),
                    ascii_delimiter);
            // Safety: skip row if header or data parsing failed
            if (ASCIIRowElemEntry.empty() || headerKeys.empty()){
                msg_string = "<OpenWQ> WARNING: SS/EWF ASCII row "
                    + std::to_string(di) + " could not be parsed (empty header or data). Skipping.";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                continue;
            }
        }

        // ###################
        // Year
        // ###################
        elemName = "Year";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(0);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"YYYY")]);}

            YYYY_json = entryVal;

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(0);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"YYYY")];}

            // Check if "all" and return flag validEntryFlag
            validEntryFlag = getArrayElem_SS(
                OpenWQ_wqconfig,OpenWQ_output,
                elemName,
                (std::string) entryVal,
                YYYY_json,
                ssf, ssi, di); // SS file, structure and row
                    
            if (!validEntryFlag){continue;}

        }
        
        // ###################
        // Month
        // ###################
        elemName = "Month";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(1);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"MM")]);}

            MM_json = entryVal;

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(1);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"MM")];}

            validEntryFlag = getArrayElem_SS(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                MM_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }
        
        // ###################
        // Day
        // ###################
        elemName = "Day";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(2);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"DD")]);}

            DD_json = entryVal;

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(2);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"DD")];}

            validEntryFlag = getArrayElem_SS(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                DD_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }

        // ###################
        // Hour
        // ###################
        elemName = "Hour";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(3);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"HH")]);}

            HH_json = entryVal;

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(3);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"HH")];}

            validEntryFlag = getArrayElem_SS(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                HH_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }        
        
        // ###################
        // Minute
        // ###################
        elemName = "Minute";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(4);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"MIN")]);}

            MIN_json = entryVal;

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(4);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"MIN")];}

            validEntryFlag = getArrayElem_SS(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                MIN_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }

        // ###################
        // Second
        // ###################
        elemName = "Second";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(5);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"SEC")]);}

            SEC_json = entryVal;

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(5);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"SEC")];}

            validEntryFlag = getArrayElem_SS(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                SEC_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }

        // chemname_ssi -> already obtained above // chemical name

        // ###################
        // ix
        // ###################
        // Always read as string (ID-based lookup).
        // Numeric indices are NOT supported -all spatial elements
        // are treated as text IDs to avoid index/ID mismatches.
        elemName = "ix";
        {
            std::string entryVal = ""; // dummy variable
            // if JSON -read as string regardless of JSON type
            if (DataFormat.compare("JSON")==0){
                // Handle both string and numeric JSON values by converting to string
                try {
                    entryVal = EWF_SS_json_sub_rowi.at(6).get<std::string>();
                } catch (...) {
                    // If it's a number in JSON, convert to string
                    try {
                        int numVal = EWF_SS_json_sub_rowi.at(6);
                        entryVal = std::to_string(numVal);
                    } catch (...) {
                        try {
                            double numVal = EWF_SS_json_sub_rowi.at(6);
                            entryVal = std::to_string(static_cast<int>(numVal));
                        } catch (...) {
                            entryVal = "";
                        }
                    }
                }
            }
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IX")];}

            // Strip surrounding quotes (CSV may preserve them)
            if (entryVal.size() >= 2
                && (entryVal.front() == '"' || entryVal.front() == '\'')
                && (entryVal.back() == '"' || entryVal.back() == '\'')) {
                entryVal = entryVal.substr(1, entryVal.size() - 2);
            }

            // Convert to uppercase for "ALL" comparison
            std::string entryVal_upper = OpenWQ_utils.ConvertStringToUpperCase(entryVal);

            // First, check if it's "ALL"
            if (entryVal_upper.compare("ALL") == 0) {
                ix_json = OpenWQ_wqconfig.get_allSS_flag();
                cell_id_lookup_used = false;  // Still need to parse iy and iz
            }
            // Otherwise, interpret as a cell_id (e.g., reach_id, hru_id)
            else {
                bool partial_match = false;
                // Cell_id lookup - this will set ix_json, iy_json, iz_json if successful
                validEntryFlag = lookupCellId_SS(
                    OpenWQ_hostModelconfig,
                    OpenWQ_wqconfig,
                    OpenWQ_output,
                    static_cast<int>(cmpi_ssi),  // compartment index
                    entryVal,                     // cell_id string
                    ix_json, iy_json, iz_json,   // output indices
                    ssf, ssi, di,                // for error reporting
                    partial_match);              // partial prefix match flag

                if (!validEntryFlag) {
                    continue;  // Skip this entry if cell_id not found
                }
                // If exact match: skip iy/iz parsing (all indices set)
                // If partial/prefix match (e.g., "1" matched "1_z1"):
                //   only ix is reliable, iy/iz should still be parsed from JSON
                cell_id_lookup_used = !partial_match;
            }
        }

        // ###################
        // iy
        // ###################
        // Skip if cell_id lookup was used (iy already set)
        if (cell_id_lookup_used) {
            // iy_json already set by lookupCellId_SS
        } else {
            // Always read as string (ID-based lookup)
            elemName = "iy";
            std::string entryVal_iy = "";
            // if JSON
            if (DataFormat.compare("JSON")==0){
                try {
                    entryVal_iy = EWF_SS_json_sub_rowi.at(7).get<std::string>();
                } catch (...) {
                    try {
                        int numVal = EWF_SS_json_sub_rowi.at(7);
                        entryVal_iy = std::to_string(numVal);
                    } catch (...) {
                        try {
                            double numVal = EWF_SS_json_sub_rowi.at(7);
                            entryVal_iy = std::to_string(static_cast<int>(numVal));
                        } catch (...) {
                            entryVal_iy = "";
                        }
                    }
                }
            }
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal_iy = ASCIIRowElemEntry[
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IY")];}

            // Strip surrounding quotes (CSV may preserve them)
            if (entryVal_iy.size() >= 2
                && (entryVal_iy.front() == '"' || entryVal_iy.front() == '\'')
                && (entryVal_iy.back() == '"' || entryVal_iy.back() == '\'')) {
                entryVal_iy = entryVal_iy.substr(1, entryVal_iy.size() - 2);
            }

            // Convert to uppercase for "ALL" comparison
            std::string entryVal_iy_upper = OpenWQ_utils.ConvertStringToUpperCase(entryVal_iy);

            if (entryVal_iy_upper.compare("ALL") == 0) {
                iy_json = OpenWQ_wqconfig.get_allSS_flag();
            } else {
                // Try ID-based lookup along iy dimension
                bool found_iy = OpenWQ_hostModelconfig.find_index_single_dim(
                    static_cast<int>(cmpi_ssi), ix_json, 1, entryVal_iy, iy_json);
                if (!found_iy) {
                    msg_string =
                        "<OpenWQ> WARNING: SS iy='" + entryVal_iy +
                        "' not found in cell_id mapping for compartment " +
                        std::to_string(cmpi_ssi) +
                        " (entry skipped): File=" + std::to_string(ssf+1) +
                        ", Sub_structure=" + std::to_string(ssi+1) +
                        ", Data_row=" + std::to_string(di + 1);
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    continue;
                }
            }
        }  // end of else block for cell_id_lookup_used (iy)

        // ###################
        // iz
        // ###################
        // Skip if cell_id lookup was used (iz already set)
        if (cell_id_lookup_used) {
            // iz_json already set by lookupCellId_SS
        } else {
            // Always read as string (ID-based lookup)
            elemName = "iz";
            std::string entryVal_iz = "";
            // if JSON
            if (DataFormat.compare("JSON")==0){
                try {
                    entryVal_iz = EWF_SS_json_sub_rowi.at(8).get<std::string>();
                } catch (...) {
                    try {
                        int numVal = EWF_SS_json_sub_rowi.at(8);
                        entryVal_iz = std::to_string(numVal);
                    } catch (...) {
                        try {
                            double numVal = EWF_SS_json_sub_rowi.at(8);
                            entryVal_iz = std::to_string(static_cast<int>(numVal));
                        } catch (...) {
                            entryVal_iz = "";
                        }
                    }
                }
            }
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal_iz = ASCIIRowElemEntry[
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IZ")];}

            // Strip surrounding quotes (CSV may preserve them)
            if (entryVal_iz.size() >= 2
                && (entryVal_iz.front() == '"' || entryVal_iz.front() == '\'')
                && (entryVal_iz.back() == '"' || entryVal_iz.back() == '\'')) {
                entryVal_iz = entryVal_iz.substr(1, entryVal_iz.size() - 2);
            }

            // Convert to uppercase for "ALL" comparison
            std::string entryVal_iz_upper = OpenWQ_utils.ConvertStringToUpperCase(entryVal_iz);

            if (entryVal_iz_upper.compare("ALL") == 0) {
                iz_json = OpenWQ_wqconfig.get_allSS_flag();
            } else {
                // Try ID-based lookup along iz dimension
                bool found_iz = OpenWQ_hostModelconfig.find_index_single_dim(
                    static_cast<int>(cmpi_ssi), ix_json, 2, entryVal_iz, iz_json);
                if (!found_iz) {
                    msg_string =
                        "<OpenWQ> WARNING: SS iz='" + entryVal_iz +
                        "' not found in cell_id mapping for compartment " +
                        std::to_string(cmpi_ssi) +
                        " (entry skipped): File=" + std::to_string(ssf+1) +
                        ", Sub_structure=" + std::to_string(ssi+1) +
                        ", Data_row=" + std::to_string(di + 1);
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    continue;
                }
            }
        }  // end of else block for cell_id_lookup_used (iz)

        // Reset cell_id_lookup_used for next iteration
        cell_id_lookup_used = false;

        // ###################
        // SS sink/source load or EWF conc
        // cannot have negative values
        // ###################

        double entryVal = 0.0f; // dummy variable
        // if JSON
        if (DataFormat.compare("JSON")==0){
            entryVal = EWF_SS_json_sub_rowi.at(9);}
        // if ASCII
        else if (DataFormat.compare("ASCII")==0){
            entryVal = std::stod(ASCIIRowElemEntry[ 
            OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"LOAD")]);}

        ss_data_json = entryVal;

        // skip if negative value in SS load/sink or EWF conc 
        // throw warning msg
        if(ss_data_json < 0.0f){
            // Create Warning Message
            msg_string = 
                    "<OpenWQ> WARNING: SS/EWF '" 
                    " load/sink/conc cannot be negative (entry skipped): File=" 
                    + std::to_string(ssf+1)
                    + ", Sub_structure=" + std::to_string(ssi+1)
                    + ", Data_row=" + std::to_string(di + 1);  
            // Print it (Console and/or Log file)
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            validEntryFlag = false;
            if (!validEntryFlag){continue;}
        }

        // ###################
        // Load scheme (discrete, continuous)
        // ###################

        elemName = "Load/Sink Scheme";

        ss_units_json = ss_units_json_mass_base;
        
        try{

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(10);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"LOAD_TYPE")];}

            loadScheme_str = entryVal;

            // Strip surrounding quotes (CSV may preserve them)
            if (loadScheme_str.size() >= 2
                && (loadScheme_str.front() == '"' || loadScheme_str.front() == '\'')
                && (loadScheme_str.back() == '"' || loadScheme_str.back() == '\'')) {
                loadScheme_str = loadScheme_str.substr(1, loadScheme_str.size() - 2);
            }

            // Convert to uppercase for case-insensitive comparison
            std::transform(loadScheme_str.begin(), loadScheme_str.end(),
                          loadScheme_str.begin(), ::toupper);

            // loading scheme only needed if SS
            // EWF is in concentration and associated with fluxes
            if (inputType.compare("ss")==0){

                // Set loadScheme_id
                // 1) discrete
                // 2) continuous (needs time units)
                if (loadScheme_str.compare("DISCRETE") == 0) loadScheme_id = 0;
                else if (loadScheme_str.compare("CONTINUOUS") == 0 && SEC_json != -1){ 
                    // continuous option needs SEC = 'all' to allow a minimum continuous period
                    loadScheme_id = 0;
                    // Create Warning Message
                    msg_string = 
                            "<OpenWQ> WARNING: SS/EWF '" 
                            + elemName 
                            + "' was defaulted to 'discrete'. The 'continuous' option is only valid with SEC set as 'ALL' to allow a minimum continuous load period (otherwise, it becomes a discrete load): File=" 
                            + std::to_string(ssf+1)
                            + ", Sub_structure=" + std::to_string(ssi+1)
                            + ", Data_row=" + std::to_string(di + 1);  
                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                }else if (loadScheme_str.compare("CONTINUOUS") == 0 && SEC_json == -1){
                    // continuous option needs SEC = 'all' (otherwise it's discrete input)
                    loadScheme_id = 1;
                    // get time units
                    try{

                        std::string entryVal = ""; // dummy variable
                        // if JSON
                        if (DataFormat.compare("JSON")==0){
                            entryVal = EWF_SS_json_sub_rowi.at(11);}
                        // if ASCII
                        else if (DataFormat.compare("ASCII")==0){
                            entryVal = ASCIIRowElemEntry[ 
                                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"TIME_UNITS")];}

                        contDt_str = entryVal;

                        // Strip surrounding quotes (CSV may preserve them)
                        if (contDt_str.size() >= 2
                            && (contDt_str.front() == '"' || contDt_str.front() == '\'')
                            && (contDt_str.back() == '"' || contDt_str.back() == '\'')) {
                            contDt_str = contDt_str.substr(1, contDt_str.size() - 2);
                        }

                        // Convert to uppercase for case-insensitive comparison
                        std::transform(contDt_str.begin(), contDt_str.end(),
                                      contDt_str.begin(), ::toupper);

                        // Concatenate the time units to the load
                        ss_units_json += "/";
                        ss_units_json += contDt_str;

                    }catch(...){ 
                        // Create Warning Message
                        msg_string = 
                            "<OpenWQ> WARNING: SS/EWF '" 
                            + elemName 
                            + "'continuous' requires an additional array element with the load/sink time units (entry skipped): File=" 
                            + std::to_string(ssf+1)
                            + ", Sub_structure=" + std::to_string(ssi+1)
                            + ", Data_row=" + std::to_string(di + 1); 
                        // Print it (Console and/or Log file)
                        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                        continue; // skip entry
                    }
                }else{

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                    continue;
                }

            }else{
                // if EWF, send a warning message saying that the 
                // load scheme and period are not used in EWF inputs

                // Create Warning Message
                msg_string = 
                    "<OpenWQ> WARNING: EWF '" 
                    + elemName 
                    + "' is not used in EWF entries because these are" 
                    " concentrations associated with external fluxes (entry ignored): File=" 
                    + std::to_string(ssf+1)
                    + ", Sub_structure=" + std::to_string(ssi+1)
                    + ", Data_row=" + std::to_string(di + 1); 

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

            }

        }catch(...){  
            
            // load scheme is needed for SS, but not for EWF
            if (inputType.compare("ss")==0){

                // Create Warning Message
                // only printed if entry is not valid
                msg_string = 
                    "<OpenWQ> WARNING: SS '" 
                    + elemName 
                    + "' is not valid. It can only be 'discrete' or 'continuous' (entry skipped): File=" 
                    + std::to_string(ssf+1)
                    + ", Sub_structure=" + std::to_string(ssi+1)
                    + ", Data_row=" + std::to_string(di + 1); 

                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                continue;
            }
        }

        // Convert SS units
        // Source/sink units (g -> default model mass units)
        // 1) Calculate unit multiplers
        std::vector<std::string> units;          // units (numerator and denominator)
        OpenWQ_units.Calc_Unit_Multipliers(
            OpenWQ_wqconfig,
            OpenWQ_output,
            unit_multiplers,    // multiplers (numerator and denominator)
            ss_units_json,      // input units
            units,
            true);              // direction of the conversion: 
                                // to native (true) or 
                                // from native to desired output units (false)

        // 2) Calculate value with new units
        OpenWQ_units.Convert_Units(
            ss_data_json,       // value passed by reference so that it can be changed
            unit_multiplers);   // units

        // Get the vector with the data
        row_data_col = {
            (double) chem_ssi,
            (double) cmpi_ssi,
            (double) sinksource_ssi,
            (double) YYYY_json,
            (double) MM_json,
            (double) DD_json,
            (double) HH_json,
            (double) MIN_json,
            (double) SEC_json,
            (double) ix_json,
            (double) iy_json,
            (double) iz_json,
            ss_data_json,
            loadScheme_id,  // only for SS: load scheme (0) not applicable, (1) discrete or (2) continuous
            0,0,0,0,0,0     // field to specify the number of times it has been used aleady
            };              // in the case of and "all" element (YYYY, MM, DD, HH, MIN, SEC)
                            // it starts with 0 (zero), meaning that has not been used
                            // if not an "all" element, then it's set to -1

        // Add new row to SinkSource_FORC or ExtFlux_FORC_jsonAscii
        AppendRow_SS_EWF_FORC_jsonAscii(
            OpenWQ_wqconfig,
            inputType,
            row_data_col);

    }

 }

/* #################################################
 // Prepare SS and EWF input data for use at running time: 
 // Case: JSON and ASCII format
 ################################################# */
void OpenWQ_extwatflux_ss::Set_EWF_h5(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output,
    json EWF_SS_json_sub,  // relevant sub-json
    std::string inputType,
    bool foundflag){
    
    // Local variables
    std::string msg_string;                 // console/logFile message
    std::size_t it;                         // interactor
    std::string ewf_h5_folderPath;          // json input: folder with EWF input *.h5 files
    std::string ewf_filenamePath;           // full path for *.h5 files
    std::string ewf_h5_units;               // json input: units of EWF h5 files
    std::string ewf_h5_units_file;          // Same as above, but "\" replaced by "|" for file search
    double conc_h5_rowi;                    // iteractive concentration extracted for h5 files
    std::vector<std::string> units;         // Vector with units information
    std::string external_compartName;       // json input: External compartment name
    std::string external_waterFluxName;     // EWF name (needs to match name hard coded in hydrolink)
    double external_waterFluxName_id;       // index of external_waterFluxName
    std::string chemname;                   // iteractive chemical name during h5 search
    arma::mat xyzEWF_h5;                    // h5 xyz field 
    arma::mat domain_EWF_h5;                // h5 information about external compartment domain (nx,ny,nz)
    arma::mat dataEWF_h5;                   // h5 data
    std::vector<double> unit_multiplers;    // multiplers (numerator and denominator)
    std::vector<std::string> tSamp_valid;   // save valid time stamps
    time_t tSamp_valid_i_time_t;            // iteractive time_t from h5 timestamp
    int x_externModel, y_externModel, z_externModel;        // iteractive x,y,z info from h5 files
    int x_interface_h5, y_interface_h5, z_interface_h5;     // iteractive x,y,z from interface requested in h5 ewf files
    int nx_interface_h5, ny_interface_h5, nz_interface_h5;  // nx,ny,nz of external compartment from h5 ewf file
    int ewfName_nx, ewfName_ny, ewfName_nz;                 // nx, ny and nz of EWF associated
    std::string ss_cmp_recipient_name;                      // name of EWF recipient
    json interaction_interface_json;                        // json substructure for interface info
    int index_i;
    bool validEntryFlag;                    // flag for valid entries
    bool foundTimeStmps;                    // flag to record (un)success in finding timestamps
    bool h5_entry_found;                    // flag for successful finding of ewf h5 file
    std::string errorMsgIdentifier;         // Start/head of error message of json key not found
    std::vector<int> valid_interfaceH5rows; // vector with indexes of relevant h5 rows that contain interface data
    int rowi_val;                           // iteractive row number from valid_interfaceH5rows 
    int point_print_n;                      // iterative trackking of "." console prints (each timeStamp) for asthetics
    bool flag_newJSON_h5Request = true;      // flag for new json block for ewf-h5
    bool flag_newChem = true;               // flag for new chem from json ewf-h5 clock
    int h5EWF_request_index;                // Index of ewf-h5 index

    // Get request index
    h5EWF_request_index = (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time).size();

    // ################################
    // Get JSON info
    // ################################

    errorMsgIdentifier = inputType + " json block with DataFormat=HDF5" ;

    // h5 IC folder path
    ewf_h5_folderPath = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "FOLDERPATH",
        errorMsgIdentifier,
        true);
    // h5 ic units
    ewf_h5_units = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "UNITS",
        errorMsgIdentifier,
        true);
    // get external compartment name (needed for both ss and ewf)
    external_compartName = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "EXTERNAL_COMPARTMENT_NAME",
        errorMsgIdentifier,
        true);
    // Get spatial coupling mode (LUMPED or DISTRIBUTED)
    //  - LUMPED:      the EWF h5 carries a single column (one HRU/reach); its
    //                 concentration is broadcast to every receiving reach.
    //  - DISTRIBUTED: the EWF h5 carries one column per receiving reach plus a
    //                 numeric 'cellID' dataset; columns are matched to reaches
    //                 by id (any mismatch => console error + abort).
    // (Replaces the former INTERACTION_INTERFACE dimension-matching, which
    //  assumed the external model shared this model's discretization.)
    std::string spatial_mode = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "SPATIAL_MODE",
        errorMsgIdentifier,
        true);
    // Get external flux name
    external_waterFluxName = OpenWQ_utils.RequestJsonKeyVal_str(
        OpenWQ_wqconfig, OpenWQ_output,
        EWF_SS_json_sub, "EXTERNAL_INPUTFLUX_NAME",
        errorMsgIdentifier,
        true);
    // Get interpolation method
    OpenWQ_wqconfig.set_h5EWF_interpMethod(
                OpenWQ_utils.RequestJsonKeyVal_str(
                    OpenWQ_wqconfig, OpenWQ_output,
                    EWF_SS_json_sub, "INTERPOLATION",
                    errorMsgIdentifier,
                    true)) ;

    // ################################
    // Some pre-processing
    // ################################

    // replace "/" by "|" is needed because "/" is not compatible with directory full paths
    ewf_h5_units_file = ewf_h5_units;
    it = (int) ewf_h5_units_file.find("/");
    if (it <= ewf_h5_units_file.size()){
        ewf_h5_units_file.replace(it,1, "|");
    }

    // Get unit conversion multipliers (return value unused; call populates unit_multiplers)
    OpenWQ_units.Calc_Unit_Multipliers(
        OpenWQ_wqconfig,
        OpenWQ_output,
        unit_multiplers,    // multiplers (numerator and denominator)
        ewf_h5_units,       // input units
        units,
        true);              // direction of the conversion: 
                            // to native (true) or 
                            // from native to desired output units (false)

    // Get corresponding id
    external_waterFluxName_id = 
        (double)OpenWQ_utils.FindStrIndexInVectStr(
            OpenWQ_hostModelconfig.get_HydroExtFlux_names(),
            external_waterFluxName);

    // If external compartment not found in internal list of EWF, then
    // throw warning msg and skip entry
    // Otherwise save it in (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_ewfCompID)
    if(external_waterFluxName_id==-1.0f){
        msg_string = 
            "<OpenWQ> WARNNING SS json key EXTERNAL_INPUTFLUX_NAME= "
            + external_waterFluxName
            + " not valid for this host-model coupling (entry ignored).";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        return;
    }else{
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_ewfCompID).push_back(external_waterFluxName_id);
        // Open a (dense) chem-id list for this request. Each species that is
        // actually found+loaded below appends its true BGC id here, so the
        // dense storage index maps back to the correct species at apply time.
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_chemID).push_back(std::vector<int>{});
    }

    // Get num of interface elements
    ewfName_nx = OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_x_at(external_waterFluxName_id);
    ewfName_ny = OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_y_at(external_waterFluxName_id);
    ewfName_nz = OpenWQ_hostModelconfig.get_HydroExtFlux_num_cells_z_at(external_waterFluxName_id);
    
    // Generate arma::cube of compartment ewfi size
    // And reset ExtFlux_FORC_data_tStep for dimensions of 
    // ewf of index external_waterFluxName_id
    arma::Cube<double> ewfi_domain_xyz(ewfName_nx, ewfName_ny, ewfName_nz);
    (*OpenWQ_wqconfig.ExtFlux_FORC_data_tStep) = ewfi_domain_xyz;

    // Get valid time steps from the logFile of EWF simulation
    // if timeStamps not found, then return 
    // warning message alerady printed in GetTimeStampsFromLogFile
    foundTimeStmps = OpenWQ_utils.GetTimeStampsFromLogFile(
        OpenWQ_wqconfig,
        OpenWQ_output,
        ewf_h5_folderPath,
        "<OpenWQ> Output export successful (HDF5): ", // Substring of the output to search
        tSamp_valid,
        "SS/EWF load/sink/conc H5 supporting logFile"); // Logfile errMsg identifier
    if (foundTimeStmps==false) return;

    // Throw console message to say it's processing the h5 interface data
    msg_string = "<OpenWQ> EWF HDF5 interface requested.\n"
                 "         Processing and checking interface...";
    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, false);


    // ################################
    // Loop over EWF h5 files
    // ################################
    // Each chemical species is in different files

    unsigned int num_chem;
    if ((OpenWQ_wqconfig.CH_model->BGC_module).compare("NATIVE_BGC_FLEX") == 0) {
        num_chem = OpenWQ_wqconfig.CH_model->NativeFlex->num_chem;
    } else {
        num_chem = OpenWQ_wqconfig.CH_model->PHREEQC->num_chem;
    }
    // ---------------------------------------------------------------
    // Recipient compartment index for DISTRIBUTED id-matching.
    // Reach ids live in the RIVER_NETWORK_REACHES compartment (mizuRoute);
    // fall back to compartment 0 if that name is not present.
    // ---------------------------------------------------------------
    int recipient_comp_index = 0;
    {
        std::vector<std::string> _compNames = OpenWQ_hostModelconfig.get_HydroComp_names();
        for (unsigned int _ci = 0; _ci < _compNames.size(); _ci++){
            if (_compNames[_ci].compare("RIVER_NETWORK_REACHES") == 0){
                recipient_comp_index = (int)_ci; break;
            }
        }
    }

    // Number of output timesteps (from the EWF source logFile)
    const long unsigned int nTsteps_h5 = tSamp_valid.size();

    // Normalize SPATIAL_MODE
    std::string spatial_mode_uc = spatial_mode;
    std::transform(spatial_mode_uc.begin(), spatial_mode_uc.end(),
                   spatial_mode_uc.begin(), ::toupper);
    const bool mode_lumped      = (spatial_mode_uc.compare("LUMPED") == 0);
    const bool mode_distributed = (spatial_mode_uc.compare("DISTRIBUTED") == 0);
    if (!mode_lumped && !mode_distributed){
        msg_string =
            "<OpenWQ> ERROR: EWF '" + external_waterFluxName
            + "' has unknown SPATIAL_MODE='" + spatial_mode
            + "' (expected 'LUMPED' or 'DISTRIBUTED'). Aborting.";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
        std::exit(EXIT_FAILURE);
    }

    // ################################
    // Loop over EWF h5 files (one per chemical species)
    // ################################
    // Dense storage-slot counter: increments once per species actually loaded.
    // A species with no source h5 is skipped, so this decouples the storage
    // index from the BGC species index (mapping recorded in _chemID).
    int dense_chem_idx = 0;
    for (unsigned int chemi = 0; chemi < num_chem; chemi++){

        // Set new chem flag true
        flag_newChem = true;

        // Get chem name
        chemname = (*OpenWQ_wqconfig.cached_chem_species_list_ptr)[chemi];

        // Console update
        msg_string = "         " + external_waterFluxName + " => " + chemname
                     + " [" + spatial_mode_uc + "] ";
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, false);

        // Build h5 filename: <folder>/<EXTERNAL_COMPARTMENT_NAME>@<chem>#<units>-main.h5
        ewf_filenamePath = ewf_h5_folderPath;
        ewf_filenamePath.append("/");
        ewf_filenamePath.append(external_compartName);
        ewf_filenamePath.append("@");
        ewf_filenamePath.append(chemname);
        ewf_filenamePath.append("#");
        ewf_filenamePath.append(ewf_h5_units_file);
        ewf_filenamePath.append("-main.h5");

        // ------------------------------------------------------------
        // Load the consolidated concentration matrix (/concentrations),
        // written by the openWQ output writer as a 2D dataset [time x cells].
        // Armadillo may transpose on load, so orientation is resolved below
        // using the known number of timesteps (nTsteps_h5).
        // ------------------------------------------------------------
        arma::mat conc_h5;
        conc_h5.load(arma::hdf5_name(ewf_filenamePath, "concentrations"));

        if (conc_h5.is_empty()){
            msg_string =
                "<OpenWQ> WARNING: EWF h5 file requested=" + ewf_filenamePath
                + " was not found or has no 'concentrations' dataset "
                "(entry skipped).";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            continue;
        }

        // Species found: record its true BGC id for this dense storage slot.
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_chemID).back().push_back((int)chemi);

        // Resolve orientation: which dimension is time vs cells
        unsigned int ncells_src;
        bool time_is_rows;
        if (conc_h5.n_rows == nTsteps_h5){
            time_is_rows = true;  ncells_src = conc_h5.n_cols;   // (time x cells)
        } else if (conc_h5.n_cols == nTsteps_h5){
            time_is_rows = false; ncells_src = conc_h5.n_rows;   // (cells x time)
        } else {
            msg_string =
                "<OpenWQ> ERROR: EWF '" + external_waterFluxName
                + "' h5 file=" + ewf_filenamePath
                + " has a 'concentrations' matrix (" + std::to_string(conc_h5.n_rows)
                + "x" + std::to_string(conc_h5.n_cols)
                + ") whose dimensions do not match the number of timesteps ("
                + std::to_string(nTsteps_h5) + ") in the supporting logFile. Aborting.";
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            std::exit(EXIT_FAILURE);
        }

        // ------------------------------------------------------------
        // Map each source column -> receiving-EWF cell index (0-based).
        //  - LUMPED:      require exactly 1 source column; broadcast to all cells.
        //  - DISTRIBUTED: require #columns == #EWF cells AND that every source
        //                 cell_id matches a host reach id; abort otherwise.
        // ------------------------------------------------------------
        std::vector<int> col2ewfIx(ncells_src, -1);

        if (mode_lumped){
            // LUMPED broadcasts a SINGLE source concentration to every reach.
            // If the source h5 has MORE than one column (e.g. the SUMMA soil
            // compartment carries one column per soil layer), collapse it to a
            // single representative series = the per-timestep MEAN across the
            // source cells (ignoring noWaterConc sentinels). This lets a
            // multi-cell compartment (soil layers, etc.) be used as a lumped EWF
            // source without a separate remap step. ncells_src==1 is unchanged.
            if (ncells_src > 1){
                arma::mat collapsed;
                if (time_is_rows){          // [T x C] -> [T x 1]
                    collapsed.set_size(conc_h5.n_rows, 1);
                    for (arma::uword r = 0; r < conc_h5.n_rows; r++){
                        arma::rowvec row = conc_h5.row(r);
                        arma::uvec ok = arma::find(row > -9990.0);
                        collapsed(r, 0) = ok.n_elem ? arma::mean(row.elem(ok)) : -9999.0;
                    }
                } else {                    // [C x T] -> [1 x T]
                    collapsed.set_size(1, conc_h5.n_cols);
                    for (arma::uword c = 0; c < conc_h5.n_cols; c++){
                        arma::vec col = conc_h5.col(c);
                        arma::uvec ok = arma::find(col > -9990.0);
                        collapsed(0, c) = ok.n_elem ? arma::mean(col.elem(ok)) : -9999.0;
                    }
                }
                msg_string =
                    "<OpenWQ> NOTE: EWF '" + external_waterFluxName
                    + "' SPATIAL_MODE=LUMPED source " + ewf_filenamePath + " has "
                    + std::to_string(ncells_src)
                    + " columns; collapsed to their per-timestep MEAN (broadcast to all reaches).";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, false);
                conc_h5 = collapsed;
                ncells_src = 1;
            }
            // broadcast handled in the timestep loop below
        }
        else { // mode_distributed
            if ((int)ncells_src != ewfName_nx){
                msg_string =
                    "<OpenWQ> ERROR: EWF '" + external_waterFluxName
                    + "' SPATIAL_MODE=DISTRIBUTED requires the h5 to have as many "
                    "columns as receiving reaches (" + std::to_string(ewfName_nx)
                    + "), but " + ewf_filenamePath + " has " + std::to_string(ncells_src)
                    + " columns. Provide a distributed source with matching reach ids. "
                    "Aborting.";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                std::exit(EXIT_FAILURE);
            }
            // Read the numeric per-column reach id vector (dataset 'cellID')
            arma::mat cellID_h5;
            cellID_h5.load(arma::hdf5_name(ewf_filenamePath, "cellID"));
            if (cellID_h5.n_elem != ncells_src){
                msg_string =
                    "<OpenWQ> ERROR: EWF '" + external_waterFluxName
                    + "' SPATIAL_MODE=DISTRIBUTED requires a numeric 'cellID' dataset "
                    "with " + std::to_string(ncells_src) + " reach ids in "
                    + ewf_filenamePath + " (found " + std::to_string(cellID_h5.n_elem)
                    + "). Aborting.";
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                std::exit(EXIT_FAILURE);
            }
            for (unsigned int j = 0; j < ncells_src; j++){
                long long _idv = (long long) std::llround(cellID_h5(j));
                std::string _id_str = std::to_string(_idv);
                int _ix = -1, _iy = -1, _iz = -1; bool _partial = false;
                bool _found = OpenWQ_hostModelconfig.find_indices_from_cellid(
                    recipient_comp_index, _id_str, _ix, _iy, _iz, _partial);
                if (!_found){
                    msg_string =
                        "<OpenWQ> ERROR: EWF '" + external_waterFluxName
                        + "' SPATIAL_MODE=DISTRIBUTED: reach id '" + _id_str
                        + "' from " + ewf_filenamePath
                        + " does not match any host reach id. Aborting.";
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                    std::exit(EXIT_FAILURE);
                }
                col2ewfIx[j] = _ix;
            }
        }

        // ------------------------------------------------------------
        // Fill ExtFlux_FORC_data_tStep per timestep and append to the
        // runtime store (one cube per timestep, per chem).
        // ------------------------------------------------------------
        int point_print_n2 = 0;
        for (long unsigned int t = 0; t < nTsteps_h5; t++){

            // timestamp string -> time_t
            tSamp_valid_i_time_t = OpenWQ_units.convertTime_str2time_t(
                OpenWQ_wqconfig, tSamp_valid[t]);

            // reset scratch cube
            (*OpenWQ_wqconfig.ExtFlux_FORC_data_tStep).zeros();

            if (mode_lumped){
                double _c = time_is_rows ? conc_h5(t, 0) : conc_h5(0, t);
                if (_c < -9990.0) _c = 0.0;          // noWaterConc sentinel -> 0 solute
                else OpenWQ_units.Convert_Units(_c, unit_multiplers);
                for (int ix = 0; ix < ewfName_nx; ix++)
                    (*OpenWQ_wqconfig.ExtFlux_FORC_data_tStep)(ix, 0, 0) = _c;
            }
            else { // mode_distributed
                for (unsigned int j = 0; j < ncells_src; j++){
                    double _c = time_is_rows ? conc_h5(t, j) : conc_h5(j, t);
                    if (_c < -9990.0) _c = 0.0;
                    else OpenWQ_units.Convert_Units(_c, unit_multiplers);
                    (*OpenWQ_wqconfig.ExtFlux_FORC_data_tStep)(col2ewfIx[j], 0, 0) = _c;
                }
            }

            AppendCube_SS_EWF_FORC_h5(
                OpenWQ_wqconfig,
                h5EWF_request_index,
                dense_chem_idx,   // dense storage slot (NOT the BGC id — mapped via _chemID)
                flag_newChem,
                flag_newJSON_h5Request,
                tSamp_valid_i_time_t);

            flag_newJSON_h5Request = false;
            flag_newChem = false;

            // progress dots (one line per 80 timesteps)
            point_print_n2++;
            if (point_print_n2 == 80){ point_print_n2 = 0; std::cout << "\n         " << std::flush; }
            std::cout << "." << std::flush;
        }

        std::cout << " => TimeSteps processed: "
                  + std::to_string(nTsteps_h5) + "\n" << std::flush;

        // Advance the dense storage slot (this species was loaded successfully)
        dense_chem_idx++;
    }
}

/* #################################################
 // Get model structure indexes for compartments and chemicals
 ################################################# */
bool OpenWQ_extwatflux_ss::getModIndex(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    std::vector<std::string> &vec_list,
    std::string &obj_name,
    std::string &obj_text,
    unsigned long &vec_obj_index){
    
    // Local Variables
    bool foundflag = false;
    std::vector<std::string>::iterator find_i;  // iteractor used to store the position or searched element
    std::string msg_string;             // error/warning message string

    // Try to find index
    find_i = 
        std::find(vec_list.begin(), 
        vec_list.end(), 
        obj_name);

    // If requested index exists, then okay 
    // (otherwise, throw warning and skip entry)

    if (find_i != vec_list.end()){
        vec_obj_index =   find_i - vec_list.begin();
        foundflag = true;
    }else{
        
        // Create Message (WARNING: entry skipped)
        msg_string = 
            "<OpenWQ> WARNNING (entry skipped): " 
            + obj_text 
            + " in source-sink file unkown: " 
            + obj_name;

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

    }
    return foundflag;

}

/* #################################################
 // Get SS vector element
 // CASE IF: elemEntry as string "all"
 ################################################# */
 bool OpenWQ_extwatflux_ss::getArrayElem_SS(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    std::string elemName,
    std::__cxx11::basic_string<char> elemEntry,
    int& elemVal,
    unsigned int& file_i,
    unsigned int& struc_i,
    unsigned int& row_i){

    // Local Variable
    bool validEntryFlag = true;
    std::string msg_string;

    // Strip surrounding quotes (CSV may preserve them)
    if (elemEntry.size() >= 2
        && (elemEntry.front() == '"' || elemEntry.front() == '\'')
        && (elemEntry.back() == '"' || elemEntry.back() == '\'')) {
        elemEntry = elemEntry.substr(1, elemEntry.size() - 2);
    }

    // Convert to uppercase for case-insensitive comparison
    std::string elemEntry_upper = elemEntry;
    std::transform(elemEntry_upper.begin(), elemEntry_upper.end(), elemEntry_upper.begin(), ::toupper);

    if(elemEntry_upper.compare("ALL") == 0){

        elemVal = OpenWQ_wqconfig.get_allSS_flag();

    }else{

        // Through a warning invalid entry
        msg_string =
            "<OpenWQ> WARNING: SS '"
                            + elemName
                            + "' entry '" + elemEntry + "' is not a valid value"
                            + " (expected 'all' or valid ID). Entry skipped: File="
                            + std::to_string(file_i+1)
                            + ", Sub_structure=" + std::to_string(struc_i+1)
                            + ", Data_row=" + std::to_string(row_i + 1);

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

        validEntryFlag = false;
    }

    return validEntryFlag;

}

/* #################################################
 // Lookup cell_id (e.g., reach_id, hru_id) and return corresponding ix, iy, iz indices
 // This allows users to specify cell IDs from the host model (e.g., "1200014181")
 // instead of internal OpenWQ indices (ix, iy, iz)
 ################################################# */
bool OpenWQ_extwatflux_ss::lookupCellId_SS(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    int compartment_index,
    const std::string& cell_id_str,
    int& ix, int& iy, int& iz,
    unsigned int& file_i,
    unsigned int& struc_i,
    unsigned int& row_i,
    bool& partial_match){

    // Local variables
    std::string msg_string;
    std::string cellid_label = OpenWQ_hostModelconfig.get_cellid_to_wqlabel();

    // Try to find the cell_id in the mapping (exact or prefix match)
    bool found = OpenWQ_hostModelconfig.find_indices_from_cellid(
        compartment_index, cell_id_str, ix, iy, iz, partial_match);

    if (!found) {
        // Cell ID not found - print warning
        msg_string =
            "<OpenWQ> WARNING: SS cell_id '" + cell_id_str +
            "' (from " + cellid_label + ") not found in compartment " +
            std::to_string(compartment_index) +
            " (entry skipped): File=" + std::to_string(file_i + 1) +
            ", Sub_structure=" + std::to_string(struc_i + 1) +
            ", Data_row=" + std::to_string(row_i + 1);

        OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

        return false;
    }

    return true;
}

// Add new row to SinkSource_FORC or ExtFlux_FORC_jsonAscii
void OpenWQ_extwatflux_ss::AppendRow_SS_EWF_FORC_jsonAscii(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    std::string inputType,
    arma::vec row_data_col){

    // OPTIMIZED (perf): stage the row in a flat row-major buffer instead of
    // calling arma::insert_rows per row. insert_rows reallocates and copies the
    // entire growing matrix on every call, making ingestion O(n^2) in the number
    // of rows (minutes for sink-source files with tens of thousands of loads).
    // The buffer is converted to the FORC matrix in a single allocation by
    // Finalize_FORC_jsonAscii() once all rows have been read.
    const unsigned int ncol = row_data_col.n_elem;

    if (inputType.compare("ss")==0){
        OpenWQ_wqconfig.SinkSource_FORC_buffer.insert(
            OpenWQ_wqconfig.SinkSource_FORC_buffer.end(),
            row_data_col.begin(), row_data_col.end());
        OpenWQ_wqconfig.SinkSource_FORC_buffer_ncol = ncol;
    }
    else if (inputType.compare("ewf")==0){
        OpenWQ_wqconfig.ExtFlux_FORC_jsonAscii_buffer.insert(
            OpenWQ_wqconfig.ExtFlux_FORC_jsonAscii_buffer.end(),
            row_data_col.begin(), row_data_col.end());
        OpenWQ_wqconfig.ExtFlux_FORC_jsonAscii_buffer_ncol = ncol;
    }
}

// OPTIMIZED (perf): build the FORC matrix from the staged flat buffer in a
// single allocation. Replaces the previous per-row arma::insert_rows growth.
void OpenWQ_extwatflux_ss::Finalize_FORC_jsonAscii(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    std::string inputType){

    // Select the relevant buffer/matrix for this inputType
    std::vector<double>* buf = nullptr;
    unsigned int ncol = 0;
    arma::Mat<double>* mat = nullptr;

    if (inputType.compare("ss")==0){
        buf  = &OpenWQ_wqconfig.SinkSource_FORC_buffer;
        ncol = OpenWQ_wqconfig.SinkSource_FORC_buffer_ncol;
        mat  = OpenWQ_wqconfig.SinkSource_FORC.get();
    }
    else if (inputType.compare("ewf")==0){
        buf  = &OpenWQ_wqconfig.ExtFlux_FORC_jsonAscii_buffer;
        ncol = OpenWQ_wqconfig.ExtFlux_FORC_jsonAscii_buffer_ncol;
        mat  = OpenWQ_wqconfig.ExtFlux_FORC_jsonAscii.get();
    }
    else { return; }

    // Nothing staged (e.g. EWF provided only via HDF5) -> leave matrix as-is
    if (buf->empty() || ncol == 0) return;

    const unsigned int nrow = buf->size() / ncol;

    // The buffer is row-major (ncol values per row). Armadillo is column-major,
    // so interpret the buffer as an (ncol x nrow) matrix (each column = one row)
    // and transpose once to obtain the desired (nrow x ncol) layout.
    arma::Mat<double> staged(buf->data(), ncol, nrow);  // copies aux memory
    arma::Mat<double> staged_t = staged.t();

    if (mat->n_rows == 0){
        *mat = std::move(staged_t);
    } else {
        // Preserve rows already present (defensive; normally the matrix is empty)
        *mat = arma::join_cols(*mat, staged_t);
    }

    // Release the staging buffer
    buf->clear();
    buf->shrink_to_fit();
}

// Add new row to SinkSource_FORC or ExtFlux_FORC_jsonAscii
void OpenWQ_extwatflux_ss::AppendCube_SS_EWF_FORC_h5(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    int h5EWF_request_index,        // get request index
    int chemi,                      // chem index
    bool flag_newChem,              // flag for new timestep, push back new vector row [i]
    bool flag_newJSON_h5Request,     // new json-h5-ewf request
    time_t timestamp_time_t){        // timestamp in time_t

    
    // Push_back/Create vector<vector<arma>> for every new request (ewf-h5) 
    if(flag_newJSON_h5Request==true){
        // Time
        std::vector<time_t> newEntryArma_time; 
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time).push_back(newEntryArma_time);
        // Data
        std::vector<std::vector<arma::Cube<double>>> newEntryArma_data; 
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data).push_back(newEntryArma_data);
    }

    // Push_back/Create vector<arma> for every new chem 
    if(flag_newChem==true){
        std::vector<arma::Cube<double>> newChemArma; // create vector<arma> for every new chem
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data)[h5EWF_request_index].push_back(newChemArma);

        
    }

    // Add new timestamp to ExtFlux_FORC_HDF5vec_time
    // But only needed on first chemi pass
    if(chemi==0){
        (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_time)[h5EWF_request_index].push_back(
            timestamp_time_t);
    }

    // add timestep data to vector ExtFlux_FORC_HDF5vec_data[chemi]
    (*OpenWQ_wqconfig.ExtFlux_FORC_HDF5vec_data)[h5EWF_request_index][chemi].push_back(
        *OpenWQ_wqconfig.ExtFlux_FORC_data_tStep);

    // Reset ExtFlux_FORC_data_tStep
    (*OpenWQ_wqconfig.ExtFlux_FORC_data_tStep).zeros();


}