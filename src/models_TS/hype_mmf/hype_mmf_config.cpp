// Copyright 2020, Diogo Costa, diogo.costa@uevora.pt
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

#include "models_TS/headerfile_TS.hpp"


 void OpenWQ_TS_model::TS_HYPE_MMF_config(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
            
        
    // Local variables
    /*
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
    unsigned int ascii_numHeaderRows;           // additional information for ASCII data input
    unsigned int ascii_headerKeyRow;            // additional information for ASCII data input
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
    */

    std::string DataFormat;

    DataFormat = std::get<3>(OpenWQ_wqconfig.TS_model->HypeMMF->info_vector);

    /* ########
    // Check if DATA_FORMAT=ASCII
    // if yes, then get file structure info
    ###########*/

    /*
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

            // number of header rows
            ascii_numHeaderRows = OpenWQ_utils.RequestJsonKeyVal_int(
                OpenWQ_wqconfig, OpenWQ_output,
                EWF_SS_json_sub["DATA"], "NUMBER_OF_HEADER_ROWS",
                errorMsgIdentifier,
                true);            

            // position of key header
            ascii_headerKeyRow = OpenWQ_utils.RequestJsonKeyVal_int(
                OpenWQ_wqconfig, OpenWQ_output,
                EWF_SS_json_sub["DATA"], "HEADER_KEY_ROW",
                errorMsgIdentifier,
            true);  

            // Get number of rows in ASCII file
            num_rowdata = OpenWQ_utils.getNumLinesfromASCIIFile(ascii_FilePath)
                - ascii_numHeaderRows;

            // Open ASCII file (read mode only)
            std::ifstream asciiFile (ascii_FilePath, std::ios::in);

            // Read skip header lines
            if (asciiFile.is_open()){
                lineCount = 1;
                while(std::getline(asciiFile, rowEntryASCII)){
                    // get header row when found
                    if(lineCount==ascii_headerKeyRow){
                        headerRowASCII = rowEntryASCII;
                        // Convert string to uppercase
                        headerRowASCII = OpenWQ_utils.ConvertStringToUpperCase(headerRowASCII);
                        // Parse the headerRowASCII to get a vector with the json entries
                        // parsing based on the delimiter ascii_delimiter defined by the user
                        headerKeys = OpenWQ_utils.StringSplit(
                            headerRowASCII,
                            ascii_delimiter);}
                    // push new entry to if actual row data
                    if(lineCount>ascii_numHeaderRows) ASCIIRowdata.push_back(rowEntryASCII);
                    // update line count
                    ++lineCount;
                }
                asciiFile.close();
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

    // ########################################
    // Loop over row data in sink-source file
    ######################################## //

    for (unsigned int di=0;di<num_rowdata;di++){
        
        // Reset the size to zero (the object will have no elements)
        row_data_col.reset(); 

        // Get row-json di from EWF_SS_json_sub ["DATA"]
        // Abort if not found
        errorMsgIdentifier = inputType + " json block 'DATA', row " 
                            + std::to_string(di) + " with DataFormat=" 
                            + DataFormat;
        EWF_SS_json_sub_rowi = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            EWF_SS_json_sub["DATA"], std::to_string(di+1),
            errorMsgIdentifier,
            true);

        // If DataFormat=ASCII, then get row data 
        // and convert-to-upper-case and split it by element entry
        if (DataFormat.compare("ASCII")==0){
            ASCIIRowElemEntry = 
                OpenWQ_utils.StringSplit(
                    OpenWQ_utils.ConvertStringToUpperCase(ASCIIRowdata[di]),
                    ascii_delimiter);}

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
            validEntryFlag = getArrayElem(
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

            validEntryFlag = getArrayElem(
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

            validEntryFlag = getArrayElem(
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

            validEntryFlag = getArrayElem(
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

            validEntryFlag = getArrayElem(
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

            validEntryFlag = getArrayElem(
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
        elemName = "ix";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(6);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IX")]);}

            ix_json = entryVal;

            // Need to do "- 1" because C++ starts in zero
            ix_json--;

            // If entry in json is zero, we get here -1
            // which is wrong, so need to send warning messsage
            // and skip entry
            if (ix_json == -1){
                // Through a warning invalid entry           
                msg_string = 
                    "<OpenWQ> WARNING: SS '" 
                    + elemName 
                    + "' cannot be zero. It needs to start in one (entry skipped): File=" 
                    + std::to_string(ssf+1)
                    + ", Sub_structure=" + std::to_string(ssi+1)
                    + ", Data_row=" + std::to_string(di + 1);   

                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                continue; // skip entry
            }

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(6);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IX")];}

            validEntryFlag = getArrayElem(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                ix_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }

        // ###################
        // iy
        // ###################
        elemName = "iy";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(7);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IY")]);}

            iy_json = entryVal;

            // Need to do "- 1" because C++ starts in zero
            iy_json--;

            // If entry in json is zero, we get here -1
            // which is wrong, so need to send warning messsage
            // and skip entry
            if (iy_json == -1){
                // Through a warning invalid entry           
                msg_string = 
                    "<OpenWQ> WARNING: SS '" 
                    + elemName 
                    + "' cannot be zero. It needs to start in one (entry skipped): File=" 
                    + std::to_string(ssf+1)
                    + ", Sub_structure=" + std::to_string(ssi+1)
                    + ", Data_row=" + std::to_string(di + 1);   

                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                continue; // skip entry
            }

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(7);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IY")];}

            validEntryFlag = getArrayElem(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                iy_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }

        // ###################
        // iz
        // ###################
        elemName = "iz";
        try{

            // try entry as int
            int entryVal = 0; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(8);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = std::stoi(ASCIIRowElemEntry[ 
                OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IZ")]);}

            iz_json = entryVal;

            // Need to do "- 1" because C++ starts in zero
            iz_json--;

            // If entry in json is zero, we get here -1
            // which is wrong, so need to send warning messsage
            // and skip entry
            if (iz_json == -1){
                // Through a warning invalid entry           
                msg_string = 
                    "<OpenWQ> WARNING: SS '" 
                    + elemName 
                    + "' cannot be zero. It needs to start in one (entry skipped): File=" 
                    + std::to_string(ssf+1)
                    + ", Sub_structure=" + std::to_string(ssi+1)
                    + ", Data_row=" + std::to_string(di + 1); 

                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                continue; // skip entry
            }

        }catch(...){

            // try as string for the cases where entry is "all"
            std::string entryVal = ""; // dummy variable
            // if JSON
            if (DataFormat.compare("JSON")==0){
                entryVal = EWF_SS_json_sub_rowi.at(8);}
            // if ASCII
            else if (DataFormat.compare("ASCII")==0){
                entryVal = ASCIIRowElemEntry[ 
                    OpenWQ_utils.FindStrIndexInVectStr(headerKeys,"IZ")];}

            validEntryFlag = getArrayElem(
                OpenWQ_wqconfig,
                OpenWQ_output,
                elemName,
                (std::string) entryVal,
                iz_json,
                ssf,    // SS file
                ssi,    // SS structure
                di);    // SS row

            if (!validEntryFlag){continue;}
        }

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

    */

}