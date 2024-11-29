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
    OpenWQ_utils& OpenWQ_utils,
    OpenWQ_output& OpenWQ_output){
            
    
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
    unsigned int ascii_numHeaderRows;           // additional information for ASCII data input
    unsigned int ascii_headerKeyRow;            // additional information for ASCII data input
    unsigned int lineCount;                     // row count for ASCII file
    std::string rowEntryASCII;                  // row data
    std::string headerRowASCII;                 // row with header data
    std::vector<std::string> ASCIIRowdata;      // all row data
    std::vector<double> row_data_entry;         // vector with header keys
    std::vector<std::vector<double>> table_data_entry;         // vector with header keys
    std::vector<std::string> row_data_entry_string;
    std::vector<std::string> headerKeys;        // vector with header keys
    std::string headerWord;                     // interactive header words
    std::string elemName;                       // temporary element name
    json EWF_SS_json_sub_rowi;                  // json row of json-substructure json_subStruct 
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
    std::string model_parameter;


    // Local variables
    std::string DataFormat;
    json json_subStruct;
    json json_paramData;
    json json_paramData_row;
    unsigned int index_elm;
    double num_to_pushBack;
    std::string json_element_str;

    // Data format
    DataFormat = std::get<3>(OpenWQ_wqconfig.TS_model->HypeMMF->info_vector);

    // Get PARAMETERS JSON data
    errorMsgIdentifier = "TS_model file";
    json_subStruct = OpenWQ_utils.RequestJsonKeyVal_json(
        OpenWQ_wqconfig, OpenWQ_output,
        OpenWQ_json.TS_module, "PARAMETERS",
        errorMsgIdentifier,
        true);


    model_parameter = "COHESION";

    // Get Parameter Data as JSON structure
    json_paramData = OpenWQ_utils.RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            json_subStruct, model_parameter,
            errorMsgIdentifier,
            true);
    
    /* ########
    // Check if DATA_FORMAT=JSON
    // if yes, then get file structure info
    ###########*/

    // If JSON
    if (DataFormat.compare("JSON")==0){

        // Get number of rows of data in JSON
        errorMsgIdentifier = errorMsgIdentifier + " > json block COHESION with DataFormat=" + DataFormat;
        
        // Get header keys
        json_paramData_row = OpenWQ_utils.RequestJsonKeyVal_json(
                OpenWQ_wqconfig, OpenWQ_output,
                json_paramData, std::to_string(0),
                errorMsgIdentifier,
                true);

        // Cast json into vector of strings
            for (unsigned int el=0;el<json_paramData_row.size();el++){
                headerKeys.push_back(json_paramData_row.at(el));
        }

        // check if field exists, return if not with error message
        num_rowdata = json_paramData.size()-1;

    }
    // Get additional info for ASCII
    else if (DataFormat.compare("ASCII")==0){

        try{

            // file path
            ascii_FilePath = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                json_paramData, "FILEPATH",
                errorMsgIdentifier,
                true);

            // delimiter
            ascii_delimiter = OpenWQ_utils.RequestJsonKeyVal_str(
                OpenWQ_wqconfig, OpenWQ_output,
                json_paramData, "DELIMITER",
                errorMsgIdentifier,
                true);

            // number of header rows
            ascii_numHeaderRows = OpenWQ_utils.RequestJsonKeyVal_int(
                OpenWQ_wqconfig, OpenWQ_output,
                json_paramData, "NUMBER_OF_HEADER_ROWS",
                errorMsgIdentifier,
                true);            

            // position of key header
            ascii_headerKeyRow = OpenWQ_utils.RequestJsonKeyVal_int(
                OpenWQ_wqconfig, OpenWQ_output,
                json_paramData, "HEADER_KEY_ROW",
                errorMsgIdentifier,
            true);  

            // Get number of rows in ASCII file
            num_rowdata = OpenWQ_utils.getNumLinesfromASCIIFile(
                OpenWQ_wqconfig,
                OpenWQ_output,
                ascii_FilePath)
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
                    + ascii_FilePath;
                // Print it (Console and/or Log file)
                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
                return;
            }
    
        }catch(...){
            // If there is an issue with the ASCII input data
            // through a warning message and skip entry
            msg_string = 
                "<OpenWQ> WARNING: SS/EWF '" 
                " load/sink/conc with ASCII format has an issue with json-keys or data structure (entry skipped): File=";
            // Print it (Console and/or Log file)
            OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            return;
        }
    }

    // ########################################
    // Loop over row data in tabular file
    // ######################################## //
    
    for (unsigned int di=1;di<num_rowdata;di++){

        // Reset the size to zero (the object will have no elements)
        row_data_col.reset();

        row_data_entry_string.clear();

        // Get row-json di from json_subStruct ["COHESION"]
        // Abort if not found
        errorMsgIdentifier = " json block 'DATA', row " 
                            + std::to_string(di) + " with DataFormat=" 
                            + DataFormat;

        // Get the data as vector(std::string)
        // both cases, JSON and ASCII
        if (DataFormat.compare("JSON")==0){

            json_paramData_row = OpenWQ_utils.RequestJsonKeyVal_json(
                OpenWQ_wqconfig, OpenWQ_output,
                json_paramData, std::to_string(di),
                errorMsgIdentifier,
                true);
            
            // Cast json into vector of strings
            for (unsigned int el=0;el<json_paramData_row.size();el++){

                // Get element (convert to string)
                json_element_str = to_string(json_paramData_row.at(el));

                row_data_entry_string.push_back(json_element_str);
            }

        // If DataFormat=ASCII, then get row data 
        // and convert-to-upper-case and split it by element entry
        }else if (DataFormat.compare("ASCII")==0){
            row_data_entry_string = 
                OpenWQ_utils.StringSplit(
                    OpenWQ_utils.ConvertStringToUpperCase(ASCIIRowdata[di]),
                    ascii_delimiter);

            // Cast json into vector of strings
            
        }

        // Now, go element by element and 
        // Save number or save to -1 if entry='ALL'
        
        row_data_entry.clear();

        for (unsigned int el=0;el<row_data_entry_string.size();el++){

            try{
                num_to_pushBack = std::stod(row_data_entry_string.at(el));
                
            }catch(...){

                // Case where the input is not a number

                if (row_data_entry_string.at(el).compare(R"(\"ALL\")")){

                    num_to_pushBack = -1;

                }else{

                    msg_string = 
                    "<OpenWQ> WARNING: Invalid entry. Expected number or 'ALL' but got : " 
                    + to_string(json_paramData_row.at(el)) +
                    "(entry ignored)";

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                    // skip entry
                    continue;

                }
            }

            // add entry
            row_data_entry.push_back(num_to_pushBack); 

        }

        table_data_entry.push_back(row_data_entry);

    }
                
        
        // ###################
        // ix
        // ###################
        elemName = "IX";
        index_elm = OpenWQ_utils.FindStrIndexInVectStr(headerKeys, elemName);

        /*
    
        // Try if number
        try{

            // try entry as int
            int entryVal = stoi(row_data_entry.at(index_elm));

            ix_json = entryVal;

            // Need to do "- 1" because C++ starts in zero
            ix_json--;

            // If entry in json is zero, we get here -1
            // which is wrong, so need to send warning messsage
            // and skip entry
            if (ix_json >= -1){
                // Through a warning invalid entry           
                msg_string = 
                    "<OpenWQ> WARNING: SS '" 
                    + elemName 
                    + "' cannot be zero. It needs to start in one (entry skipped): File=";   

                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                continue; // skip entry
            }

        }catch(...){

            std::string entryVal = row_data_entry.at(index_elm);

            if(entryVal.compare("ALL") == 0){
        
                entryVal = OpenWQ_wqconfig.get_allSS_flag();

            }

            

            //validEntryFlag = getArrayElem_SS(
            //    OpenWQ_wqconfig,
            //    OpenWQ_output,
            //    elemName,
            //    (std::string) entryVal,
            //    ix_json,
            //    ssf,    // SS file
            //    ssi,    // SS structure
            //    di);    // SS row

            if (!validEntryFlag){continue;}

        }

        */
        
        /*
        // Get the vector with the data
        row_data_col = {
            (double) ix_json,
            (double) iy_json,
            (double) iz_json,
            value,
            }; 

        // Add new row to SinkSource_FORC or ExtFlux_FORC_jsonAscii
        AppendRow_SS_EWF_FORC_jsonAscii(
            OpenWQ_wqconfig,
            inputType,
            row_data_col);

        */

}
