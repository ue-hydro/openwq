
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

#include "utils/headerfile_UTILS.hpp"

// Get Model Parameter data as tuple<vector>
// Parameter data can be provided as JSON or ASCII
       
arma::Cube<double>  OpenWQ_utils::LoadCubeComprtModelParameters_asARMACUBE_JSONorASCII(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    std::string DataFormat,
    std::string model_parameter,
    unsigned int icmp,
    json json_PARAMETER_subStruct,
    json json_PARAMETER_DEFAULTS_subStruct,
    std::string errorMsgIdentifier){

    // Local variables
    std::string errorMsgIdentifier_local;               // error message composition
    double default_parameter_val;                       // default parameter value taken from JSON
    json json_paramData;                                // entire json entry for parameter
    json json_paramData_row;                            // sub-json struture with parameter entry row
    unsigned int num_rowdata;                           // number of json or ASCII-file rows
    std::vector<std::string> headerKeys;                // vector with header keys
    std::string ascii_FilePath;                         // additional information for ASCII data input
    double num_to_pushBack;                             // numerical value to add to parameter vector
    std::string json_element_str;                       // json entry casted as string
    std::string ascii_delimiter;                        // additional information for ASCII data input
    unsigned int ascii_numHeaderRows;                   // additional information for ASCII data input
    unsigned int ascii_headerKeyRow;                    // additional information for ASCII data input
    std::string headerRowASCII;                         // row with header data
    unsigned int lineCount;                             // row count for ASCII file
    std::string rowEntryASCII;                          // row data
    std::vector<std::string> ASCIIRowdata;              // all row data
    std::string msg_string;                             // error/warning message string
    arma::vec row_data_col;                             // new row data (initially as col data)
    std::vector<std::string> row_data_entry_string;     // full row entry casted as vector of strings
    std::vector<double> row_data_entry;                 // vector with header keys
    std::vector<std::vector<double>> table_data_entry;  // parameter values as vector<vector>         
    typedef std::tuple<
        std::vector<std::string>,
        arma::Cube<double>
        > TupleParameterData;                           // define data type: vector with header keys
    TupleParameterData tuple_data_entry;
    unsigned int ix_loc, iy_loc, iz_loc, val_loc;       // indexes of locations in JSON or ASCII file
    int ix_num_str, iy_num_str, iz_num_str;    // iterative ix, iy, and iz indexes
    int ix_num_end, iy_num_end, iz_num_end;    // iterative ix, iy, and iz indexes
    double value_num;                                   // iterative parameter value number for arma::cube
    unsigned int icmp_x_num, icmp_y_num, icmp_z_num;

    // get num of cells
    icmp_x_num = OpenWQ_hostModelconfig.get_HydroComp_num_cells_x_at(icmp);
    icmp_y_num = OpenWQ_hostModelconfig.get_HydroComp_num_cells_y_at(icmp);
    icmp_z_num = OpenWQ_hostModelconfig.get_HydroComp_num_cells_z_at(icmp);

    // generate armaCube_data_entry for relevant icmp
    arma::Cube<double> armaCube_data_entry(icmp_x_num, icmp_y_num, icmp_z_num);

    /* ########################################
    //  Get default parameter value 
    ######################################## */
    errorMsgIdentifier_local = errorMsgIdentifier + " > Default parmeter value";
    
    // Compulsary to be there
    default_parameter_val = RequestJsonKeyVal_double(
        OpenWQ_wqconfig, OpenWQ_output,
        json_PARAMETER_DEFAULTS_subStruct, model_parameter,
        errorMsgIdentifier_local,
        true);

    // set all ix,iy,iz parameter values in arma::cube 
    // to default parameter values
    armaCube_data_entry.fill(default_parameter_val);

    /* ########################################
    // Get entire JSON entry to parameter 
    // Not compulsary to be there
    ######################################## */
    errorMsgIdentifier_local = errorMsgIdentifier + " > Parameter entires (JSON or ASCII)";

    // Get Parameter Data as JSON structure
    errorMsgIdentifier_local = errorMsgIdentifier_local;
    json_paramData = RequestJsonKeyVal_json(
            OpenWQ_wqconfig, OpenWQ_output,
            json_PARAMETER_subStruct, model_parameter,
            errorMsgIdentifier_local,
            false);

    // If JSON info found for discrete parameter dara,
    // then overwrite the default values where applicable 
    if (!json_paramData.empty()){

        /* ########################################
        //  Get Header keys (JSON or ASCII)
        ######################################## */

        // Get number of rows of data in JSON
        errorMsgIdentifier_local = errorMsgIdentifier + " > DataFormat=" + DataFormat;

        // If JSON
        if (DataFormat.compare("JSON")==0){

            // Get header keys
            json_paramData_row = RequestJsonKeyVal_json(
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
                ascii_FilePath = RequestJsonKeyVal_str(
                    OpenWQ_wqconfig, OpenWQ_output,
                    json_paramData, "FILEPATH",
                    errorMsgIdentifier,
                    true);

                // delimiter
                ascii_delimiter = RequestJsonKeyVal_str(
                    OpenWQ_wqconfig, OpenWQ_output,
                    json_paramData, "DELIMITER",
                    errorMsgIdentifier,
                    true);

                // number of header rows
                ascii_numHeaderRows = RequestJsonKeyVal_int(
                    OpenWQ_wqconfig, OpenWQ_output,
                    json_paramData, "NUMBER_OF_HEADER_ROWS",
                    errorMsgIdentifier,
                    true);            

                // position of key header
                ascii_headerKeyRow = RequestJsonKeyVal_int(
                    OpenWQ_wqconfig, OpenWQ_output,
                    json_paramData, "HEADER_KEY_ROW",
                    errorMsgIdentifier,
                true);  

                // Get number of rows in ASCII file
                num_rowdata = getNumLinesfromASCIIFile(
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
                            headerRowASCII = ConvertStringToUpperCase(headerRowASCII);

                            // Parse the headerRowASCII to get a vector with the json entries
                            // parsing based on the delimiter ascii_delimiter defined by the user
                            headerKeys = StringSplit(
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
                        "<OpenWQ> WARNING: " + errorMsgIdentifier_local +
                        " ASCII file cannot be found (entry skipped): File=" 
                        + ascii_FilePath;

                    // Print it (Console and/or Log file)
                    OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);

                }
        
            }catch(...){

                // If there is an issue with the ASCII input data
                // through a warning message and skip entry
                msg_string = 
                    "<OpenWQ> WARNING:'" + errorMsgIdentifier_local
                    + " ASCII format has an issue with json-keys or data structure (entry skipped): File="
                    + ascii_FilePath;
                // Print it (Console and/or Log file)

                OpenWQ_output.ConsoleLog(OpenWQ_wqconfig, msg_string, true, true);
            }
        }

        /* ########################################
        //  Loop over row data in tabular file (JSON or ASCII)
        ######################################## */
        
        errorMsgIdentifier_local = errorMsgIdentifier + DataFormat; 

        for (unsigned int di=0;di<num_rowdata;di++){

            // Reset the size to zero (the object will have no elements)
            row_data_col.reset();
            row_data_entry_string.clear();

            // Get the data as vector(std::string)
            // both cases, JSON and ASCII
            if (DataFormat.compare("JSON")==0){

                json_paramData_row = RequestJsonKeyVal_json(
                    OpenWQ_wqconfig, OpenWQ_output,
                    json_paramData, std::to_string(di+1),
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
                    StringSplit(
                        ConvertStringToUpperCase(ASCIIRowdata[di]),
                        ascii_delimiter);
                
            }

            // Now, go element by element and 
            // save number or save to -1 if entry='ALL'
            
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

        /* ########################################
        // Convert 
        // IN) table_data_entry (as vector<vector<double>>)
        // OUT) armaCube_data_entry
        ######################################## */

        // get ix, iy, iz and vals indexes
        auto ix_loc_itrc = std::find( begin(headerKeys), end(headerKeys), "IX" );
        auto iy_loc_itrc = std::find( begin(headerKeys), end(headerKeys), "IY" );
        auto iz_loc_itrc = std::find( begin(headerKeys), end(headerKeys), "IZ" );
        auto val_loc_itrc = std::find( begin(headerKeys), end(headerKeys), "VALUE" );

        ix_loc = std::distance(headerKeys.begin(), ix_loc_itrc);
        iy_loc = std::distance(headerKeys.begin(), iy_loc_itrc);
        iz_loc = std::distance(headerKeys.begin(), iz_loc_itrc);
        val_loc = std::distance(headerKeys.begin(), val_loc_itrc);

        // loop over entires

        for (int row_i = 0; row_i < table_data_entry.size(); row_i++){

            // get values in row
            // need to do "-1" because c++ starts at zero, but
            // json entries start at 1
            ix_num_str = table_data_entry.at(row_i).at(ix_loc) - 1;
            ix_num_end = ix_num_str;
            iy_num_str = table_data_entry.at(row_i).at(iy_loc) - 1;
            iy_num_end = iy_num_str;
            iz_num_str = table_data_entry.at(row_i).at(iz_loc) - 1;
            iz_num_end = iz_num_str;
            value_num = table_data_entry.at(row_i).at(val_loc);

            // check if values = -1
            // this is the case of "ALL"
            // below it sayd "-2" because the variables are subtracted "-1" above
            if (ix_num_str == -2){
                ix_num_str = 0;
                ix_num_end = icmp_x_num - 1;
            }
            if (iy_num_str == -2){
                iy_num_str = 0;
                iy_num_end =  icmp_y_num - 1;
            }
            if (iz_num_str == -2){
                iz_num_str = 0;
                iz_num_end =  icmp_z_num - 1; 
            }

            armaCube_data_entry(
                arma::span(ix_num_str, ix_num_end), 
                arma::span(iy_num_str, iy_num_end), 
                arma::span(iz_num_str, iz_num_end))
                .fill(value_num);

        }

    }

    return armaCube_data_entry;

}