// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
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

// For VTK support:
// This was built for hexahedrons (see line 66) for now (cubes, Rectangular cuboid, Trigonal trapezohedron, etc)
// But line 33:33 makes determining cubes of side lenght = 1 (but this can all be changed)
// based on https://lorensen.github.io/VTKExamples/site/Cxx/IO/WriteVTU/


#include "OpenWQ_utils.hpp"


// ########################################                                                          
// Get number of lines in file
// ########################################                                                          
unsigned int OpenWQ_utils::getNumLinesfromASCIIFile(
    std::string ascii_FilePath){ // path to file

    std::ifstream aFile (ascii_FilePath);   
    std::size_t lines_count=0;
    std::string line;

    // loop to get number of lines
    while (std::getline(aFile , line))
            ++lines_count;

    // Close file
    aFile.close();

    // return number of lines
    return lines_count;

}

// ########################################
// Split string into vector of strings based on string-delimiter
// ########################################
std::vector<std::string> OpenWQ_utils::StringSplit (
    std::string stringInput, 
    std::string delimiter) {

    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> stringVect;

    // Check if it contains \r,
    // if yes, remove it
    if (stringInput.find("\r") != std::string::npos){
        stringInput.erase(stringInput.find("\r"),2);}

    while ((pos_end = stringInput.find (delimiter, pos_start)) != std::string::npos) {
        // Get token
        token = stringInput.substr (pos_start, pos_end - pos_start);
        // Removing leading and trailling whitespaces
        token = RemoveStrLeadTrailWhiteSpaces(token);
        // update pos_start location
        pos_start = pos_end + delim_len;
        // add string to string-vector
        stringVect.push_back (token);
    }

    // Last sub-string
    token = stringInput.substr(pos_start);
    // Remove whitespaces
    token = RemoveStrLeadTrailWhiteSpaces(token);
    // Add last string
    stringVect.push_back (token);
    return stringVect;
}

// ########################################
// Converts string text to lower case
// ########################################

std::string OpenWQ_utils::ConvertStringToUpperCase(
    const std::string StrEntry){

    // Local variable
    std::string NewStr = StrEntry;

    // Try to transform to upper case if text
    try{
        std::transform(
        NewStr.begin(),
        NewStr.end(),
        NewStr.begin(), // Convert to lower case to avoid issues
        [](unsigned char c) { return std::toupper(c); });
        
    }catch (...){}

    return NewStr;

}

// Get valid time stamps
bool OpenWQ_utils::GetTimeStampsFromLogFile(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    std::string logFile_folderPath,             // Path where LogFile is location
    std::string preOutputRow_substring,         // Substring of the output to search
    std::vector<std::string>& timeStamps_vec,   // Vector with timestamps
    std::string errMsg_LofFileIdentifier){      // Logfile errMsg identifier

    // Local variables
    std::string rowEntryLogfile;
    std::string ewf_sim_logFile_path;
    std::string output_str;
    std::string msg_string;
    int pos;
    bool foundTimeStmps = true;

    // Full path to LogFile of EWF source
    ewf_sim_logFile_path = logFile_folderPath;
    ewf_sim_logFile_path.append("/");
    ewf_sim_logFile_path.append(OpenWQ_wqconfig.get_LogFile_name());

    // Open LogFile_name_fullpath file to read timestamps (read mode only)
    std::ifstream logFile_ewf (ewf_sim_logFile_path, std::ios::in);

    // Read and skip irrelevant lines
    if (logFile_ewf.is_open()){

        while(std::getline(logFile_ewf, rowEntryLogfile)){
            // Check if row contains output info
            if (rowEntryLogfile.find(preOutputRow_substring) != std::string::npos) {
                pos = preOutputRow_substring.size();
                timeStamps_vec.push_back(rowEntryLogfile.substr(pos));
            }
        }
        logFile_ewf.close();

        // If output info in supporting "logFile_folderPath" file not found, 
        // through warning message (entry skipped)
        if (timeStamps_vec.empty()){

            msg_string = 
                "<OpenWQ> WARNING: '"
                + errMsg_LofFileIdentifier 
                + " =" + ewf_sim_logFile_path
                + " found, but it contains no output/print info."
                " Make sure this is the correct log file (entry skipped).";
            // Print it (Console and/or Log file)
            OpenWQ_output.ConsoleLog(
                OpenWQ_wqconfig,    // for Log file name
                msg_string,         // message
                true,               // print in console
                true);              // print in log file

            foundTimeStmps = false;

        }

    }else{

        // If supporting "logFile_folderPath" not found, 
        // through warning message (entry skipped)
        msg_string = 
            "<OpenWQ> WARNING: SS/EWF '" 
                + errMsg_LofFileIdentifier 
                + " =" + ewf_sim_logFile_path
                + " not found. Make sure to include it in the same path "
                "as the SS/EWF h5 files (entry skipped).";
        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file
        
        foundTimeStmps = false;

    }

    return foundTimeStmps;
    
}


// Read JSON array element and return -1 if "all"
bool OpenWQ_utils::Convert2NegativeOneIfAll_inputInt(
    OpenWQ_wqconfig & OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    std::string errMsg_LofFileIdentifier,
    json json_array,
    int index_of_json_array,
    int & returnValue,
    int nDomain){

    // Valid entry
    bool validEntryFlag = false;
    std::string msg_string;

    try{
        // try entry as double
        returnValue = json_array.at(index_of_json_array);
        validEntryFlag = true;
    }catch(...){

        try{
            // try as string for the cases where entry is "all"
            std::string entryVal = json_array.at(index_of_json_array);
            // Check if "all" entry
            if (entryVal.compare("ALL")==0){
                returnValue = -1.0f;
                validEntryFlag = true;
            }else if (entryVal.compare("END")==0){
                returnValue = nDomain;
                validEntryFlag = true;       
            }else {validEntryFlag = false;}
            
        }catch(...){
            validEntryFlag = false;
        }
    }

    // Warning message if needed
    if(!validEntryFlag){

        msg_string = 
            "<OpenWQ> WARNING:'" 
                + errMsg_LofFileIdentifier 
                + " (entry skipped).";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file
    }

    return validEntryFlag;

}

// Read JSON keyVal with type: Int
int OpenWQ_utils::RequestJsonKeyVal_int(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    json json_struct,
    std::string jsonKey,
    std::string msgIndetifier,
    bool abort_flag){

    // Local variables
    int jsonVal_int;
    std::string varType = "integer";

    try{

        // Try to access json key
        jsonVal_int = json_struct[jsonKey]; 

    }catch(...){

        // Abort and through error message
        RequestJsonKeyVal_errorAbort(
            OpenWQ_wqconfig, 
            OpenWQ_output,
            jsonKey,
            msgIndetifier,
            varType,
            abort_flag);
        
    } 

    // If jsonVal found, return it
    return jsonVal_int;

}

// Read JSON keyVal with type: Int
double OpenWQ_utils::RequestJsonKeyVal_double(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    json json_struct,
    std::string jsonKey,
    std::string msgIndetifier,
    bool abort_flag){

    // Local variables
    double jsonVal_double;
    std::string varType = "integer";

    try{

        // Try to access json key
        jsonVal_double = json_struct[jsonKey]; 

    }catch(...){

        // Abort and through error message
        RequestJsonKeyVal_errorAbort(
            OpenWQ_wqconfig, 
            OpenWQ_output,
            jsonKey,
            msgIndetifier,
            varType,
            abort_flag);
        
    } 

    // If jsonVal found, return it
    return jsonVal_double;

}

// Read JSON keyVal with type: json
json OpenWQ_utils::RequestJsonKeyVal_json(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_output& OpenWQ_output,
    json json_struct,
    std::string jsonKey,
    std::string msgIndetifier,
    bool abort_flag){

    // Local variables
    json jsonVal_json;
    std::string varType = "json/sub-json structure";

    // Try to access json key
    jsonVal_json = json_struct[jsonKey]; 
        
    // If json empty,
    // abort and through error message
    if(jsonVal_json.empty()){
        
        RequestJsonKeyVal_errorAbort(
            OpenWQ_wqconfig, 
            OpenWQ_output,
            jsonKey,
            msgIndetifier,
            varType,
            abort_flag);

    }

    // If jsonVal found, return it
    return jsonVal_json;

}


// ########################################
// Removing leading and trailling whitespaces
// ########################################
std::string OpenWQ_utils::RemoveStrLeadTrailWhiteSpaces(
    std::string String2RemWhiteSpace){

    // Remove Leading Spaces
    String2RemWhiteSpace.erase(
        String2RemWhiteSpace.begin(), 
        std::find_if(String2RemWhiteSpace.begin(), 
        String2RemWhiteSpace.end(), 
        std::bind1st(std::not_equal_to<char>(), 
        ' ')));

    // Remove trailling spaces
    String2RemWhiteSpace.erase(
        std::find_if(
            String2RemWhiteSpace.rbegin(), 
            String2RemWhiteSpace.rend(), 
            std::bind1st(std::not_equal_to<char>(), 
            ' ')).base(), String2RemWhiteSpace.end());

    return String2RemWhiteSpace;

}