
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


#ifndef OPENWQ_UTILS_INCLUDED
#define OPENWQ_UTILS_INCLUDED

#include <sys/stat.h>

#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "solver/OpenWQ_solver.hpp"
#include "output/OpenWQ_output.hpp"

class OpenWQ_utils{

    private:
         // General JSON key null error
        const std::string jsonKeyNull_msg_start_abort = "<OpenWQ> Execution ABORTED!\nExpected json value for key=";
        const std::string jsonKeyNull_msg_end_abort = " but not found! Revise the JSON files.";
        const std::string jsonKeyNull_msg_start_NOabort = "<OpenWQ> WARNING: Expected json value for key=";
        const std::string jsonKeyNull_msg_end_NOabort = " but not found! The entry has been zeroed. Make sure this was intended!";
        
        
        
        // Methods only needed by this class
        /*
         * JSON key null error
        */
        std::string get_jsonKeyNull_msg_start_abort();
        std::string get_jsonKeyNull_msg_end_abort();
        std::string get_jsonKeyNull_msg_start_NOabort();
        std::string get_jsonKeyNull_msg_end_NOabort();

        // Removing leading and trailling whitespaces
        std::string RemoveStrLeadTrailWhiteSpaces(
                std::string String2RemWhiteSpace);

        
        void RequestJsonKeyVal_errorAbort(
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_output& OpenWQ_output,
            std::string jsonKey,
            std::string msgIndetifier,
            std::string varType,
            bool abort_flag);

    public:

    // Get number of lines in file
    unsigned int getNumLinesfromASCIIFile(
        std::string ascii_FilePath); // path to file

    // Split string into vector of strings based on string-delimiter
    std::vector<std::string> StringSplit (
        std::string stringInput, 
        std::string delimiter);

    // Find index of string in vector::string
    int FindStrIndexInVectStr(
        std::vector<std::string> VectString,
        std::string Str2Find);

    // Convert string to upper case
    std::string ConvertStringToUpperCase(
        const std::string StrEntry);

    int getNumberOfDaysInMonthYear(
        const unsigned int YYYY_check,          // json: Year 
        const unsigned int MM_check);           // json: Month

    // Get timestamps from logFile
    bool GetTimeStampsFromLogFile(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_output& OpenWQ_output,
        std::string logFile_folderPath,
        std::string preOutputRow_substring,
        std::vector<std::string>& timeStamps_vec,
        std::string errMsg_LofFileIdentifier); 

    bool Convert2NegativeOneIfAll_inputInt(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_output& OpenWQ_output,
        std::string errMsg_LofFileIdentifier,
        json json_array,
        int index_of_json_array,
        int & returnValue,
        int nDomain);

    std::string RequestJsonKeyVal_str(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_output& OpenWQ_output,
        json json_struct,
        std::string jsonKey,
        std::string msgIndetifier,
        bool abort_flag);

    int RequestJsonKeyVal_int(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_output& OpenWQ_output,
        json json_struct,
        std::string jsonKey,
        std::string msgIndetifier,
        bool abort_flag);

    double RequestJsonKeyVal_double(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_output& OpenWQ_output,
        json json_struct,
        std::string jsonKey,
        std::string msgIndetifier,
        bool abort_flag);

    json RequestJsonKeyVal_json(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_output& OpenWQ_output,
        json json_struct,
        std::string jsonKey,
        std::string msgIndetifier,
        bool abort_flag);



};

#endif