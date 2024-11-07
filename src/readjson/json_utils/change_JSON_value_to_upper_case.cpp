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

#include "readjson/headerfile_readjson.hpp"

// ########################################
// Change JSON value
// ########################################
void OpenWQ_readjson::change_JSON_value_to_upper_case(
    OpenWQ_utils &OpenWQ_utils,
    json &object,
    std::string new_jsonkey_layer_1){
    
    // Local variables
    long num_vals;
    bool array_flag = false;


    // Vector of words to keep with native upper/lower case provided by user
    std::vector<std::string> exclude_vales;
    exclude_vales.push_back("FILEPATH");
    exclude_vales.push_back("FOLDERPATH");
    exclude_vales.push_back("MODULE_CONFIG_FILEPATH");
    unsigned int num_exclude_vales = exclude_vales.size();

    // if FILEPATH or FOLDERPATH, do not change to upper case
    for (unsigned int w=0;w<num_exclude_vales;w++){

        if (new_jsonkey_layer_1.find(exclude_vales[w]) != std::string::npos){
                return;
        }

    }

    // check if object is an array
    if(object[new_jsonkey_layer_1].type() == json::value_t::array)
        array_flag = true;        
    
    // ########################################
    // Replace elements by upper case version
    // ########################################

    // if string (not array or tuple) - most of the cases
    if (array_flag != true){ // check if value is array

        // Local variables
        
        std::string new_value;
        std::string old_value;

        num_vals = 1;
        old_value = object[new_jsonkey_layer_1];

        // Convert to upper case
        new_value = OpenWQ_utils.ConvertStringToUpperCase(old_value);    

        // replace value
        object[new_jsonkey_layer_1].swap(new_value);         

    // If value is an array or tuple 
    // (loop over the array elements to replace lower case string to upper case strings
    }else{

        // initiate json arrays and get existing one 
        json::array_t old_value_array = object[new_jsonkey_layer_1];
        json::array_t new_value_array;
        
        // size of array
        num_vals = old_value_array.size();

        // Loop over number of values
        // if NOT FILEPATH or FOLDERPATH, change string to Upper Case
        for (unsigned int i=0;i<num_vals;i++){
                
            // Get element in array/tuple
            auto old_value = old_value_array.at(i);
            
            // if string
            try{         
                std::string new_value_str;

                // Convert to upper case
                // Convert to upper case
                new_value_str = OpenWQ_utils.ConvertStringToUpperCase(old_value);

                // add string to new_value_array
                new_value_array.push_back(new_value_str);

            // if num
            }catch(...){
                double new_value_num = old_value;
                
                // add num to new_value_array
                new_value_array.push_back(new_value_num);
            }

        }                        

        // replace value
        object[new_jsonkey_layer_1].swap(new_value_array); 

    }

}