
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

#include "readjson/headerfile_RJSON.hpp"


/* ########################################
// Convert JSON text to Upper Case
######################################### */
void OpenWQ_readjson::ConvertJSONtext_2upperCase(
    OpenWQ_utils &OpenWQ_utils,
    json &jsondata){

    // Local Variables
    std::string old_jsonkey_layer_1;
    std::string old_jsonkey_layer_2;
    std::string old_jsonkey_layer_3;
    std::string old_jsonkey_layer_4;
    std::string old_jsonkey_layer_5;
    std::string old_jsonkey_layer_6;
    std::string old_jsonkey_layer_7;

    std::string new_jsonkey_layer_1;
    std::string new_jsonkey_layer_2;
    std::string new_jsonkey_layer_3;
    std::string new_jsonkey_layer_4;
    std::string new_jsonkey_layer_5;
    std::string new_jsonkey_layer_6;
    std::string new_jsonkey_layer_7;
    
    // ########################################
    // Identify all keys in JSON tree data and replace ny upper case keys

    // ------------> First layer of keys
    for (auto& x1 : jsondata.items()){

        try{
                
            old_jsonkey_layer_1 = x1.key();
            
            // Save key to o if a key
            if (!old_jsonkey_layer_1.empty()){ // if not null
                OpenWQ_readjson::change_JSON_key_to_upper_case(
                    OpenWQ_utils,
                    jsondata, 
                    old_jsonkey_layer_1,
                    new_jsonkey_layer_1);
            }else{
                continue;
            }
            
            // -------------> Second possible layer of keys
            for (auto& x2 : jsondata[new_jsonkey_layer_1].items()){
                    
                try{
                    
                    old_jsonkey_layer_2 = x2.key();

                    if (!old_jsonkey_layer_2.empty()&& // if not null
                        jsondata[new_jsonkey_layer_1].type() != json::value_t::array){ // and not array
                        
                        OpenWQ_readjson::change_JSON_key_to_upper_case(
                            OpenWQ_utils,
                            jsondata[new_jsonkey_layer_1], 
                            old_jsonkey_layer_2,
                            new_jsonkey_layer_2);
                    }else{
                        OpenWQ_readjson::change_JSON_value_to_upper_case(
                            OpenWQ_utils,
                            jsondata,
                            new_jsonkey_layer_1);
                        continue;
                    }

                    // -------------> Third possible layer of keys
                    for (auto& x3 : jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2].items()){
                                
                        try{
                        
                            old_jsonkey_layer_3 = x3.key();

                            if (!old_jsonkey_layer_3.empty()&& // if not null
                                jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2].type() != json::value_t::array){ // and not array
                                
                                OpenWQ_readjson::change_JSON_key_to_upper_case(
                                    OpenWQ_utils,
                                    jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2], 
                                    old_jsonkey_layer_3,
                                    new_jsonkey_layer_3);
                            }else{
                                OpenWQ_readjson::change_JSON_value_to_upper_case(
                                    OpenWQ_utils,
                                    jsondata[new_jsonkey_layer_1],
                                    new_jsonkey_layer_2);
                                continue;
                            }

                            // -------------> Fourth possible layer of keys
                            for (auto& x4 : jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3].items()){
                                    
                                try{
                                
                                    old_jsonkey_layer_4 = x4.key();

                                    if (!old_jsonkey_layer_4.empty() && // if not null
                                        jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3].type() != json::value_t::array){ // and not array
                                        
                                        OpenWQ_readjson::change_JSON_key_to_upper_case(
                                            OpenWQ_utils,
                                            jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3], 
                                            old_jsonkey_layer_4,
                                            new_jsonkey_layer_4);
                                    }else{
                                        OpenWQ_readjson::change_JSON_value_to_upper_case(
                                        OpenWQ_utils,
                                        jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2],
                                        new_jsonkey_layer_3);
                                    continue;
                                    }   

                                    // -------------> Fifth possible layer of keys
                                    for (auto& x5 : jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4].items()){
                            
                                        try{
                                        
                                            old_jsonkey_layer_5 = x5.key();

                                            if (!old_jsonkey_layer_5.empty() && // if not null
                                                jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4].type() != json::value_t::array){ // and not array
                                                
                                                OpenWQ_readjson::change_JSON_key_to_upper_case(
                                                    OpenWQ_utils,
                                                    jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4], 
                                                    old_jsonkey_layer_5,
                                                    new_jsonkey_layer_5);
                                            }else{
                                                OpenWQ_readjson::change_JSON_value_to_upper_case(
                                                OpenWQ_utils,
                                                jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3],
                                                new_jsonkey_layer_4);
                                                continue;
                                            }

                                            // -------------> Sixth possible layer of keys
                                            for (auto& x6 : jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5].items()){
                                
                                                try{
                                                
                                                    old_jsonkey_layer_6 = x6.key();

                                                    if (!old_jsonkey_layer_6.empty() && // if not null
                                                        jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5].type() != json::value_t::array){ // and not array
                                                        
                                                        OpenWQ_readjson::change_JSON_key_to_upper_case(
                                                            OpenWQ_utils,
                                                            jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5], 
                                                            old_jsonkey_layer_6,
                                                            new_jsonkey_layer_6);
                                                    }else{
                                                        OpenWQ_readjson::change_JSON_value_to_upper_case(
                                                        OpenWQ_utils,
                                                        jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4],
                                                        new_jsonkey_layer_5);
                                                        continue;
                                                    }

                                                    // -------------> Seven possible layer of keys
                                                    for (auto& x7 : jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5][new_jsonkey_layer_6].items()){
                                        
                                                        try{
                                                        
                                                            old_jsonkey_layer_7 = x7.key();

                                                            if (!old_jsonkey_layer_7.empty() && // if not null
                                                                jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5][new_jsonkey_layer_6].type() != json::value_t::array){ // and not array
                                                                
                                                                OpenWQ_readjson::change_JSON_key_to_upper_case(
                                                                    OpenWQ_utils,
                                                                    jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5][new_jsonkey_layer_6], 
                                                                    old_jsonkey_layer_7,
                                                                    new_jsonkey_layer_7);
                                                            }else{
                                                                OpenWQ_readjson::change_JSON_value_to_upper_case(
                                                                OpenWQ_utils,
                                                                jsondata[new_jsonkey_layer_1][new_jsonkey_layer_2][new_jsonkey_layer_3][new_jsonkey_layer_4][new_jsonkey_layer_5],
                                                                new_jsonkey_layer_6);
                                                                continue;
                                                            }

                                                        }catch (...){}

                                                    }

                                                }catch (...){}
                                                
                                            }

                                        }catch (...){}
                                    }                                   
                                        
                                }catch (...){}
                            }

                        }catch (...){}
                    }
                            
                }catch (...){}
            }

        }catch (...){}

    }
    

}