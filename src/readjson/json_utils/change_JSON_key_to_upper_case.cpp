
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

#include "readjson/headerfile_RJSON.hpp"

// ########################################
// Change JSON key
// ########################################
void OpenWQ_readjson::change_JSON_key_to_upper_case(
    OpenWQ_utils &OpenWQ_utils,
    json &object, 
    const std::string& old_key,
    std::string& new_key)
{

    try{

        // Convert to upper case
        new_key = OpenWQ_utils.ConvertStringToUpperCase(old_key);

        // get iterator to old key; TODO: error handling if key is not present
        json::iterator it = object.find(old_key);

        // create null value for new key and swap value from old key
        std::swap(object[new_key], it.value());

    }catch(const std::exception &e){}

}