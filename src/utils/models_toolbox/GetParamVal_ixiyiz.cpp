
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

// Get relevant parameter vals at ix, iy, iz
double OpenWQ_utils::GetParamVal_ixiyiz(
    TupleParameterData Param_entryTuple,
    unsigned ix, 
    unsigned iy, 
    unsigned iz
  ){

  // Local variables
  unsigned int num_rows;
  std::vector<std::string> HeaderKeys;

  HeaderKeys = std::get<0>(Param_entryTuple);
  num_rows = std::get<1>(Param_entryTuple).size(); 

    

}