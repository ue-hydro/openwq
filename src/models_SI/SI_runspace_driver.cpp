
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

#include "models_SI/headerfile_SI.hpp"


void OpenWQ_SI_model::SI_driver_run(
    OpenWQ_json& OpenWQ_json,
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_units& OpenWQ_units,
    OpenWQ_output& OpenWQ_output){
    
    // Local variables
    std::string msg_string;
    
    // #################################################
    // ST module
    // Sediment Transport
    // #################################################

    if (OpenWQ_wqconfig.SI->SI_module.compare("LANGMUIR") == 0)
    {
        
        langmuir(
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_wqconfig,
            OpenWQ_hostModelconfig,
            OpenWQ_output);

    // if TD_module != NONE (and any of the others)
    }else if (OpenWQ_wqconfig.SI->SI_module.compare("FREUNDLICH") == 0)
    {
        freundlich(
            OpenWQ_json,
            OpenWQ_vars,
            OpenWQ_wqconfig,
            OpenWQ_hostModelconfig,
            OpenWQ_output);
    }
    
}