
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

#include "models_LE/headerfile_LE.hpp"


void OpenWQ_LE_model::LE_driver_run(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_output& OpenWQ_output,                             // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r, const double wmass_source){
    
    // Local variables
    std::string msg_string;
    
    // #################################################
    // LE module
    // Lateral Exchange model (mobile material only)
    // #################################################

    if (OpenWQ_wqconfig.LE_module.compare("NATIVE_LE_BOUNDMIX") == 0)
    {
        // Boundary Mixing due to velocity gradients
        // due to turbulence and cross-boarder eddies
        // only apply if fluxe between cells in same compartment          
        native_BoundMix(
            OpenWQ_hostModelconfig, OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);

    // if TD_module != NONE (and any of the others)
    }else if ((OpenWQ_wqconfig.LE_module).compare("NONE")!= 0)
    {
        // Create Message
        msg_string = 
            "<OpenWQ> ERROR: No LE_module found or unkown";

        // Print it (Console and/or Log file)
        OpenWQ_output.ConsoleLog(
            OpenWQ_wqconfig,    // for Log file name
            msg_string,         // message
            true,               // print in console
            true);              // print in log file
        
        // Abort (Fatal error)
        exit(EXIT_FAILURE);

    }

    /* IntMOB ----> maybe recover in the future
    // Internal mobilization of immobile pools
    // Erosion and weathering
    OpenWQ_TD_model.IntMob(
        OpenWQ_vars, OpenWQ_wqconfig,
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source);
    */

}