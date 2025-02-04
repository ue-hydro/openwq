

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

#include "models_TD/headerfile_td.hpp"


void OpenWQ_TD_model::TD_driver_run(
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_output& OpenWQ_output,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r, 
    const double wmass_source){
    
    // Local variables
    std::string msg_string;

    // if TD_module = OPENWQ_NATIVE_TD_ADVDISP
    if ((OpenWQ_wqconfig.TD_model->TD_module).compare("OPENWQ_NATIVE_TD_ADVDISP") == 0)
    {
        AdvDisp(
            OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);

    // if TD_module = NATIVE_TD_ADV
    }else if ((OpenWQ_wqconfig.TD_model->TD_module).compare("NATIVE_TD_ADV") == 0)
    {
        Adv(
            OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);
    }

}