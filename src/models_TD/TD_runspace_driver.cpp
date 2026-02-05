

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

#include "models_TD/headerfile_TD.hpp"


void OpenWQ_TD_model::TD_driver_run(
    OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
    OpenWQ_wqconfig& OpenWQ_wqconfig,
    OpenWQ_vars& OpenWQ_vars,
    OpenWQ_output& OpenWQ_output,
    const int source, const int ix_s, const int iy_s, const int iz_s,
    const int recipient, const int ix_r, const int iy_r, const int iz_r,
    const double wflux_s2r,
    const double wmass_source){

    // OPTIMIZED: use cached bools instead of string comparisons
    if (OpenWQ_wqconfig.is_TD_advdisp)
    {
        AdvDisp(
            OpenWQ_hostModelconfig,
            OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);

    }else if (OpenWQ_wqconfig.is_TD_adv)
    {
        Adv(
            OpenWQ_vars, OpenWQ_wqconfig,
            source, ix_s, iy_s, iz_s,
            recipient, ix_r, iy_r, iz_r,
            wflux_s2r, wmass_source);
    }

}