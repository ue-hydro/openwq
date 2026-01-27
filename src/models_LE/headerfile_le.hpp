

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


#ifndef OPENWQ_HEADERFILE_LEH_INCLUDED
#define OPENWQ_HEADERFILE_LEH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>
#include "global/exprtk.hpp"
#include <cstdio>

#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_vars.hpp"
#include "readjson/headerfile_RJSON.hpp"
#include "extwatflux_ss/headerfile_EWF_SS.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_hostModelConfig.hpp"
#include "utils/headerfile_UTILS.hpp"
#include "units/headerfile_units.hpp"
#include "initiate/headerfile_INIT.hpp"
#include "output/headerfile_OUT.hpp"
#include "compute/headerfile_compute.hpp"

#include "models_CH/headerfile_CH.hpp"
#include "models_TD/headerfile_TD.hpp"


class OpenWQ_LE_model{

    public:

        // ##########################
        // MAIN DRIVER 
        // => Run LE_model model
        // ##########################
        
        void LE_driver_run(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_output& OpenWQ_output,                               // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r, const double wmass_source);


        // ##########################
        // MODEL OPTIONS 
        // ##########################

        // Boundary Mixing due to velocity gradients
        // due to turbulence and cross-boarder eddies
        void native_BoundMix(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, // no need for ix_r, iy_r, iz_r (because mobilization occurs at the adjacent cell of the upper compartment)
            double wflux_s2r, 
            double wmass_source);

         // Internal mobilization of immobile pools
        // Erosion and weathering
        //void IntMob(
        //    OpenWQ_vars& OpenWQ_vars,
        //    OpenWQ_wqconfig& OpenWQ_wqconfig,
        //   const int source, const int ix_s, const int iy_s, const int iz_s,
        //   const int recipient, const int ix_r, const int iy_r, const int iz_r,
        //    double wflux_s2r, 
        //    double wmass_source);

};

#endif