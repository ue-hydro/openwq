

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


#ifndef OPENWQ_HEADERFILE_STH_INCLUDED
#define OPENWQ_HEADERFILE_STH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>
#include "global/exprtk.hpp"
#include <cstdio>

#include "global/openwq_json.hpp"
#include "global/openwq_vars.hpp"
#include "readjson/headerfile_RJSON.hpp"
#include "extwatflux_ss/headerfile_EWF_SS.hpp"
#include "global/openwq_wqconfig.hpp"
#include "global/openwq_hostmodelconfig.hpp"
#include "utils/headerfile_UTILS.hpp"
#include "units/headerfile_UNITS.hpp"
#include "initiate/headerfile_INIT.hpp"
#include "output/headerfile_OUT.hpp"
#include "compute/headerfile_compute.hpp"


class OpenWQ_ST_model{

    public:

        // ##########################
        // MAIN DRIVER 
        // => Run ST model
        // ##########################
        
        void ST_driver_run(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_json& OpenWQ_json,                       // create OpenWQ_json object
            OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
            OpenWQ_units& OpenWQ_units,                     // functions for unit conversion
            OpenWQ_utils& OpenWQ_utils,                       // utility methods/functions
            OpenWQ_readjson& OpenWQ_readjson,               // read json files
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,     // sink and source modules)
            OpenWQ_compute& OpenWQ_compute,
            OpenWQ_output& OpenWQ_output,
            time_t simtime,                                 // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r, const double wmass_source);


        // ##########################
        // MODEL OPTIONS 
        // ##########################

        // HBV sediments as implemented in HYPE
        void hype_hbvsed(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            double wflux_s2r, 
            double wmass_source);

};

#endif