

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


#ifndef OPENWQ_HEADERFILE_TSH_INCLUDED
#define OPENWQ_HEADERFILE_TSH_INCLUDED

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


class OpenWQ_TS_model{

    public:

        // ##########################
        // MAIN DRIVERS
        // ##########################

         // ###########################     
        // => Run Config method

        void TS_driver_config(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // ###########################     
        // => Run ST model
       
        void TS_driver_run(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,               // create OpenWQ_wqconfig object
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output,                                   // simulation time in seconds since seconds since 00:00 hours, Jan 1, 1970 UTC
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r, const double wmass_source,
            std::string TS_type);

    private:
    
        // ##########################
        // MODEL OPTIONS 
        // ##########################

        // ##########################
        // CONFIG

        // HYPE_MMF
        void TS_HYPE_MMF_config(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output);

        // HYPE_HBVSED
        void TS_HYPE_HBVSED_config(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_output& OpenWQ_output);

        // ##########################
        // RUN
        
        // HBV sediments as implemented in HYPE
        void hbvsed_hype_erosion_run(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, // / no need for ix_r, iy_r, iz_r (because mobilization occurs at the adjacent cell of the upper compartment)
            double wflux_s2r, 
            double wmass_source,
            std::string TS_type);
            
         // Based on the Morgan-Morgan-Finney erosion model implemented in HYPE
        void mmf_hype_erosion_run(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, // / no need for ix_r, iy_r, iz_r (because mobilization occurs at the adjacent cell of the upper compartment)
            double wflux_s2r, 
            double wmass_source,
            std::string TS_type);

        

};

#endif