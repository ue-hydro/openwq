 

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


#ifndef OPENWQ_TD_HEADERH_INCLUDED
#define OPENWQ_TD_HEADERH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>
#include "global/exprtk.hpp"
#include <cstdio>

#include "global/OpenWQ_vars.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_hostModelConfig.hpp"
#include "utils/headerfile_UTILS.hpp"
#include "output/headerfile_OUT.hpp"

#include "models_CH/headerfile_CH.hpp"
#include "models_TD/headerfile_td.hpp"


class OpenWQ_TD_model{

    public:

        // ##########################
        // MAIN DRIVER 
        // => Run TD_model model
        // ##########################

        void TD_driver_run(
            OpenWQ_wqconfig& OpenWQ_wqconfig, 
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_output& OpenWQ_output,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r, 
            const double wmass_source);


        // ##########################
        // MODEL OPTIONS 
        // ##########################

        // ##########################
        // NATIVE Advection
        void Adv(
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            double wflux_s2r, 
            double wmass_source);
    
        // ##########################
         // NATIVE Advection and Dispersion
        void AdvDisp(
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            double wflux_s2r, 
            double wmass_source);

        // ##########################
        // SPECIAL CASE => INPUT FROM EWF (linked to precipitation inputs) 
        // ##########################

        // INPUT FROM EXTERNAL FORCES
        void EWF_flux_driver_run(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig, 
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_utils& OpenWQ_utils,
            OpenWQ_output& OpenWQ_output,
            std::string source_EWF_name,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r);

        void EWF_flux_IN(
            OpenWQ_vars& OpenWQ_vars, 
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r, // water flux (m3/s)
            int chemi,              // chemical id
            double ewf_conc);
};

#endif