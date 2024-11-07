

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


#ifndef OPENWQ_TD_HEADERH_INCLUDED
#define OPENWQ_TD_HEADERH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>
#include "global/exprtk.hpp"
#include <cstdio>
//#include "utility.h"
//#include "ADE_solver.h"

#include "global/openwq_vars.hpp"
#include "global/openwq_wqconfig.hpp"
#include "global/openwq_hostmodelconfig.hpp"
//#include "OpenWQ_initiate.h"
//#include "OpenWQ_print.h"

class OpenWQ_TD_model{

    public:

        // Advection and Dispersion
        void AdvDisp(
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            double wflux_s2r, 
            double wmass_source);

        // Advection and Dispersion
        // General case (flux exchanges within the model domain)
        void Adv(
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            double wflux_s2r, 
            double wmass_source);

        // Special case (in-fluxes from external water flux sources)
        void Adv_IN(
            OpenWQ_vars& OpenWQ_vars, 
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
            const double wflux_s2r,     // water flux (m3/s)
            int chemi,                  // chemical id
            double ewf_conc);           // concentration of chemical

};

#endif