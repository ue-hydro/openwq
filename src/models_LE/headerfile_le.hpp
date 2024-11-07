

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


#ifndef OPENWQ_HEADERFILE_LEH_INCLUDED
#define OPENWQ_HEADERFILE_LEH_INCLUDED

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

class OpenWQ_LE_model{

    public:

        // Boundary Mixing due to velocity gradients
        // due to turbulence and cross-boarder eddies
        void native_BoundMix(
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            const int source, const int ix_s, const int iy_s, const int iz_s,
            const int recipient, const int ix_r, const int iy_r, const int iz_r,
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