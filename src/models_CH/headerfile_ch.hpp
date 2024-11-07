
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


#ifndef OPENWQ_CHEM_HEADERFILEH_INCLUDED
#define OPENWQ_CHEM_HEADERFILEH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>
#include "global/exprtk.hpp"
#include <cstdio>
#include <tuple>
#include <vector>

// #include "utility.h"
#include "global/openwq_json.hpp"
using json = nlohmann::json;

#include "global/openwq_json.hpp"
#include "global/openwq_vars.hpp"
#include "global/openwq_wqconfig.hpp"
#include "global/openwq_hostmodelconfig.hpp"
#include "output/headerfile_output.hpp"
#include "units/headerfile_units.hpp"

// Biogeochemistry

class OpenWQ_chem{

    public:
        // Parse biogeochemical expressions
        // chem_decl.cpp
        void setBGCexpressions(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_units& OpenWQ_units,
            OpenWQ_output& OpenWQ_output);

        // Run Chemistry (main call)
        // chem_runtimestart.cpp
        void Run(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_output& OpenWQ_output);
    
    private:
        // Perform transformations
        // chem_runtimestart.cpp
        void BGC_Transform(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_output& OpenWQ_output,
            unsigned int icmp);

};

#endif