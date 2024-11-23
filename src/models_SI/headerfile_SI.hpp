
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


#ifndef OPENWQ_SI_HEADERFILEH_INCLUDED
#define OPENWQ_SI_HEADERFILEH_INCLUDED

#include <armadillo>
#include <string>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <vector>

// #include "utility.h"
#include "global/openwq_json.hpp"

#include "global/openwq_json.hpp"
#include "global/openwq_vars.hpp"
#include "global/openwq_wqconfig.hpp"
#include "global/openwq_hostmodelconfig.hpp"
#include "output/headerfile_OUT.hpp"
#include "units/headerfile_UNITS.hpp"

// Biogeochemistry

class OpenWQ_SI_model{

    public:

        // ##########################
        // MAIN DRIVER 
        // => Run ST model
        // ##########################

        void SI_driver_run(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_output& OpenWQ_output);

    private:

        // ##########################
        // MODEL OPTIONS 
        // ##########################

        // langmuir
        void langmuir(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_output& OpenWQ_output);

        // freundlich
        void freundlich(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_output& OpenWQ_output);

    private:

        // ##################
        // BGC_FLEX model
        // ##################

        // Perform transformations
        // chem_runtimestart.cpp
        void bgc_flex_transform(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_output& OpenWQ_output,
            unsigned int icmp);

        // ##################
        // PHREEQC model
        // ##################

        // add here

};

#endif