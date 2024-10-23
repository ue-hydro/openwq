
// Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
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


#ifndef OPENWQ_SOLVERH_INCLUDED
#define OPENWQ_SOLVERH_INCLUDED

#include "global/OpenWQ_vars.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_hostModelconfig.hpp"
#include "global/OpenWQ_json.hpp"
class OpenWQ_output;
class OpenWQ_chem;

class OpenWQ_solver{

    private:
    
    // Solver using Sundials
    void Solve_with_CVode(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_json& OpenWQ_json,
        OpenWQ_output& OpenWQ_output,
        OpenWQ_chem& OpenWQ_chem);

    // Solver using Euler method
    void Solve_with_BE(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_json& OpenWQ_json,
        OpenWQ_output& OpenWQ_output,
        OpenWQ_chem& OpenWQ_chem);

    public:

    // Generic Numerical Solver
    void Numerical_Solver(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_json& OpenWQ_json,
        OpenWQ_output& OpenWQ_output,
        OpenWQ_chem& OpenWQ_chem);

    // Reset derivatives (before each time iteraction)
    void Reset_Deriv(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars,
        bool inst_deriv_flag,
        bool cum_deriv_flag);

    // Reset EWF conc (before each time iteraction)
    void Reset_EWFconc(
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_vars& OpenWQ_vars);

};

#endif