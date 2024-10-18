
// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
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

#ifndef OPENWQ_OUTPUT_INCLUDED
#define OPENWQ_OUTPUT_INCLUDED

#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_vars.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_hostModelconfig.hpp"
#include "solver/OpenWQ_solver.hpp"

class OpenWQ_output{

    public:

        // main print function
        int writeResults(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_vars& OpenWQ_vars,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            OpenWQ_solver& OpenWQ_solver,
            time_t simtime);
        
        // Print in console and log file
        void ConsoleLog(
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            std::string& msg_string,
            bool print_console,
            bool print_logFile);
    private:
        // Print output in CSV
        int writeCSV(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            std::unique_ptr<arma::field<arma::field<arma::cube>>>& OpenWQ_var2print,
            std::string& output_file_label,
            std::string timestr,            // time step (in seconds)
            int icmp);
        
    // Print output in HDF5
        int writeHDF5(
            OpenWQ_json& OpenWQ_json,
            OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
            OpenWQ_wqconfig& OpenWQ_wqconfig,
            std::unique_ptr<arma::field<arma::field<arma::cube>>>& OpenWQ_var2print,
            std::string& output_file_label,
            std::string timestr,            // time step (in seconds)
            int icmp);

    // append data to HDF5 file 
        bool appendData_to_HDF5_file(
            hid_t file,
            arma::mat& data,
            std::string name);

};

#endif