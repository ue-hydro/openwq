// Copyright 2026, Diogo Costa
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

#include "OpenWQ_hostModelConfig.hpp"


// Constructor
OpenWQ_hostModelconfig::OpenWQ_hostModelconfig() 
{

    // Water volumes from hostmodel
    this->cellid_to_wq = std::unique_ptr<
        std::vector<std::vector<std::vector<std::vector<std::string>>>>>(
            new std::vector<std::vector<std::vector<std::vector<std::string>>>>);

    // Water volumes from hostmodel
    this->waterVol_hydromodel = std::unique_ptr<
        std::vector<                            // compartments
        arma::Cube<                             // ix, iy, iz
        double>>>(new std::vector<arma::cube>); 

    // Dependencies from hostmodel 
    // (to be available for BGC)
    this->dependVar = std::unique_ptr<
        std::vector<                            // dependency variable
        arma::Cube<                             // ix, iy, iz
        double>>>(new std::vector<arma::cube>);

    this->dependVar_scalar = std::unique_ptr<
        std::vector<
        double>>(new std::vector<double>);
}

// Destructor
OpenWQ_hostModelconfig::~OpenWQ_hostModelconfig() 
{
}

/******** 
 * Methods 
 * *********/

// cellid_to_wq methods to map OpenWQ elements in the outputs
// set cellid_to_wq label
void OpenWQ_hostModelconfig::set_cellid_to_wqlabel(const std::string& label) {
    cellid_to_wqlabel = label;
}

// Const version - returns const reference to unique_ptr
std::string OpenWQ_hostModelconfig::get_cellid_to_wqlabel() const {
    return cellid_to_wqlabel;
}

void OpenWQ_hostModelconfig::set_cellid_to_wq_size(arma::Cube<double> domain_xyz) 
{
    std::vector<std::vector<std::vector<std::string>>> cube;
    cube.resize(domain_xyz.n_rows);
    for (arma::uword i = 0; i < domain_xyz.n_rows; ++i) {
        cube[i].resize(domain_xyz.n_cols);
        for (arma::uword j = 0; j < domain_xyz.n_cols; ++j) {
            cube[i][j].resize(domain_xyz.n_slices);
            for (arma::uword k = 0; k < domain_xyz.n_slices; ++k) {
                cube[i][j][k] = std::string();
            }
        }
    }
    (*this->cellid_to_wq).push_back(std::move(cube));
}

const std::unique_ptr<std::vector<std::vector<std::vector<std::vector<std::string>>>>>& OpenWQ_hostModelconfig::get_cellid_to_wq() const {
    return cellid_to_wq;
}

// Non-const version - returns non-const reference to unique_ptr
std::unique_ptr<std::vector<std::vector<std::vector<std::vector<std::string>>>>>& OpenWQ_hostModelconfig::get_cellid_to_wq() {
    return cellid_to_wq;
}

// set cellid_to_wq - ids fron hostmodel to map OpenWQ elements in the outputs
void OpenWQ_hostModelconfig::set_cellid_to_wq_at(int index, int ix, int iy, int iz, const std::string& value)
{
    (*this->cellid_to_wq)[index][ix][iy][iz] = value;
}
std::string OpenWQ_hostModelconfig::get_cellid_to_wq_at(int index, int ix, int iy, int iz)
{
    return (*this->cellid_to_wq)[index][ix][iy][iz];
}

// Reverse lookup: find (ix, iy, iz) from cell_id string
// Searches through the cellid_to_wq structure for the given compartment
// Returns true if found, false otherwise. If found, ix, iy, iz are set to the indices.
// partial_match output: set to true if only a prefix match was found (ix set, iy/iz not definitive)
bool OpenWQ_hostModelconfig::find_indices_from_cellid(int compartment_index,
                                                       const std::string& cell_id,
                                                       int& ix, int& iy, int& iz,
                                                       bool& partial_match)
{
    partial_match = false;

    // Check if cellid_to_wq is initialized and compartment index is valid
    if (!this->cellid_to_wq ||
        compartment_index < 0 ||
        compartment_index >= static_cast<int>((*this->cellid_to_wq).size())) {
        return false;
    }

    // Get the 3D structure for this compartment
    const auto& compartment_data = (*this->cellid_to_wq)[compartment_index];

    // First pass: try exact match
    for (size_t i = 0; i < compartment_data.size(); ++i) {
        for (size_t j = 0; j < compartment_data[i].size(); ++j) {
            for (size_t k = 0; k < compartment_data[i][j].size(); ++k) {
                if (compartment_data[i][j][k] == cell_id) {
                    // Found exact match - set all output indices
                    ix = static_cast<int>(i);
                    iy = static_cast<int>(j);
                    iz = static_cast<int>(k);
                    return true;
                }
            }
        }
    }

    // Second pass: try prefix match (e.g., cell_id "1" matches "1_z1", "1_z2", etc.)
    // This handles cases like SUMMA where cell_ids are "{hruId}_z{layer}"
    // but the SS JSON only specifies the hruId part
    std::string prefix = cell_id + "_";
    for (size_t i = 0; i < compartment_data.size(); ++i) {
        for (size_t j = 0; j < compartment_data[i].size(); ++j) {
            for (size_t k = 0; k < compartment_data[i][j].size(); ++k) {
                const std::string& stored = compartment_data[i][j][k];
                if (stored.substr(0, prefix.size()) == prefix) {
                    // Found prefix match - set ix only
                    // iy and iz should be handled by their own "all" or specific parsing
                    ix = static_cast<int>(i);
                    iy = static_cast<int>(j);
                    iz = static_cast<int>(k);
                    partial_match = true;
                    return true;
                }
            }
        }
    }

    // Not found
    return false;
}


// Single-dimension lookup: given a compartment and a known ix, find iy or iz
// by searching cell_ids that contain the search_str as a substring.
// dim: 1 = search along iy (vary j, fix i=known_ix), 2 = search along iz (vary k, fix i=known_ix)
// Returns the found dimension index in dim_result. Returns true if found.
bool OpenWQ_hostModelconfig::find_index_single_dim(
    int compartment_index,
    int known_ix,
    int dim,               // 1=iy, 2=iz
    const std::string& search_str,
    int& dim_result)
{
    if (!this->cellid_to_wq ||
        compartment_index < 0 ||
        compartment_index >= static_cast<int>((*this->cellid_to_wq).size())) {
        return false;
    }

    const auto& compartment_data = (*this->cellid_to_wq)[compartment_index];

    if (known_ix < 0 || known_ix >= static_cast<int>(compartment_data.size())) {
        return false;
    }

    const auto& ix_data = compartment_data[known_ix];

    if (dim == 1) {
        // Search along iy dimension
        for (size_t j = 0; j < ix_data.size(); ++j) {
            for (size_t k = 0; k < ix_data[j].size(); ++k) {
                const std::string& stored = ix_data[j][k];
                // Check if search_str matches as substring (e.g., "1" in "1_z1")
                if (stored == search_str || stored.find("_" + search_str) != std::string::npos
                    || stored.find(search_str + "_") != std::string::npos) {
                    dim_result = static_cast<int>(j);
                    return true;
                }
            }
        }
    } else if (dim == 2) {
        // Search along iz dimension
        for (size_t j = 0; j < ix_data.size(); ++j) {
            for (size_t k = 0; k < ix_data[j].size(); ++k) {
                const std::string& stored = ix_data[j][k];
                // Check for suffix match (e.g., "z1" or "1" matching "_z1" part of "1_z1")
                // Try exact suffix "_z{search_str}" first (SUMMA layer format)
                std::string z_suffix = "_z" + search_str;
                if (stored.size() >= z_suffix.size() &&
                    stored.compare(stored.size() - z_suffix.size(), z_suffix.size(), z_suffix) == 0) {
                    dim_result = static_cast<int>(k);
                    return true;
                }
                // Then try direct substring match
                if (stored == search_str || stored.find(search_str) != std::string::npos) {
                    dim_result = static_cast<int>(k);
                    return true;
                }
            }
        }
    }

    return false;
}

// Add a compartment to the vector of compartments
void OpenWQ_hostModelconfig::add_HydroComp(int index, std::string name, int num_cells_x, 
                    int num_cells_y, int num_cells_z) 
{
    this->HydroComp.push_back(OpenWQ_hostModelconfig::hydroTuple(
        index, name, num_cells_x, num_cells_y, num_cells_z));
}
// Add an external flux to the vector of external fluxes
void OpenWQ_hostModelconfig::add_HydroExtFlux(int index, std::string name, int num_cells_x, 
                    int num_cells_y, int num_cells_z) 
{
    this->HydroExtFlux.push_back(OpenWQ_hostModelconfig::hydroTuple(
        index, name, num_cells_x, num_cells_y, num_cells_z));
}
// Add a dependency to the vector of dependencies
void OpenWQ_hostModelconfig::add_HydroDepend(int index, std::string name, int num_cells_x, 
                    int num_cells_y, int num_cells_z) 
{
    this->HydroDepend.push_back(OpenWQ_hostModelconfig::hydroTuple(
        index, name, num_cells_x, num_cells_y, num_cells_z));
}

/**
 * HydroTuple methods
*/
unsigned int OpenWQ_hostModelconfig::get_num_HydroComp()
{
    return this->HydroComp.size();
}
unsigned int OpenWQ_hostModelconfig::get_num_HydroExtFlux()
{
    return this->HydroExtFlux.size();
}
unsigned int OpenWQ_hostModelconfig::get_num_HydroDepend()
{
    return this->HydroDepend.size();
}

// get names of compartments and external fluxes
// for a specific index
std::string OpenWQ_hostModelconfig::get_HydroComp_name_at(int index)
{
    return std::get<1>(this->HydroComp[index]);
}
std::string OpenWQ_hostModelconfig::get_HydroDepend_name_at(int index)
{
    return std::get<1>(this->HydroDepend[index]);
}

// get index of compartment name
int OpenWQ_hostModelconfig::get_HydroComp_index(
    std::string cmp_name,
    std::string msg_string,
    bool abort_if_noexist){

    int icmp = 0;

    for (auto &i : this->HydroComp)
    {
        if (cmp_name.compare(std::get<1>(i))==0)
            break;
        icmp += 1;
    }

    // if not found, set to -1
    if (icmp == get_num_HydroComp()){
        icmp = -1;
    }

    // abort if compartment not found
    if (icmp == -1 && abort_if_noexist){

        msg_string += "> Could not find compartment: " + cmp_name;

        // Print in console 
        std::cout << msg_string << std::endl;

        // Abort
        exit(EXIT_FAILURE);
    }

    return icmp;

}

// Get vector of compartment names and external fluxes names
std::vector<std::string> OpenWQ_hostModelconfig::get_HydroComp_names()
{
    std::vector<std::string> names;
    for (auto &i : this->HydroComp)
    {
        names.push_back(std::get<1>(i));
    }
    return names;
}
std::vector<std::string> OpenWQ_hostModelconfig::get_HydroExtFlux_names()
{
    std::vector<std::string> names;
    for (auto &i : this->HydroExtFlux)
    {
        names.push_back(std::get<1>(i));
    }
    return names;
}

// Get number of cells in x, y and z directions - HydroComp
unsigned int OpenWQ_hostModelconfig::get_HydroComp_num_cells_x_at(int index)
{
    return std::get<2>(this->HydroComp[index]);
}
unsigned int OpenWQ_hostModelconfig::get_HydroComp_num_cells_y_at(int index)
{
    return std::get<3>(this->HydroComp[index]);
}
unsigned int OpenWQ_hostModelconfig::get_HydroComp_num_cells_z_at(int index)
{
    return std::get<4>(this->HydroComp[index]);
}

// Get number of cells in x, y and z directions - HydroExtFlux
unsigned int OpenWQ_hostModelconfig::get_HydroExtFlux_num_cells_x_at(int index)
{
    return std::get<2>(this->HydroExtFlux[index]);
}
unsigned int OpenWQ_hostModelconfig::get_HydroExtFlux_num_cells_y_at(int index)
{
    return std::get<3>(this->HydroExtFlux[index]);
}
unsigned int OpenWQ_hostModelconfig::get_HydroExtFlux_num_cells_z_at(int index)
{
    return std::get<4>(this->HydroExtFlux[index]);
}

// Get number of cells in x, y and z directions - HydroDepend
unsigned int OpenWQ_hostModelconfig::get_HydroDepend_num_cells_x_at(int index)
{
    return std::get<2>(this->HydroDepend[index]);
}
unsigned int OpenWQ_hostModelconfig::get_HydroDepend_num_cells_y_at(int index)
{
    return std::get<3>(this->HydroDepend[index]);
}
unsigned int OpenWQ_hostModelconfig::get_HydroDepend_num_cells_z_at(int index)
{
    return std::get<4>(this->HydroDepend[index]);
}

/*
* waterVol methods
*/
void OpenWQ_hostModelconfig::add_waterVol_hydromodel(arma::Cube<double> waterVol) 
{
    (*this->waterVol_hydromodel).push_back(waterVol);
}
double OpenWQ_hostModelconfig::get_waterVol_hydromodel_at(int index, int ix, int iy, int iz) 
{
    return (*this->waterVol_hydromodel)[index](ix,iy,iz);
}
void OpenWQ_hostModelconfig::set_waterVol_hydromodel_at(int index, int ix, int iy, int iz, double value) 
{
    (*this->waterVol_hydromodel)[index](ix,iy,iz) = value;
}
// waterVol_minlim method
double OpenWQ_hostModelconfig::get_watervol_minlim()
{
    return this->watervol_minlim;
}
// waterflux_minlim method
double OpenWQ_hostModelconfig::get_waterflux_minlim()
{
    return this->waterflux_minlim;
}

/*
* dependVar methods
*/
void OpenWQ_hostModelconfig::add_dependVar(arma::Cube<double> dependVar) 
{
    (*this->dependVar).push_back(dependVar);
}
double OpenWQ_hostModelconfig::get_dependVar_at(int index, int ix, int iy, int iz) 
{
    return (*this->dependVar)[index](ix,iy,iz);
}
void OpenWQ_hostModelconfig::set_dependVar_at(int index, int ix, int iy, int iz, double value) 
{
    (*this->dependVar)[index](ix,iy,iz) = value;
}
// dependVar_scalar methods
void OpenWQ_hostModelconfig::add_dependVar_scalar(double dependVar_scalar)
{
    (*this->dependVar_scalar).push_back(dependVar_scalar);
}
double OpenWQ_hostModelconfig::get_dependVar_scalar_at(int index)
{
    return (*this->dependVar_scalar)[index];
}
void OpenWQ_hostModelconfig::set_dependVar_scalar_at(int index, double value)
{
    (*this->dependVar_scalar)[index] = value;
}


/**************
 interaction_step methods
**************/
void OpenWQ_hostModelconfig::increment_interaction_step()
{
    this->interaction_step++;
}

bool OpenWQ_hostModelconfig::is_first_interaction_step()
{
    return (this->interaction_step == 1);
}


/**************
time methods
**************/
bool OpenWQ_hostModelconfig::is_time_step_0()
{
    return (this->time_step == 0);
}
void OpenWQ_hostModelconfig::set_time_step(double time_step)
{
    this->time_step = time_step;
}
double OpenWQ_hostModelconfig::get_time_step()
{
    return this->time_step;
}

