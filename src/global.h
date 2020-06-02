// Copyright 2020, Diogo Costa
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

#ifndef GLOBALH_INCLUDED
#define GLOBALH_INCLUDED

#include <string>
#include <armadillo>
#include <memory> 
#include "jnlohmann/json.h"
using json = nlohmann::json;


// Project Information
class JSONfiles
{
    public:

    json Master;
    json H2O;
    json BGC;

};


// General information about the project
class Prj_StateVar
{
    public:
    Prj_StateVar(){

    }
    Prj_StateVar(size_t numcmp){

        this-> numcmp = numcmp;

        try{
            char cmpt_names[numcmp]; // compartment names

            wflux = std::unique_ptr<arma::field<arma::Cube<double>>>(new arma::field<arma::cube>(numcmp)); // 1 field: flow
            wmass = std::unique_ptr<arma::field<arma::Cube<double>>>(new arma::field<arma::cube>(numcmp)); // 1 field: water mass
            chemass = std::unique_ptr<arma::field<arma::field<arma::Cube<double>>>>(new arma::field<arma::field<arma::cube>>(numcmp));  // multiple fields: one for eacg chem

        }catch(int e){
            std::cout << "An exception occurred creating the domain: ERR " << std::to_string(e) << std::endl;
        }

    }
    size_t numcmp;

    std::unique_ptr<arma::field<arma::Cube<double>>> wflux, wmass;
    std::unique_ptr<arma::field<arma::field<arma::Cube<double>>>> chemass; 

};

#endif