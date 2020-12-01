
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


#include "DEMOS_OpenWQ_print.h"

int DEMOS_OpenWQ_print::writeVTU(DEMOS_OpenWQ_json& DEMOS_OpenWQ_json,
    DEMOS_OpenWQ_vars& DEMOS_OpenWQ_vars,
    int numcmp, 
    int tmpst_i)
{

    // This was built for hexahedrons (see line 66) for now (cubes, Rectangular cuboid, Trigonal trapezohedron, etc)
    // But line 33:33 makes determining cubes of side lenght = 1 (but this can all be changed)
    // based on https://lorensen.github.io/VTKExamples/site/Cxx/IO/WriteVTU/

    // loop over compartments
    for (int icmp=0;icmp<numcmp;icmp++){

        std::string filename = DEMOS_OpenWQ_json.Master["export_results_folder"];
        filename.append("/");
        filename.append(DEMOS_OpenWQ_json.H2O["compartments"][std::to_string(icmp+1)]);
        filename.append("_");
        filename.append(std::to_string(tmpst_i)); // time stamp
        filename.append(".vtu");
        std::string chemname;

        // get domain dimensions
        int nx = DEMOS_OpenWQ_json.H2O[std::to_string(icmp+1)]["nx"];
        int ny = DEMOS_OpenWQ_json.H2O[std::to_string(icmp+1)]["ny"];
        int nz = DEMOS_OpenWQ_json.H2O[std::to_string(icmp+1)]["nz"];

        // get number of species for compartment icmp
        int numspec = DEMOS_OpenWQ_json.WQ["compartments"][std::to_string(icmp+1)]["chem_species"].size();

        // total number of vertices and hexahedrons
        int numvert = (nx+1)*(ny+1)*(nz+1);
        int nnumel = nx * ny * nz;
        
        // Determine the all the vertices (assuming a spacing of 1 m for now)
        double x[numvert][3];
        int i = 0;
        for (int iz=0;iz<=nz;iz++){   
                for (int ix=0;ix<=nx;ix++){
                    for (int iy=0;iy<=ny;iy++){
                    x[i][0] = ix-0.5;
                    x[i][1] = iy-0.5;
                    x[i][2] = iz-0.5;
                    i++;
                }
            }
        }
    
    // Determine the faces of each element/hexahedron
        vtkIdType pts[nnumel][8];
        i = 0;
        for (int iz=0;iz<nz;iz++){   
                for (int ix=0;ix<nx;ix++){
                    for (int iy=0;iy<ny;iy++){
                    pts[i][0] = (nx+1)*(ny+1)*(iz) + (ny+1)*(ix) + iy;
                    pts[i][1] = (nx+1)*(ny+1)*(iz) + (ny+1)*(ix+1) + iy;
                    pts[i][2] = (nx+1)*(ny+1)*(iz) + (ny+1)*(ix+1) + iy+1;
                    pts[i][3] = (nx+1)*(ny+1)*(iz) + (ny+1)*(ix) + iy+1;
                    pts[i][4] = (nx+1)*(ny+1)*(iz+1) + (ny+1)*(ix) + iy;
                    pts[i][5] = (nx+1)*(ny+1)*(iz+1) + (ny+1)*(ix+1) + iy;
                    pts[i][6] = (nx+1)*(ny+1)*(iz+1) + (ny+1)*(ix+1) + iy+1;
                    pts[i][7] = (nx+1)*(ny+1)*(iz+1) + (ny+1)*(ix) + iy+1;
                i++;
                }
            }
        }

        // Determine vtkPoints for vtkUnstructuredGrid build
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for (i=0; i<numvert; i++) points->InsertPoint(i,x[i]);

        // Initiate vtkUnstructuredGrid and add all the hexahedrons
        vtkSmartPointer<vtkUnstructuredGrid> ugrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
        ugrid->Allocate(100);
        for (i=0; i<nnumel; i++){
            ugrid->InsertNextCell(VTK_HEXAHEDRON, 8, pts[i]);
        }

        ugrid->SetPoints(points);

        // Add information to the unstructured grid: compartment data
    
        std::vector<std::string> chemnames = 
            DEMOS_OpenWQ_json.WQ["compartments"][std::to_string(icmp+1)]["chem_species"]; //chem species # within compartment icmp (from JSON.BGQ)

        for (int ichem=0;ichem<numspec;ichem++){ // all chemical species

            vtkSmartPointer<vtkDoubleArray> varexpot = vtkSmartPointer<vtkDoubleArray>::New();
            varexpot->SetNumberOfValues(numvert);

            // Set name of array (chem variable name)
            chemname = chemnames[ichem];
            varexpot->SetName(chemname.c_str());

            i = 0;
            for (int iz=0;iz<=nz;iz++){   
                    for (int ix=0;ix<=nx;ix++){
                        for (int iy=0;iy<=ny;iy++){
                            if(iz!=nz && iy!=ny && ix!=nx){
                                varexpot->SetValue(i, (*DEMOS_OpenWQ_vars.chemass)(icmp)(ichem)(ix,iy,iz));
                            }else{
                                varexpot->SetValue(i, 0);
                            }
                    i++;
                    }
                }
            }
            ugrid->GetPointData()->AddArray(varexpot);
        }

        // Write file
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();  
        writer->SetFileName(filename.c_str());
        
        writer->SetInputData(ugrid);
        writer->Write();

    }

    return EXIT_SUCCESS;

}