
# Copyright 2020, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


##########################################
# Example on how to generate SS json from data in CSV file
##########################################

import sys
sys.path.insert(0, 'SS_lib')

import Generate_SS_Files_from_CSVdata as csv_sslib

# Metadata to go on JSON file
METADATA_Comment = 'Test'
METADATA_Source = 'Test'

# Configuration and data sources
chemNames = ["NO3", "NH4"]
comptNames = ["SWE", "SWE"]
ssType = ["source", "source"] # sink or source
ssUnits = ["kg", "kg"]
Data_CSV_sources = ["examples/NO3_fertilizer.csv", "examples/NH4_fertilizer.csv"]

# full path to export ss json file
new_SS_JSONfile = 'examples/OpenWQ_source_example.json'

csv_sslib.GenSSjson_fromCSV_driver(
    new_SS_JSONfile=new_SS_JSONfile,
    METADATA_Comment=METADATA_Comment,
    METADATA_Source=METADATA_Source,
    chemNames=chemNames,
    comptNames=comptNames,
    ssType=ssType,
    ssUnits=ssUnits,
    Data_CSV_sources=Data_CSV_sources
    )

##########################################
# Example on how to generate SS json from copernicus LULC
##########################################

import sys
sys.path.insert(0, 'SS_lib')

import Generate_SS_Files_from_LULC_Copernicus as copernicus_sslib

# Input data
hru_shp = "/Users/diogocosta/Documents/esa_cci_landcover/basins"

# Copernicus raster datasets (all years)
copernicus_lulc_dir = "/Users/diogocosta/Documents/esa_cci_landcover"

copernicus_sslib.GenSSjson_fromCopernicus_driver(
    hru_shp=hru_shp,
    copernicus_lulc_dir=copernicus_lulc_dir
)