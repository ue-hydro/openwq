# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
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


#################################################
# Example on how to use read, plot and map model results
# using the supporting libraries in "hdf5_support_lib"
#
# Set sediment_as_well=True to also read, plot and map
# sediment transport results (requires export_sediment=True
# in template_model_config.py and TS module enabled).
#################################################
# Loading libraries
import sys
sys.path.insert(0, 'hdf5_support_lib')

############################################
# General Mapping Info:
# 1) shapefile
# 2) hydromodel results
# 3) openwq results

shpfile_info={
        'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2/mizuroute_in/shapefiles/mizuSegId/mizuSegId.shp',
        'mapping_key': 'SegId'
    }

hydromodel_info={
        'path_to_results': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2/mizuroute_out/allvarsL5.h.1961-01-01-00000.nc',
        'mapping_key': 'reachID'
    }

openwq_info= {
    "path_to_results": '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2/openwq_out',
    "mapping_key": "reachID"
    }

######################
# Read requested data: Read_h5_driver
# sediment_as_well=True will also read the COMPARTMENT@Sediment-main.h5 file
#######################

import Read_h5_driver as h5_rlib

openwq_results = h5_rlib.Read_h5_driver(
            openwq_info= openwq_info,
            output_format='HDF5',   # don't change
            debugmode=False,         # if True, will read also d_output_dt_chemistry, d_output_dt_transport, d_output_ss, d_output_ewf, d_output_ic
            cmp=['RIVER_NETWORK_REACHES'],
            space_elem='all',       # what cells to read
            chemSpec=["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"],
            chemUnits="MG/L",
            noDataFlag=-9999,
            sediment_as_well=True)  # also read sediment transport HDF5 output

############################################
# Creating GIF MAPS with results
# sediment_as_well=True will also map sediment in addition to the chemSpec listed
"""
import Map_h5_driver as h5_mplib

h5_mplib.Map_h5_driver(
    # 1) What results to map?
    what2map='openwq',  # hostmodel or openwq
    hostmodel='mizuroute',  # mizuroute or summa
    # 2) shapefile info
    shpfile_info=shpfile_info,
    # 3) openwq info (used if what2map=openwq)
    openwq_results=openwq_results,
    chemSpec=["NO3-N","N_ORG_fresh"],
    sediment_as_well=True,  # also map sediment transport
    # 4) hostmodel info (used if what2map=hostmodel)
    hydromodel_info=hydromodel_info,
    hydromodel_var2print='DWroutedRunoff', #e.g., 'basRunoff' or 'DWroutedRunoff
    # 5) output config
    output_html_path="/Users/diogocosta/Documents/openwq_code/test_demonstration/mizuroute_cslm_openwq/route/build/openwq/openwq/bin/openwq_out/",
    create_gif=True,
    timeframes=30,
    gif_duration=30
)
"""

############################################
# Plotting time evolution at specific locations
# sediment_as_well=True will also plot sediment in addition to the chemSpec listed

import Plot_h5_driver as h5_plib

h5_plib.Plot_h5_driver(
    # 1) What results to plot?
    what2map='openwq',  # hostmodel or openwq
    hostmodel='mizuroute',
    mapping_key_values=[1200014181, 200014181],
    # 2) openwq info (used if what2map=openwq)
    openwq_results=openwq_results,
    chemSpec=["NO3-N","N_ORG_active"],
    sediment_as_well=True,  # also plot sediment transport
    debugmode=False, # if True, will read also d_output_dt_chemistry, d_output_dt_transport, d_output_ss, d_output_ewf, d_output_ic
    # 3) hostmodel info (used if what2map=hostmodel)
    hydromodel_info=hydromodel_info,
    hydromodel_var2print='DWroutedRunoff',
    # 4) output config
    output_path='/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/bin/openwq_out/plotSeries.png'
    )
