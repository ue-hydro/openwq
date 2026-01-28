#####
# Reading Openwq outputs and plot
#####

import sys
import os
sys.path.insert(0, 'hdf5_support_lib')

import Read_h5_driver as h5_rlib

######################
# General Mapping Info
######################

shpfile_info={
        'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/shapefiles/mizuSegId/mizuSegId.shp',
        'mapping_key': 'SegId'
    }

hydromodel_info={
        'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_out/allvarsL5.h.1961-01-01-00000.nc',
        'mapping_key': 'reachID' # 'basRunoff' or 'DWroutedRunoff
    }

openwq_info= {
    "path_to_results": '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out',
    "mapping_key": "reachID"}

# Read requested data: Read_h5_driver
openwq_results = h5_rlib.Read_h5_driver(
            openwq_info= openwq_info,
            output_format='HDF5',
            debugmode=False,
            cmp=['RIVER_NETWORK_REACHES'],
            space_elem='all',
            chemSpec=["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"],
            chemUnits="MG/L",
            noDataFlag=-9999)


##############################
# Creating GIF MAPS with results
#############################

"""
import Map_h5_driver as h5_mplib

h5_mplib.Map_h5_driver(
    # 1) What results to map?
    what2map='openwq',  # hostmodel or openwq
    hostmodel='mizuroute',  # mizuroute or summa
    # 2) shapefile info
    shpfile_info=shpfile_info,
    # 3) openwq info (👉🏼 used if what2map=openwq)
    openwq_results=openwq_results,
    chemSpec=["NO3-N"],
    # 4) hostmodel info (👉🏼 used if what2map=hostmodel)
    hydromodel_info=hydromodel_info,
    hydromodel_var2print='DWroutedRunoff',
    # 5) output config
    output_html_path="/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/mapResults.gif",
    create_gif=True,
    timeframes=200,
    gif_duration=100
)
"""

##############################
# Creating GIF MAPS with results
#############################

import Plot_h5_driver as h5_mplib

h5_mplib.Plot_h5_driver(
    # 1) What results to map?
    what2map='openwq',
    hostmodel='mizuroute',
    mapping_key_values=[1200014181, 200014181],
    # 2) openwq info (👉🏼 used if what2map=openwq)
    openwq_results=openwq_results,
    chemSpec=["NO3-N"],
    # 3) hostmodel info (👉🏼 used if what2map=hostmodel)
    hydromodel_info=hydromodel_info,
    hydromodel_var2print='DWroutedRunoff',
    # 4) output config
    output_path='/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/plotSeries.png'
    )
