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

openwq_mapKey_dic= {
    "path_to_results": '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out',
    "mapping_key": "reachID"}

shpfile_fullpath_mapKey={
        'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/shapefiles/mizuSegId/mizuSegId.shp',
        'mapping_key': 'SegId'
    }

hydromodel_out_fullpath={
        'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_out/allvarsL5.h.1961-01-01-00000.nc',
        'mapping_key': 'reachID' # 'basRunoff' or 'DWroutedRunoff
    }

# Read requested data: Read_h5_driver
openwq_results = h5_rlib.Read_h5_driver(
            openwq_mapKey_dic= openwq_mapKey_dic,
            output_format='HDF5',
            debugmode=False,
            cmp=['RIVER_NETWORK_REACHES'],
            space_elem='all',
            chemSpec=["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"],
            chemUnits="MG/L",
            noDataFlag=-9999)

######
# Creating GIF with results
# mapping results
#####

import Map_h5_driver as h5_mplib

h5_mplib.Map_h5_driver(
    # What to print?
    what2map='openwq',  # hostmodel or openwq
    hostmodel='mizuroute',  # mizuroute or summa
    # shapefile info
    shpfile_fullpath_mapKey=shpfile_fullpath_mapKey,
    # openwq info (👉🏼 used if what2map=openwq)
    openwq_results=openwq_results,
    chemSpec=["NO3-N"],
    # hostmodel info (👉🏼 used if what2map=hostmodel)
    hydromodel_out_fullpath=hydromodel_out_fullpath,
    hydromodel_var2print='DWroutedRunoff',
    # gif config
    output_html_path="/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/mapResults.gif",
    create_gif=True,
    timeframes=200,
    gif_duration=100
)

"""
import Plot_h5_driver as h5_mplib

h5_mplib.Plot_h5_driver(
    hydromodel_out_fullpath={
        'path_to_shp': 'summa_output.nc',
        'mapping_key': 'hruId',
        'var2print': 'scalarSWE',
        'feature_ids': [1, 2, 3, 4, 5]
    },
    output_path='/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/plotSeries.png',
    hostmodel='mizuroute',
    what2map='hostmodel')
"""