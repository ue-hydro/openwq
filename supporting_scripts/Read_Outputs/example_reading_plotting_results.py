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
        'mapping_key': 'reachID',
        'var2print': 'DWroutedRunoff' # 'basRunoff' or 'DWroutedRunoff
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
#####

# mapping results

import Map_h5_driver as h5_mplib

h5_mplib.Map_h5_driver(
    openwq_results=openwq_results,
    shpfile_fullpath_mapKey=shpfile_fullpath_mapKey,
    hydromodel_out_fullpath=hydromodel_out_fullpath,
    output_html_path="/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/mapResults.gif",
    chemSpec=["NO3-N"],
    hostmodel='mizuroute',  # mizuroute or summa
    what2map='openwq',  # hostmodel or openwq
    create_gif=True,
    timeframes=200,
    gif_duration=100
)


"""
# test
import read_ncdf as read_ncdf

shpfile_fullpath = "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/10_Univ_Evora/3_research/03_ACTIVE/1_MAIN_AUTHOR/9_openWQ/2_initial_develop_synthetic_tests_GMD_paper/code/mizuRoute/case_studies/synthetic_tests/2_nrTrans_instS_PorMedia/mizuroute/mizuroute_in/shapefile/shapefiles/Flowline_CO_14_cameo.shp"
mizuroute_out_fullpath = '/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/10_Univ_Evora/3_research/03_ACTIVE/1_MAIN_AUTHOR/9_openWQ/2_initial_develop_synthetic_tests_GMD_paper/code/mizuRoute/case_studies/synthetic_tests/4_nrTrans_contS_PorMedia/mizuroute/mizuroute_out/v1.2_case1.h.1950-01-02-43200.nc'
folderSaveFigs_fullpath = "/Users/diogocosta/Downloads/"
gif_name_fullpath = 'results_animation.html'
print_step = 10

read_ncdf.MapGeoPandas(shpfile_fullpath,
                openwq_results,
                mizuroute_out_fullpath,
                folderSaveFigs_fullpath,
                gif_name_fullpath,
                print_step)

"""