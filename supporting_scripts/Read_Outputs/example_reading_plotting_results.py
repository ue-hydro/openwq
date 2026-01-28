#####
# Reading Openwq outputs and plot
#####

import sys
import os
sys.path.insert(0, 'hdf5_support_lib')


###############
# CONFIGURATION
##############

# Openwq output folder
openwq_resDir_mapKey = {
    "path_to_results": '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out',
    "mapping_key": "reachID"}

# Output specs
openwq_output_format = 'HDF5' # needs to be HDF5 (CSV disabled)
openwq_output_debugmode = False

# What to print?
chemical_species = ["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"]
concentration_units = "MG/L"
compartments = ['RIVER_NETWORK_REACHES']
cells_spatial = 'all'

# Where to save plots?
plots_save_dir = openwq_resDir_mapKey["path_to_results"]

"""
# SUMMA
# Openwq output folder
openwq_resDir = '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out'

# Output specs
openwq_output_format = 'HDF5' # needs to be HDF5 (CSV disabled)
openwq_output_debugmode = False

# What to print?
chemical_species = ["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"]
concentration_units = "MG/L"
compartments = ['RIVER_NETWORK_REACHES']
cells_spatial = 'all'

# Where to save plots?
plots_save_dir = '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out'
"""

###############
# CALLING THE NECESSARY SUPPORTING SCRIPS
# 1) read
# 2) plot or/and map
##############

import Read_h5_driver as h5_rlib

# Read requested data: Read_h5_driver
openwq_results = h5_rlib.Read_h5_driver(
            resDir_mapKey_dic=openwq_resDir_mapKey,
            output_format=openwq_output_format,
            debugmode=openwq_output_debugmode,
            cmp=compartments,
            space_elem=cells_spatial,
            chemSpec=chemical_species,
            chemUnits=concentration_units,
            noDataFlag=-9999)

######
# Plot results
#####


import Plot_h5_driver as h5_plib

# Example 1: Temporal plot using coordinates
"""
figs1 = h5_plib.Plot_h5_driver(
        openwq_results,
        plot_type='temporal',
        cells_to_plot=[(0, 0, 0), (100, 0, 0), (200, 0, 0)],
        output_dir=plots_save_dir,
        debug_mode=False
    )
"""
"""
# Example 2: Spatial plot - snapshot at specific date
figs4 = h5_plib.Plot_h5_driver(
        openwq_results,
        plot_type='spatial',
        time_slice=50,
        # time_slice=15
        spatial_axis='x',
        fixed_coords={'y': 0, 'z': 0},
        output_dir=plots_save_dir,
        debug_mode=False
    )
"""

# mapping results

import Map_h5_driver as h5_mplib

shpfile_fullpath_mapKey = {
    "path_to_shp": "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/shapefiles/mizuSegId/mizuSegId.shp",
    "mapping_key": "SegId"}
hydromodel_out_fullpath =  "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_out/allvarsL5.h.1961-01-01-00000.nc"
output_html_path = plots_save_dir + "/mapResults.gif"
hostmodel="mizuroute"
timestep = 40
what2map='hostmodel' # hostmodel or openwq

h5_mplib.Map_h5_driver(
    shpfile_fullpath_mapKey=shpfile_fullpath_mapKey,
    openwq_results=openwq_results,
    hydromodel_out_fullpath=hydromodel_out_fullpath,
    output_html_path=output_html_path,
    chemical_species=chemical_species[1],
    file_extension='main',
    timesteps=timestep,
    hostmodel=hostmodel,
    what2map="openwq",
    create_gif=True,
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