#####
# Reading Openwq outputs and plot
#####

import sys
import os
sys.path.insert(0, 'hdf5_support_lib')


###############
# CONFIGURATION
##############

"""
# Openwq output folder
openwq_resDir = '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out'

# Output specs
OpenWQ_output_format = 'HDF5' # needs to be HDF5 (CSV disabled)
OpenWQ_output_debugmode = False

# What to print?
Chemical_species = ["SPECIES_A","SPECIES_B"]
Concentration_units = "MG/L"
Compartments = ['RIVER_NETWORK_REACHES']
Cells_spatial = 'all'

# Where to save plots?
plots_save_dir = '/Users/diogocosta/Downloads/'
"""

# Openwq output folder
openwq_resDir = '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out'

# Output specs
OpenWQ_output_format = 'HDF5' # needs to be HDF5 (CSV disabled)
OpenWQ_output_debugmode = False

# What to print?
Chemical_species = ["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"]
Concentration_units = "MG/L"
Compartments = ['RIVER_NETWORK_REACHES']
Cells_spatial = 'all'

# Where to save plots?
plots_save_dir = '/Users/diogocosta/Downloads/'

"""
# SUMMA
# Openwq output folder
openwq_resDir = '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out'

# Output specs
OpenWQ_output_format = 'HDF5' # needs to be HDF5 (CSV disabled)
OpenWQ_output_debugmode = False

# What to print?
Chemical_species = ["NO3-N","NH4-N","N_ORG_fresh","N_ORG_stable","N_ORG_active"]
Concentration_units = "MG/L"
Compartments = ['RIVER_NETWORK_REACHES']
Cells_spatial = 'all'

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
            openwq_resDir,
            OpenWQ_output_format,
            OpenWQ_output_debugmode,
            Compartments,
            Cells_spatial,
            Chemical_species,
            Concentration_units,
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

shpfile_fullpath = "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/10_Univ_Evora/3_research/03_ACTIVE/1_MAIN_AUTHOR/9_openWQ/2_initial_develop_synthetic_tests_GMD_paper/code/mizuRoute/case_studies/synthetic_tests/2_nrTrans_instS_PorMedia/mizuroute/mizuroute_in/shapefile/shapefiles/Flowline_CO_14_cameo.shp"
hydromodel_out_fullpath =  "/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/10_Univ_Evora/3_research/03_ACTIVE/1_MAIN_AUTHOR/9_openWQ/2_initial_develop_synthetic_tests_GMD_paper/code/mizuRoute/case_studies/synthetic_tests/2_nrTrans_instS_PorMedia/mizuroute/mizuroute_out/v1.2_case1.h.2017-07-28-44100.nc"
output_html_path = "/Users/diogocosta/Downloads/results_animation.html"
#timestep = "all"
timestep = [0, 100]

h5_mplib.Map_h5_driver(
    shpfile_fullpath=shpfile_fullpath,
    openwq_results = openwq_results,
    hydromodel_out_fullpath = hydromodel_out_fullpath,
    output_html_path = output_html_path,
    timestep = timestep,
    model="mizuroute"
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