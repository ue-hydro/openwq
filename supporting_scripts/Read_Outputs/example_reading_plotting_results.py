#####
# Reading Openwq outputs and plot
#####

import sys
sys.path.insert(0, 'hdf5_support_lib')
import Read_h5_driver as h5_rlib
import Plot_h5_driver as h5_plib


###############
# CONFIGURATION
##############

res_dir = '/Users/diogocosta/Library/CloudStorage/OneDrive-impactblue-scientific.com/10_Univ_Evora/3_research/03_ACTIVE/1_MAIN_AUTHOR/9_openWQ/2_initial_develop_synthetic_tests_GMD_paper/code/mizuRoute/case_studies/synthetic_tests/2_nrTrans_instS_PorMedia/mizuroute/Output_OpenWQ'

# Output specs
OpenWQ_output_format = 'HDF5' # needs to be HDF5 (CSV disabled)
OpenWQ_output_debugmode = True

# What to print?
Chemical_species = ["SPECIES_A","SPECIES_B"]
Concentration_units = "MG/L"
Compartments = ['RIVER_NETWORK_REACHES']
Cells_spatial = 'all'

# Where to save plots?
plots_save_dir = '/Users/diogocosta/Downloads/'

###############
# CALLING THE NECESSARY SUPPORTING SCRIPS
##############

# Read requested data: Read_h5_driver
openwq_results = h5_rlib.Read_h5_driver(
            res_dir,
            OpenWQ_output_format,
            OpenWQ_output_debugmode,
            Compartments,
            Cells_spatial,
            Chemical_species,
            Concentration_units)

######
# Plot results
#####

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

# Example 2: Spatial plot - snapshot at specific date
figs4 = h5_plib.Plot_h5_driver(
        openwq_results,
        plot_type='spatial',
        time_slice=0,
        spatial_axis='x',
        fixed_coords={'y': 0, 'z': 0},
        output_dir=plots_save_dir,
        debug_mode=True
    )
