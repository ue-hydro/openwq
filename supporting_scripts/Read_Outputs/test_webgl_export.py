# Test script for WebGL 3D viewer export
# Uses test_case_2 data

import sys
sys.path.insert(0, 'hdf5_support_lib')

############################################
# General Mapping Info
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
# Read requested data
#######################
import Read_h5_driver as h5_rlib

print("=" * 70)
print("STEP 1: Reading OpenWQ HDF5 data...")
print("=" * 70)

openwq_results = h5_rlib.Read_h5_driver(
    openwq_info=openwq_info,
    output_format='HDF5',
    debugmode=False,
    cmp=['RIVER_NETWORK_REACHES'],
    space_elem='all',
    chemSpec=["NO3-N","NH4-N"],
    chemUnits="MG/L",
    noDataFlag=-9999,
    sediment_as_well=False)

print("\nAvailable species keys:")
for key in openwq_results.keys():
    print(f"  - {key}")

############################################
# Export WebGL viewer
############################################
print("\n" + "=" * 70)
print("STEP 2: Exporting WebGL viewer...")
print("=" * 70)

import WebGL_h5_driver as h5_wlib

result = h5_wlib.WebGL_h5_driver(
    what2map='openwq',
    hostmodel='mizuroute',
    shpfile_info=shpfile_info,
    openwq_results=openwq_results,
    chemSpec=["NO3-N","NH4-N"],
    sediment_as_well=False,
    hydromodel_info=hydromodel_info,
    hydromodel_var2print='DWroutedRunoff',
    output_dir='/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/bin/openwq_out/openwq_webgl_viewer',
    timeframes=20,
    dem_path=None,
    n_particles=65536,
    river_width_cells=3,
    grid_resolution=None
)

print("\nResult:", result)
