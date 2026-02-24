import netCDF4
import h5py
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import contextily as cx
import pyogrio
from datetime import datetime
import os
from shutil import copyfile

ncfile = '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/input/all_vars_E5evap.nc'

# Open original file (read-only)
ncfile_data = netCDF4.Dataset(ncfile, 'r')
allVars = list(ncfile_data.variables)

# Read var55
precp = np.array(ncfile_data["var55"][:])
print(f"Original var55 shape: {precp.shape}")
print(f"Original var55 range: [{np.min(precp):.3e}, {np.max(precp):.3e}]")

# Close the original file
ncfile_data.close()

# Create output filename
output_ncfile = ncfile.replace('.nc', '_var55modified.nc')

# Copy the entire file
print(f"\nCopying NetCDF file...")
copyfile(ncfile, output_ncfile)

# Open the copied file in read-write mode
print(f"Modifying var55 in: {output_ncfile}")
ncfile_modified = netCDF4.Dataset(output_ncfile, 'r+')

# MODIFY var55 - Choose one option:

# Option 1: Set to zero
#ncfile_modified.variables['var55'][:] = 0.0

# Option 2: Multiply by a factor (comment out option 1 and uncomment this)
ncfile_modified.variables['var55'][:] = 0.005

# Option 3: Set to custom values (comment out option 1 and uncomment this)
# custom_values = np.ones_like(precp) * 0.001
# ncfile_modified.variables['var55'][:] = custom_values

# Verify the change
modified_var55 = np.array(ncfile_modified.variables['var55'][:])
print(f"\nModified var55 range: [{np.min(modified_var55):.3e}, {np.max(modified_var55):.3e}]")

# Close the modified file
ncfile_modified.close()

print(f"\n✓ Modified NetCDF saved to: {output_ncfile}")