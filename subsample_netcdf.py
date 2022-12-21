"""
Author: Mikhail Schee
Created: 2022-04-11

This script is set up to subsample Arctic Ocean profile data that has been
formatted into netcdfs by the `make_netcdf` function, adding them to those files.

"""

import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

import scipy.ndimage as ndimage
from scipy import interpolate

################################################################################

# Select the subsampling scheme
# ss_scheme = 'AIDJEX PDF'
ss_scheme = 'Keep every'
# ss_scheme = 'Interpolate'
# Select integer x to Keep every x points
n_pts = 4
# Define the spacing (in dbar) for interpolation
interp_spacing = 2.0

################################################################################

def keep_every(vert, n_pts):
    """
    Sub-samples one profile, keeping only every <n_pts> points, skipping the rest
    Returns a mask which, when applied to the given vertical array, results in a
        subsampled vertical profile. Doesn't change any of the values

    vert    An array of vertical values to subsample (pressure or depth)
    n_pts   An integer so that only every n_pts point is not masked
    """
    # Get the length of the vertical array
    vert_len = len(vert)
    # Remove all nan values from vert assuming they will only be on the end
    vert = [x for x in vert if np.isnan(x)==False]
    # Check to make sure there is data to use
    if vert_len < 1:
        return [None]*vert_len
    else:
        # Make a blank array in which to put the mask
        ss_mask = np.array([None]*vert_len)
    #
    # Some really convenient syntax [::x] which returns a view of every x point
    #   Not a new array, changing it will change the original
    keep_these = ss_mask[::n_pts]
    # Unmask every n_pts point
    keep_these[:] = 1
    return ss_mask

################################################################################
# Main execution
################################################################################
# Select the netcdf to modify
ncs_to_modify = [
                 # 'netcdfs/ITP_1.nc',
                 'netcdfs/ITP_2.nc',
                 # 'netcdfs/ITP_3.nc'
                 ]

# Loop through the netcdfs to modify
for my_nc in ncs_to_modify:
    print('Reading',my_nc)
    # Load in with xarray
    ds = xr.load_dataset(my_nc)

    gattrs_to_print = ['Last modified', 'Last modification', 'Sub-sample scheme']

    print('')
    for attr in gattrs_to_print:
        print('\t',attr+':',ds.attrs[attr])

    print('making changes')

    ## Update the subsample mask
    # Find the original vertical variable (pressure or depth)
    vert_var = ds.attrs['Original vertical measure']

    # Loop over all the profiles
    if ss_scheme == 'AIDJEX_PDF':
        for i in range(len(ds[vert_var].values)):
            # Add subsample mask for this profile
            ds['ss_mask'][i] = AIDJEX_PDF(ds[vert_var].values[i])
    elif ss_scheme == 'Keep every':
        for i in range(len(ds[vert_var].values)):
            # Add subsample mask for this profile
            ds['ss_mask'][i] = keep_every(ds[vert_var].values[i], n_pts)
        ss_scheme += ' '+str(n_pts)+' point'
    # print('ss_mask:',ds['ss_mask'].values)

    # Update the global variables:
    ds.attrs['Last modified'] = str(datetime.now())
    ds.attrs['Last modification'] = 'Modified sub-sample scheme'
    ds.attrs['Sub-sample scheme'] = ss_scheme

    # Write out to netcdf
    print('Writing data to',my_nc)
    ds.to_netcdf(my_nc, 'w')

    # Load in with xarray
    ds2 = xr.load_dataset(my_nc)

    for attr in gattrs_to_print:
        print('\t',attr+':',ds2.attrs[attr])
