"""
Author: Mikhail Schee
Created: 2022-06-17

This script is set up to take the moving average of Arctic Ocean profile data
that has been formatted into netcdfs by the `make_netcdf` function, adding
the variables 'ma_iT', 'ma_CT', 'ma_SP', 'ma_SA', and 'ma_sigma' to those files.

"""

import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

import scipy.ndimage as ndimage
from scipy import interpolate

# Import the Thermodynamic Equation of Seawater 2010 (TEOS-10) from GSW
# For calculating density anomaly
import gsw

# The moving average window in dbar
c3 = 100

################################################################################
# Main execution
################################################################################
# Select the netcdfs to modify
ncs_to_modify = [
                #  'netcdfs/ITP_1.nc',
                 'netcdfs/ITP_2.nc',
                 'netcdfs/ITP_3.nc',
                #  'netcdfs/ITP_13.nc'
                 ]

# Loop through the netcdfs to modify
for my_nc in ncs_to_modify:
    print('Reading',my_nc)
    # Load in with xarray
    ds = xr.load_dataset(my_nc)

    gattrs_to_print = ['Last modified', 'Last modification', 'Moving average window']

    print('')
    # See the variables before
    for attr in gattrs_to_print:
        print('\t',attr+':',ds.attrs[attr])

    print('making changes')

    ## Get the moving average profiles
    # Convert a subset of the dataset to a dataframe
    df = ds[['press','iT','CT','PT','SP','SA']].squeeze().to_dataframe()
    # Get the original dimensions of the data to reshape the arrays later
    len0 = len(df.index.get_level_values(0).unique())
    len1 = len(df.index.get_level_values(1).unique())
    # Use the pandas `rolling` function to get the moving average
    #   center=True makes the first and last window/2 of the profiles are masked
    #   win_type='boxcar' uses a rectangular window shape
    #   on='press' means it will take `press` as the index column
    #   .mean() takes the average of the rolling
    df1 = df.rolling(window=c3, center=True, win_type='boxcar', on='press').mean()
    # Put the moving average profiles for temperature, salinity, and density into the dataset
    ds['ma_iT'].values = df1['iT'].values.reshape((len0, len1))
    ds['ma_CT'].values = df1['CT'].values.reshape((len0, len1))
    ds['ma_PT'].values = df1['PT'].values.reshape((len0, len1))
    ds['ma_SP'].values = df1['SP'].values.reshape((len0, len1))
    ds['ma_SA'].values = df1['SA'].values.reshape((len0, len1))
    ds['ma_sigma'].values= gsw.sigma1(ds['ma_SP'], ds['ma_CT'])

    # Update the global variables:
    ds.attrs['Last modified'] = str(datetime.now())
    ds.attrs['Last modification'] = 'Added moving averages with '+str(c3)+' dbar window'
    ds.attrs['Moving average window'] = str(c3)+' dbar'

    # Write out to netcdf
    print('Writing data to',my_nc)
    ds.to_netcdf(my_nc, 'w')

    # Load in with xarray
    ds2 = xr.load_dataset(my_nc)

    # See the variables after
    for attr in gattrs_to_print:
        print('\t',attr+':',ds2.attrs[attr])
