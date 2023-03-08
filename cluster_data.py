"""
Author: Mikhail Schee
Created: 2023-03-06

This script is set up to run the HDBSCAN clustering algorithm on a set of ITP 
data given certain inputs. It will fill in values in the 'cluster' and 
'clst_prob' columns, as well as for the global attributes:
'Last clustered'
'Clustering x-axis'
'Clustering y-axis'
'Clustering m_pts'
'Clustering ranges'
'Clustering DBCV'

"""

import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

# Import the Thermodynamic Equation of Seawater 2010 (TEOS-10) from GSW
# For converting from depth to pressure
import gsw

# For custom analysis functions
import analysis_helper_functions as ahf

################################################################################
# Main execution
################################################################################
# Select the netcdf to modify
my_nc = 'netcdfs/ITP_2.nc'

print('Reading',my_nc)
# Load in with xarray
ds = xr.load_dataset(my_nc)

gattrs_to_print =  ['Last clustered',
                    'Clustering x-axis',
                    'Clustering y-axis',
                    'Clustering m_pts',
                    'Clustering filters',
                    'Clustering DBCV']

print('')
# See the variables before
for attr in gattrs_to_print:
    print('\t',attr+':',ds.attrs[attr])

print('making changes')

## Reproducing figures from Timmermans et al. 2008
ITP2_all  = {'ITP_2':'all'}
ITP2_p_range = [185,300]
T2008_m_pts = 100
dfs0 = ahf.Data_Filters()
ds_ITP2_all = ahf.Data_Set(ITP2_all, dfs0)
pfs_T2008  = ahf.Profile_Filters(p_range=ITP2_p_range)
## The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_T2008_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=False)
group_T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_clstr)
ahf.make_figure([group_T2008_clstr])

exit(0)

## Get the moving average profiles
# Convert a subset of the dataset to a dataframe
vars_to_keep = ['entry']
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
