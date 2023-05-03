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

# print('before')
# print(ds['cluster'].loc[dict(Time='2004-08-20 00:00:01', Vertical=1000)].values)

## Put the clustering results back into the netcdf file
# Convert a subset of the dataset to a dataframe
og_df = ds[['cluster', 'clst_prob']].squeeze().to_dataframe()
# Get the original dimensions of the data to reshape the arrays later
len0 = len(og_df.index.get_level_values(0).unique())
len1 = len(og_df.index.get_level_values(1).unique())

print('Dataframe from netcdf, selected columns:')
print(og_df)
print('')

print('Dataframe from netcdf, non-NaN rows:')
print(og_df[~og_df.isnull().any(axis=1)])

# og_df.at[('2004-08-20 00:00:01', 1000), 'cluster'] = 3

# ds['cluster'].values = og_df['cluster'].values.reshape((len0, len1))

# print('after')
# print(ds['cluster'].loc[dict(Time='2004-08-20 00:00:01', Vertical=1000)].values)

# exit(0)

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

dfs0 = ahf.Data_Filters()

## Reproducing figures from Timmermans et al. 2008
ITP2_clst_dict = {'sources_dict':{'ITP_2':'all'},
                  'data_filters':dfs0,
                  'SP_range':[34.05,34.75],
                  'cl_x_var':'SP',
                  'cl_y_var':'la_CT',
                  'm_pts':170
                 }


ds_ITP2_all = ahf.Data_Set(ITP2_clst_dict['sources_dict'], ITP2_clst_dict['data_filters'])
pfs_T2008  = ahf.Profile_Filters(SP_range=ITP2_clst_dict['SP_range'])
## The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_T2008_clstr = ahf.Plot_Parameters(x_vars=[ITP2_clst_dict['cl_x_var']], y_vars=[ITP2_clst_dict['cl_y_var']], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':ITP2_clst_dict['cl_x_var'], 'cl_y_var':ITP2_clst_dict['cl_y_var'], 'm_pts':ITP2_clst_dict['m_pts']}, legend=False)

group_T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_clstr)

# print('Before:')
# print(group_T2008_clstr.data_frames[0])
# print('DBCV:',group_T2008_clstr.data_set.arr_of_ds[0].attrs['Clustering DBCV'])

ahf.make_figure([group_T2008_clstr])

# print('After:')
# print(group_T2008_clstr.data_frames[0])
# print('DBCV:',group_T2008_clstr.data_set.arr_of_ds[0].attrs['Clustering DBCV'])

new_df = group_T2008_clstr.data_frames[0][['cluster','clst_prob']]
print('Dataframe from plotting, selected columns:')
print(new_df)
print('')

print('Dataframe from plotting, non-NaN rows:')
print(new_df[~new_df.isnull().any(axis=1)])

# Update the original dataframe from the netcdf with the filtered dataframe from plotting
og_df.update(new_df)

print('Dataframe from netcdf, updated:')
print(og_df)
print('')

print('Dataframe from netcdf, non-NaN rows:')
print(og_df[~og_df.isnull().any(axis=1)])

# exit(0)
# Put the clustering variables back into the dataset
ds['cluster'].values   = og_df['cluster'].values.reshape((len0, len1))
ds['clst_prob'].values = og_df['clst_prob'].values.reshape((len0, len1))

# Update the global variables:
ds.attrs['Last modified'] = str(datetime.now())
ds.attrs['Last modification'] = 'Updated clustering'
ds.attrs['Last clustered'] = str(datetime.now())
ds.attrs['Clustering x-axis'] = ITP2_clst_dict['cl_x_var']
ds.attrs['Clustering y-axis'] = ITP2_clst_dict['cl_y_var']
ds.attrs['Clustering m_pts'] = ITP2_clst_dict['m_pts']
ds.attrs['Clustering filters'] = 'SP_range: '+str(ITP2_clst_dict['SP_range'])
ds.attrs['Clustering DBCV'] = group_T2008_clstr.data_set.arr_of_ds[0].attrs['Clustering DBCV']

# Write out to netcdf
print('Writing data to',my_nc)
ds.to_netcdf(my_nc, 'w')

# Load in with xarray
ds2 = xr.load_dataset(my_nc)

# See the variables after
for attr in gattrs_to_print:
    print('\t',attr+':',ds2.attrs[attr])
