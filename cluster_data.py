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

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions in source code must retain the accompanying copyright notice, this list of conditions, and the following disclaimer.
    2. Redistributions in binary form must reproduce the accompanying copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3. Names of the copyright holders must not be used to endorse or promote products derived from this software without prior written permission from the copyright holders.
    4. If any files are modified, you must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

gattrs_to_print =  ['Last clustered',
                    'Clustering x-axis',
                    'Clustering y-axis',
                    'Clustering m_pts',
                    'Clustering filters',
                    'Clustering DBCV']

dfs0 = ahf.Data_Filters()

## Reproducing figures from Timmermans et al. 2008
ITP2_clstr_dict = {'netcdf_to_load':'netcdfs/ITP_2.nc',
                   'sources_dict':{'ITP_2':'all'},
                   'data_filters':dfs0,
                   'SP_range':[34.05,34.75],
                   'cl_x_var':'SP',
                   'cl_y_var':'la_CT',
                   'm_pts':170
                   }

## Reproducing figures from Lu et al. 2022
ITP3_clstr_dict = {'netcdf_to_load':'netcdfs/ITP_3.nc',
                   'sources_dict':{'ITP_3':'all'},
                   'data_filters':dfs0,
                   'SP_range':[34.21,34.82],
                   'cl_x_var':'SP',
                   'cl_y_var':'la_CT',
                   'm_pts':580
                   }

for clstr_dict in [ITP2_clstr_dict, ITP3_clstr_dict]:
    # Find the netcdf to use
    my_nc = clstr_dict['netcdf_to_load']
    print('Reading',my_nc)
    # Load in with xarray
    ds = xr.load_dataset(my_nc)

    # Convert a subset of the dataset to a dataframe
    og_df = ds[['cluster', 'clst_prob']].squeeze().to_dataframe()
    # Get the original dimensions of the data to reshape the arrays later
    len0 = len(og_df.index.get_level_values(0).unique())
    len1 = len(og_df.index.get_level_values(1).unique())
    # print('Dataframe from netcdf, selected columns:')
    # print(og_df)
    # print('')
    # print('Dataframe from netcdf, non-NaN rows:')
    # print(og_df[~og_df.isnull().any(axis=1)])
    # print('')
    # See the variables before
    for attr in gattrs_to_print:
        print('\t',attr+':',ds.attrs[attr])

    # Create data set object
    ds_ITP2_all = ahf.Data_Set(clstr_dict['sources_dict'], clstr_dict['data_filters'])
    # Create profile filter object
    pfs_T2008  = ahf.Profile_Filters(SP_range=clstr_dict['SP_range'])
    # Create plot parameters object
    pp_T2008_clstr = ahf.Plot_Parameters(x_vars=[clstr_dict['cl_x_var']], y_vars=[clstr_dict['cl_y_var']], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':clstr_dict['cl_x_var'], 'cl_y_var':clstr_dict['cl_y_var'], 'm_pts':clstr_dict['m_pts']}, legend=False)
    # Create analysis group
    group_T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_clstr)
    # Make a figure to run clustering algorithm and check results
    ahf.make_figure([group_T2008_clstr])

    print('making changes')
    # Pull the now modified dataframe from the analysis group
    new_df = group_T2008_clstr.data_frames[0][['cluster','clst_prob']]
    # print('Dataframe from plotting, selected columns:')
    # print(new_df)
    # print('')
    # print('Dataframe from plotting, non-NaN rows:')
    # print(new_df[~new_df.isnull().any(axis=1)])
    # Update the original dataframe from the netcdf with the filtered dataframe from plotting
    og_df.update(new_df)
    # print('Dataframe from netcdf, updated:')
    # print(og_df)
    # print('')
    # print('Dataframe from netcdf, non-NaN rows:')
    # print(og_df[~og_df.isnull().any(axis=1)])
    # Put the clustering variables back into the dataset
    ds['cluster'].values   = og_df['cluster'].values.reshape((len0, len1))
    ds['clst_prob'].values = og_df['clst_prob'].values.reshape((len0, len1))
    # Update the global variables:
    ds.attrs['Last modified'] = str(datetime.now())
    ds.attrs['Last modification'] = 'Updated clustering'
    ds.attrs['Last clustered'] = str(datetime.now())
    ds.attrs['Clustering x-axis'] = clstr_dict['cl_x_var']
    ds.attrs['Clustering y-axis'] = clstr_dict['cl_y_var']
    ds.attrs['Clustering m_pts'] = clstr_dict['m_pts']
    ds.attrs['Clustering filters'] = 'SP_range: '+str(clstr_dict['SP_range'])
    ds.attrs['Clustering DBCV'] = group_T2008_clstr.data_set.arr_of_ds[0].attrs['Clustering DBCV']
    # Write out to netcdf
    print('Writing data to',my_nc)
    ds.to_netcdf(my_nc, 'w')
    # Load in with xarray
    ds2 = xr.load_dataset(my_nc)
    # See the variables after
    for attr in gattrs_to_print:
        print('\t',attr+':',ds2.attrs[attr])
