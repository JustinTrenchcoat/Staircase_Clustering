"""
Author: Mikhail Schee
Created: 2022-03-31

This script contains helper functions which filter data and create figures

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions in source code must retain the accompanying copyright notice, this list of conditions, and the following disclaimer.
    2. Redistributions in binary form must reproduce the accompanying copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3. Names of the copyright holders must not be used to endorse or promote products derived from this software without prior written permission from the copyright holders.
    4. If any files are modified, you must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# For making insets in plots
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# For getting different marker styles
from matplotlib.markers import MarkerStyle as mplms
# For storing figure objects in files (can use `pickle` instead if need be)
import dill as pl
# For formatting data into dataframes
import pandas as pd
# For matching regular expressions
import re
# For formatting date objects
from datetime import datetime
# For reading netcdf files
import xarray as xr
# Import the Thermodynamic Equation of Seawater 2010 (TEOS-10) from GSW
# For finding alpha and beta values
import gsw
# For adding subplot labels a, b, c, ...
import string
# For taking moving averages
import scipy.ndimage as ndimage
# For getting zscores to find outliers
from scipy import stats
# For making clusters
import hdbscan
# For calculating the distance between pairs of (latitude, longitude)
from geopy.distance import geodesic
# For calculating Orthogonal Distance Regression for Total Least Squares
from orthoregress import orthoregress

"""
To install Cartopy and its dependencies, follow:
https://scitools.org.uk/cartopy/docs/latest/installing.html#installing

Relevent command:
$ conda install -c conda-forge cartopy
"""
try:
    import cartopy.crs as ccrs
    import cartopy.feature
except:
    print('Failed to import cartopy\nIf you do not need to plot maps, comment out these import statements')
    exit(0)

science_data_file_path = '/Users/Grey/Documents/Research/Science_Data/'

# This list gets filled in with the names of all the available columns for the
#   given data during the `apply_data_filters()` function
available_variables_list = []

################################################################################
# Declare variables for plotting
################################################################################
# Be sure to set dark_mode similiarly in the other scripts
dark_mode = False

# Colorblind-friendly palette by Krzywinski et al. (http://mkweb.bcgsc.ca/biovis2012/)
#   See the link below for a helpful color wheel:
#   https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
jackson_clr = np.array(["#000000",  #  0 black
                        "#004949",  #  1 dark olive
                        "#009292",  #  2 teal
                        "#ff6db6",  #  3 hot pink
                        "#ffb6db",  #  4 light pink
                        "#490092",  #  5 dark purple
                        "#006ddb",  #  6 royal blue
                        "#b66dff",  #  7 violet
                        "#6db6ff",  #  8 sky blue
                        "#b6dbff",  #  9 pale blue
                        "#920000",  # 10 dark red
                        "#924900",  # 11 brown
                        "#db6d00",  # 12 dark orange
                        "#24ff24",  # 13 neon green
                        "#ffff6d"]) # 14 yellow

# Enable dark mode plotting
if dark_mode:
    plt.style.use('dark_background')
    std_clr = 'w'
    alt_std_clr = 'yellow'
    noise_clr = 'w'
    cnt_clr = '#db6d00'
    clr_ocean = 'k'
    clr_land  = 'grey'
    clr_lines = 'w'
    clstr_clrs = jackson_clr[[2,14,4,6,7,9,13]]
    bathy_clrs = ['#000040','k']
else:
    std_clr = 'k'
    alt_std_clr = 'olive'
    noise_clr = 'k'
    cnt_clr = "#ffff6d"
    clr_ocean = 'w'
    clr_land  = 'grey'
    clr_lines = 'k'
    clstr_clrs = jackson_clr[[12,1,10,6,5,2,13]]
    bathy_clrs = ['w','#b6dbff']

# Define bathymetry colors
n_bathy = 5
cm0 = mpl.colors.LinearSegmentedColormap.from_list("Custom", bathy_clrs, N=n_bathy)
cm1 = mpl.colors.LinearSegmentedColormap.from_list("Custom", bathy_clrs[::-1], N=n_bathy)
# Initialize colormap to get ._lut attribute
cm0._init()
cm1._init()
# There are 3 extra colors at the end for some weird reason, so eliminate them
rgbas0 = cm0._lut[0:-3]
# Convert to a list of hex values in strings
bathy_clrs = [mpl.colors.rgb2hex(x) for x in rgbas0]
# Make an object for color bars
bathy_smap = plt.cm.ScalarMappable(cmap=cm1)
bathy_smap.set_array([])

# Set some plotting styles
font_size_plt = 12
font_size_labels = 14
font_size_ticks = 12
font_size_lgnd = 10
mpl.rcParams['font.size'] = font_size_plt
mpl.rcParams['axes.labelsize'] = font_size_labels
mpl.rcParams['xtick.labelsize'] = font_size_ticks
mpl.rcParams['ytick.labelsize'] = font_size_ticks
mpl.rcParams['legend.fontsize'] = font_size_lgnd
mrk_size      = 0.5
mrk_alpha     = 0.3
noise_alpha   = 0.2
grid_alpha    = 0.3
pf_alpha      = 0.4
map_alpha     = 0.7
lgnd_mrk_size = 60
map_mrk_size  = 25
big_map_mrkr  = 80
sml_map_mrkr  = 5
cent_mrk_size = 50
pf_mrk_size   = 15
std_marker = '.'
map_marker = '.' #'x'
map_ln_wid = 0.5
l_cap_size = 3.0

# The magnitude limits for axis ticks before they use scientific notation
sci_lims = (-2,3)

#   Get list of standard colors
distinct_clrs = clstr_clrs
#   Make list of marker and fill styles
unit_star_3 = mpl.path.Path.unit_regular_star(3)
unit_star_4 = mpl.path.Path.unit_regular_star(4)
mpl_mrks = [mplms('o',fillstyle='left'), mplms('o',fillstyle='right'), 'x', unit_star_4, '*', unit_star_3, '1','+']
# Define array of linestyles to cycle through
l_styles = ['-', '--']#, '-.', ':']

# A list of variables for which the y-axis should be inverted so the surface is up
y_invert_vars = ['press', 'pca_press', 'ca_press', 'cmm_mid', 'pca_depth', 'ca_depth', 'sigma', 'ma_sigma', 'pca_sigma', 'ca_sigma', 'pca_iT', 'ca_iT', 'pca_CT', 'ca_CT', 'pca_PT', 'ca_PT', 'pca_SP', 'ca_SP', 'pca_SA', 'ca_SA']
# A list of the per profile variables
pf_vars = ['entry', 'prof_no', 'BL_yn', 'dt_start', 'dt_end', 'lon', 'lat', 'region', 'up_cast', 'R_rho']
# A list of the variables on the `Vertical` dimension
vertical_vars = ['press', 'depth', 'iT', 'CT', 'PT', 'SP', 'v1_SP', 'SA', 'sigma', 'alpha', 'beta', 'aiT', 'aCT', 'aPT', 'BSP', 'BSA', 'ss_mask', 'ma_iT', 'ma_CT', 'ma_PT', 'ma_SP', 'ma_SA', 'ma_sigma', 'la_iT', 'la_CT', 'v2_CT', 'la_PT', 'la_SP', 'la_SA', 'la_sigma']
# Make lists of clustering variables
pca_prefix = 'pca_' # profile cluster average
pca_vars   = [ f'{pca_prefix}{var}' for var in vertical_vars]
pcs_prefix = 'pcs_' # profile cluster span
pcs_vars   = [ f'{pcs_prefix}{var}' for var in vertical_vars]
cmc_prefix = 'cmc_' # cluster mean-centered
cmc_vars   = [ f'{cmc_prefix}{var}' for var in vertical_vars]
ca_prefix  = 'ca_'  # cluster average
ca_vars    = [ f'{ca_prefix}{var}'  for var in vertical_vars]
cs_prefix  = 'cs_'  # cluster span
cs_vars    = [ f'{cs_prefix}{var}'  for var in vertical_vars]
cmm_prefix = 'cmm_' # cluster min/max
cmm_vars   = [ f'{cmm_prefix}{var}'  for var in vertical_vars]
nir_prefix = 'nir_' # normalized inter-cluster range, with max-min instead of stdv
nir_vars   = [ f'{nir_prefix}{var}'  for var in vertical_vars]
# Make a complete list of cluster-related variables
clstr_vars = ['cluster', 'cRL', 'cRl'] + pca_vars + pcs_vars + cmc_vars + ca_vars + cs_vars + cmm_vars + nir_vars
# For parameter sweeps of clustering
#   Independent variables
clstr_ps_ind_vars = ['m_pts', 'n_pfs', 'ell_size']
#   Dependent variables
clstr_ps_dep_vars = ['DBCV', 'n_clusters']
clstr_ps_vars = clstr_ps_ind_vars + clstr_ps_dep_vars

################################################################################
# Declare classes for custom objects
################################################################################

class Data_Filters:
    """
    A set of filters to select certain profiles
    Note: Only filters to a selection of profiles. Does not filter individual
        profiles in any manner

    keep_black_list     True/False whether to keep profiles marked on black list
    cast_direction      'up','down',or 'all' to keep certain directions of casts
    geo_extent          None to keep all profiles regardless of geographic region
                        or 'CS' for Chukchi Sea, 'SBS' for Southern Beaufort Sea,
                        'CB' for Canada Basin, 'MB' for Makarov Basin, 'EB' for
                        Eurasian Basin, or 'BS' for Barents Sea
    date_range          ['start_date','end_date'] where the dates are strings in
                        the format 'YYYY/MM/DD' or None to keep all profiles
    """
    def __init__(self, keep_black_list=False, cast_direction='up', geo_extent=None, date_range=None):
        self.keep_black_list = keep_black_list
        self.cast_direction = cast_direction
        self.geo_extent = geo_extent
        self.date_range = date_range

################################################################################

class Data_Set:
    """
    Loads data from netcdfs and returns a list of xarrays with filters applied
    Note: Only filters to a selection of profiles. Does not filter individual
        profiles in any manner

    sources_dict    A dictionary with the following format:
                    {'ITP_1':[13,22,32],'ITP_2':'all'}
                    where the keys are the netcdf filenames without the extension
                    and the values are lists of profiles to include or 'all'
    data_filters    A custom Data_Filters object that contains the filters to apply
    """
    def __init__(self, sources_dict, data_filters):
        # Load just the relevant profiles into the xarrays
        self.sources_dict = sources_dict
        self.data_filters = data_filters
        xarrs, self.var_attr_dicts = list_xarrays(sources_dict)
        self.arr_of_ds = apply_data_filters(xarrs, data_filters)

################################################################################

class Profile_Filters:
    """
    A set of filters to apply to each profile. If None is passed for any filter,
    that filter will not be applied. All these are ignored when plotting 'by_pf'

    p_range             [p_min, p_max] where the values are floats in dbar
    d_range             [d_min, d_max] where the values are floats in m
    *T_range            [T_min, T_max] where the values are floats in degrees C
    S*_range            [S_min, S_max] where the values are floats in g/kg
    subsample           True/False whether to apply the subsample mask to the profiles
    regrid_TS           [1st_var_str, Delta_1st_var, 2nd_var_str, Delta_2nd_var], a pair of 
                            [var, Delta_var] where you specify the variable then the value
                            of the spacing to regrid that value to
    m_avg_win           The value in dbar of the moving average window to take for ma_ variables
                            This is multiplied by 4 to get the number of rows to average
    """
    def __init__(self, p_range=None, d_range=None, iT_range=None, CT_range=None, PT_range=None, SP_range=None, SA_range=None, subsample=False, regrid_TS=None, m_avg_win=None):
        self.p_range = p_range
        self.d_range = d_range
        self.iT_range = iT_range
        self.CT_range = CT_range
        self.PT_range = PT_range
        self.SP_range = SP_range
        self.SA_range = SA_range
        self.subsample = subsample
        self.regrid_TS = regrid_TS
        self.m_avg_win = m_avg_win

################################################################################

class Plot_Parameters:
    """
    A set of parameters to determine what kind of plot to make

    plot_type       A string of the kind of plot to make from these options:
                        'xy', 'map', 'profiles'
    plot_scale      A string specifying either 'by_pf' or 'by_vert'
                    which will determine whether to add one data point per
                    profile or per vertical observation, respectively
    x_vars          A list of strings of the variables to plot on the x axis
    y_vars          A list of strings of the variables to plot on the y axis
                    Note: 'xy' and 'profiles' expect both x_vars and y_vars and
                      'map' ignores this variable as it always plots lon vs. lat
                    Default plot is 'SP' vs. 'iT'
                    Accepted 'by_pf' vars: 'entry', 'prof_no', 'BL_yn', 'dt_start',
                      'dt_end', 'lon', 'lat', 'region', 'up_cast', 'R_rho', 'hist'
                    Accepted 'by_vert' vars: all 'by_pf' vars plus 'press',
                      'depth', 'iT', 'CT', 'PT', 'SP', 'SA', 'sigma', 'alpha', 'beta',
                      'aiT', 'aCT', 'aPT', 'BSP', 'BSA', 'ma_iT', 'ma_CT', 'ma_PT', 
                      'ma_SP', 'ma_SA', 'ma_sigma', 'ss_mask', 'la_iT', 'la_CT', 
                      'la_PT', 'la_SA','la_sigma'
    clr_map         A string to determine what color map to use in the plot
                    'xy' can use 'clr_all_same', 'clr_by_source',
                      'clr_by_instrmt', 'density_hist', 'cluster', or any
                      regular variable
                    'map' can use 'clr_all_same', 'clr_by_source',
                      'clr_by_instrmt', or any 'by_pf' regular variable
                    'profiles' can use 'clr_all_same' or any regular variable
    first_dfs       A list of booleans of whether to take the first differences
                      of the plot_vars before analysis. Defaults to all False
                    'xy' and 'profiles' expect a list of 2, ex: [True, False]
                    'map' and ignores this parameter
    finit_dfs       A list of booleans of whether to take the finite differences
                      of the plot_vars before analysis. Defaults to all False
                    'xy' and 'profiles' expect a list of 2, ex: [True, False]
                    'map' and ignores this parameter
    legend          True/False whether to add a legend on this plot. Default:True
    add_grid        True/False whether to add a grid on this plot. Default:True
    ax_lims         An optional dictionary of the limits on the axes for the final plot
                      Ex: {'x_lims':[x_min,x_max], 'y_lims':[y_min,y_max]}
    extra_args      A general use argument for passing extra info for the plot
                    'map' expects a dictionary following this format:
                        {'map_extent':'Canada_Basin'}
                    'density_hist' expects a dictionary following this format:
                        {'clr_min':0, 'clr_max':20, 'clr_ext':'max', 'xy_bins':250}
                        All 4 of those arguments must be provided for it to work
                        If no extra_args is given, those default values are used
                    If 'hist' is one of the plot variables, you can add the
                        number of histogram bins to use: {'n_h_bins':25}, if not
                        given, it will use that default value
                    'profiles' plots will accept optional arguments such as:
                        {'pfs_to_plot':[183,185,187]}
                        Note: the code narrows to these profiles just before
                        plotting, so use this argument instead of narrowing earlier
                        in the Data_Set object if you want to plot by 'cluster'
                    If plotting anything that has to do clustering, need these:
                        {'cl_x_var':var0, 'cl_y_var':var1} where the two variables
                        can be any that work for an 'xy' plot
                        Optional 'b_a_w_plt':True/False for box and whisker plots,
                        'm_pts':90 / 'min_samp':90 to specify m_pts or min pts per
                        cluster, 'plot_slopes' to plot lines showing the least
                        squares slope of each cluster
                    If doing a parameter sweep of clustering, expects the following:
                        {'cl_x_var':var0, 'cl_y_var':var1, 'cl_ps_tuple':[100,410,50]}
                        where var0 and var1 are as specified above, 'cl_ps_tuple' is
                        the [start,stop,step] for the xvar can be 'm_pts', 'min_samps'
                        'n_pfs' or 'ell_size', and yvar(s) must be 'DBCV' or 'n_clusters'
                        Optional: {'z_var':'m_pts', 'z_list':[90,120,240]} where
                        z_var can be any variable that var0 can be
                    To plot isopycnal contour lines add {'isopycnals':X} where X is the 
                        value in dbar to which the isopycnals are referenced or True 
                        which will set the reference to the median pressure
                    To manually place the inline labels for isopycnals, add this:
                        {'place_isos':'manual'}, otherwise they will be placed 
                        automatically
    """
    def __init__(self, plot_type='xy', plot_scale='by_vert', x_vars=['SP'], y_vars=['CT'], clr_map='clr_all_same', first_dfs=[False, False], finit_dfs=[False, False], legend=True, add_grid=True, ax_lims=None, extra_args=None):
        # Add all the input parameters to the object
        self.plot_type = plot_type
        self.plot_scale = plot_scale
        self.x_vars = x_vars
        self.y_vars = y_vars
        self.xlabels = [None,None]
        self.ylabels = [None,None]
        self.clr_map = clr_map
        self.clabel = None
        self.first_dfs = first_dfs
        self.finit_dfs = finit_dfs
        self.legend = legend
        self.add_grid = add_grid
        self.ax_lims = ax_lims
        # If trying to plot a map, make sure x and y vars are None
        if plot_type == 'map':
            x_vars = None
            y_vars = None
        # If trying to plot a map, make sure `map_extent` exists
        if plot_type == 'map' and isinstance(extra_args, type(None)):
            self.extra_args = {'map_extent':'Canada_Basin'}
        else:
            self.extra_args = extra_args

################################################################################

class Analysis_Group:
    """
    Takes in a list of xarray datasets and applies filters to individual profiles
    Returns a list of pandas dataframes, one for each source

    data_set            A Data_Set object
    profile_filters     A custom Profile_Filters object with filters to apply to
                        all individual profiles in the xarrays
    plt_params          A custom Plot_Parameters object
    plot_title          A string to use as the title for this subplot
    """
    def __init__(self, data_set, profile_filters, plt_params, plot_title=None):
        self.data_set = data_set
        self.vars_available = list(data_set.arr_of_ds[0].keys())
        self.profile_filters = profile_filters
        self.plt_params = get_axis_labels(plt_params, data_set.var_attr_dicts)
        self.plot_title = plot_title
        self.vars_to_keep = find_vars_to_keep(plt_params, profile_filters, self.vars_available)
        # Load just the relevant profiles into the xarrays
        self.data_frames = apply_profile_filters(data_set.arr_of_ds, self.vars_to_keep, profile_filters, plt_params)

################################################################################
# Define class functions #######################################################
################################################################################

def list_xarrays(sources_dict):
    """
    Returns a list of xarrays, one for each data source as specified by the
    input dictionary

    sources_dict    A dictionary with the following format:
                    {'ITP_1':[13,22,32],'ITP_2':'all'}
                    where the keys are the netcdf filenames without the extension
                    and the values are lists of profiles to include or 'all'
    """
    # Make an empty list
    xarrays = []
    var_attr_dicts = []
    for source in sources_dict.keys():
        ds = xr.load_dataset('netcdfs/'+source+'.nc')
        # Build the dictionary of netcdf attributes, variables, units, etc.
        var_attrs = {}
        for this_var in list(ds.keys()):
            var_attrs[ds[this_var].name] = ds[this_var].attrs
        var_attr_dicts.append(var_attrs)
        # Convert the datetime variables from dtype `object` to `datetime64`
        ds['Time'] = pd.DatetimeIndex(ds['Time'].values)
        # Check whether to only take certain profiles
        pf_list = sources_dict[source]
        if isinstance(pf_list, list):
            # Start a blank list for all the profile-specific xarrays
            temp_list = []
            for pf in pf_list:
                temp_list.append(ds.where(ds.prof_no==pf, drop=True).squeeze())
            # Concatonate all those xarrays together
            ds1 = xr.concat(temp_list, 'Time')
            xarrays.append(ds1)
        # Take all profiles
        elif pf_list=='all':
            xarrays.append(ds)
        else:
            print(pf_list,'is not a valid list of profiles')
            exit(0)
        #
    # print(var_attr_dicts)
    return xarrays, var_attr_dicts

################################################################################

def apply_data_filters(xarrays, data_filters):
    """
    Returns a list of the same xarrays as provided, but with the filters applied

    xarrays         A list of xarray dataset objects
    data_filters    A custom Data_Filters object that contains the filters to apply
    """
    # Make an empty list
    output_arrs = []
    for ds in xarrays:
        ## Filters on a per-profile basis
        ##      I turned off squeezing because it drops the `Time` dimension if
        ##      you only pass in one profile per dataset
        #   Filter based on the black list
        if data_filters.keep_black_list == False:
            ds = ds.where(ds.BL_yn==False, drop=True)#.squeeze()
        #   Filter based on the cast direction
        if data_filters.cast_direction == 'up':
            ds = ds.where(ds.up_cast==True, drop=True)#.squeeze()
        elif data_filters.cast_direction == 'down':
            ds = ds.where(ds.up_cast==False, drop=True)#.squeeze()
        #   Filter based on the geographical region
        if data_filters.geo_extent == 'CB':
            ds = ds.where(ds.region=='CB', drop=True)#.squeeze()
        #   Filter based on the date range
        if not isinstance(data_filters.date_range, type(None)):
            # Allow for date ranges that do or don't specify the time
            try:
                start_date_range = datetime.strptime(data_filters.date_range[0], r'%Y/%m/%d %H:%M:%S')
            except:
                start_date_range = datetime.strptime(data_filters.date_range[0], r'%Y/%m/%d')
            try:
                end_date_range   = datetime.strptime(data_filters.date_range[1], r'%Y/%m/%d %H:%M:%S')
            except:
                end_date_range   = datetime.strptime(data_filters.date_range[1], r'%Y/%m/%d')
            ds = ds.sel(Time=slice(start_date_range, end_date_range))
        #
        output_arrs.append(ds)
    return output_arrs

################################################################################

def find_vars_to_keep(pp, profile_filters, vars_available):
    """
    Returns a list the variables to keep for an analysis group based on the plot
    parameters and profile filters given

    pp                  A custom Plot_Parameters object
    profile_filters     A custom Profile_Filters object with filters to apply to
                        all individual profiles in the xarrays
    vars_available      A list of variables available for this data
    """
    # If plot type is 'summary', then keep all the variables
    if pp.plot_type == 'summary':
        # Check the plot scale
        #   This will include all the measurements within each profile
        if plt_params.plot_scale == 'by_vert':
            vars_to_keep = ['entry', 'prof_no', 'dt_start', 'dt_end', 'lon', 'lat', 'R_rho', 'press', 'depth', 'iT', 'CT', 'PT', 'SP', 'SA', 'sigma', 'alpha', 'beta', 'ss_mask', 'ma_iT', 'ma_CT', 'ma_PT', 'ma_SP', 'ma_SA', 'ma_sigma']
        #   This will only include variables with 1 value per profile
        elif plt_params.plot_scale == 'by_pf':
            vars_to_keep = ['entry', 'prof_no', 'dt_start', 'dt_end', 'lon', 'lat', 'R_rho']
        # print('vars_to_keep:')
        # print(vars_to_keep)
        return vars_to_keep
    else:
        # Always include the entry and profile numbers
        vars_to_keep = ['entry', 'prof_no']
        # Add vars based on plot type
        if pp.plot_type == 'map':
            vars_to_keep.append('lon')
            vars_to_keep.append('lat')
        # Add all the plotting variables
        re_run_clstr = False
        plot_vars = pp.x_vars+pp.y_vars+[pp.clr_map]
        # Check the extra arguments
        if not isinstance(pp.extra_args, type(None)):
            for key in pp.extra_args.keys():
                if key in ['cl_x_var', 'cl_y_var']:
                    plot_vars.append(pp.extra_args[key])
                    re_run_clstr = True
                if key == 'z_var':
                    if pp.extra_args[key] == 'ell_size':
                        # Make sure the parameter sweeps run correctly 
                        vars_to_keep.append('press')
                # If adding isopycnals, make sure to keep press, SA, and iT to use with the 
                #   function gsw.pot_rho_t_exact(SA,t,p,p_ref)
                if key == 'isopycnals':
                    if not isinstance(pp.extra_args[key], type(None)):
                        if not 'press' in vars_to_keep:
                            vars_to_keep.append('press')
                        if not 'SA' in vars_to_keep:
                            vars_to_keep.append('SA')
                        if not 'iT' in vars_to_keep:
                            vars_to_keep.append('iT')
                        #
                    #
                #
            #
        # print('plot_vars:',plot_vars)
        # print('vars_available:',vars_available)
        for var in plot_vars:
            if var in vars_available:
                vars_to_keep.append(var)
            if var == 'cluster':
                vars_to_keep.append('clst_prob')
            if var == 'aiT':
                vars_to_keep.append(var)
                vars_to_keep.append('alpha_iT')
                vars_to_keep.append('iT')
            elif var == 'aCT':
                vars_to_keep.append(var)
                vars_to_keep.append('alpha')
                vars_to_keep.append('CT')
            elif var == 'aPT':
                vars_to_keep.append(var)
                vars_to_keep.append('alpha_PT')
                vars_to_keep.append('PT')
            elif var == 'BSP':
                vars_to_keep.append(var)
                vars_to_keep.append('beta')
                vars_to_keep.append('SP')
            elif var == 'BSt':
                vars_to_keep.append(var)
                vars_to_keep.append('beta_PT')
                vars_to_keep.append('SP')
            elif var == 'BSA':
                vars_to_keep.append(var)
                vars_to_keep.append('beta')
                vars_to_keep.append('SA')
            elif var == 'cRL':
                vars_to_keep.append('alpha')
                vars_to_keep.append('CT')
                vars_to_keep.append('beta')
                vars_to_keep.append('SP')
            elif var == 'cRl':
                vars_to_keep.append('alpha_PT')
                vars_to_keep.append('PT')
                vars_to_keep.append('beta_PT')
                vars_to_keep.append('SP')
            elif var == 'ell_size':
                vars_to_keep.append('press')
            # Check for variables that have underscores, ie. profile cluster
            #   average variables and local anomalies
            if '_' in var:
                # Split the prefix from the original variable (assumes an underscore split)
                split_var = var.split('_', 1)
                prefix = split_var[0]
                var_str = split_var[1]
                # If plotting a cluster variable, add the original variable (without
                #   the prefix) to the list of variables to keep
                if var in clstr_vars:
                    vars_to_keep.append(var_str)
                # Add very specific variables directly to the list, including prefixes
                if prefix in ['la', 'v1', 'v2']:
                    vars_to_keep.append(var)
                    vars_to_keep.append(var_str)
                #
            #
        # If re-running the clustering, remove 'cluster' from vars_to_keep
        if re_run_clstr:
            try:
                vars_to_keep.remove('cluster')
                vars_to_keep.remove('clst_prob')
            except:
                foo = 2
        # Add vars for the profile filters, if applicable
        scale = pp.plot_scale
        if scale == 'by_vert':
            if not isinstance(profile_filters.p_range, type(None)):
                vars_to_keep.append('press')
            if not isinstance(profile_filters.d_range, type(None)):
                vars_to_keep.append('depth')
            if not isinstance(profile_filters.iT_range, type(None)):
                vars_to_keep.append('iT')
            if not isinstance(profile_filters.CT_range, type(None)):
                vars_to_keep.append('CT')
            if not isinstance(profile_filters.PT_range, type(None)):
                vars_to_keep.append('PT')
            if not isinstance(profile_filters.SP_range, type(None)):
                vars_to_keep.append('SP')
            if not isinstance(profile_filters.SA_range, type(None)):
                vars_to_keep.append('SA')
            if profile_filters.subsample:
                vars_to_keep.append('ss_mask')
            #
        #
        # Remove duplicates
        vars_to_keep = list(set(vars_to_keep))
        # print('vars_to_keep:')
        # print(vars_to_keep)
        # exit(0)
        return vars_to_keep

################################################################################

def apply_profile_filters(arr_of_ds, vars_to_keep, profile_filters, pp):
    """
    Returns a list of pandas dataframes, one for each array in arr_of_ds with
    the filters applied to all individual profiles

    arr_of_ds           An array of datasets from a custom Data_Set object
    vars_to_keep        A list of variables to keep for the analysis
    profile_filters     A custom Profile_Filters object that contains the filters to apply
    pp                  A custom Plot_Parameters object that contains at least:
    """
    print('- Applying profile filters')
    plot_scale = pp.plot_scale
    # Add all the plotting variables
    plot_vars = pp.x_vars+pp.y_vars+[pp.clr_map]
    if not isinstance(pp.extra_args, type(None)):
        for key in pp.extra_args.keys():
            if key in ['cl_x_var', 'cl_y_var']:
                plot_vars.append(pp.extra_args[key])
    # Make an empty list
    output_dfs = []
    # What's the plot scale?
    if plot_scale == 'by_vert':
        for ds in arr_of_ds:
            # Find extra variables, if applicable
            ds = calc_extra_vars(ds, vars_to_keep)
            # Convert to a pandas data frame
            df = ds[vars_to_keep].to_dataframe()
            # Find average variables, if applicable
            # df = calc_avg_vars(df, vars_to_keep)
            # Add a notes column
            df['notes'] = ''
            #   If the m_avg_win is not None, take the moving average of the data
            if not isinstance(profile_filters.m_avg_win, type(None)):
                df = take_m_avg(df, profile_filters.m_avg_win, vars_to_keep)
            #   True/False, apply the subsample mask to the profiles
            if profile_filters.subsample:
                # `ss_mask` is null for the points that should be masked out
                # so apply mask to all the variables that were kept
                df.loc[df['ss_mask'].isnull(), vars_to_keep] = None
                # Set a new column so the ss_scheme can be found later for the title
                df['ss_scheme'] = ds.attrs['Sub-sample scheme']
            # Remove rows where the plot variables are null
            for var in plot_vars:
                if var in vars_to_keep:
                    df = df[df[var].notnull()]
            # Add expedition and instrument columns
            df['source'] = ds.Expedition
            df['instrmt'] = ds.Instrument
            ## Filters on each profile separately
            df = filter_profile_ranges(df, profile_filters, 'press', 'depth', iT_key='iT', CT_key='CT',PT_key='PT', SP_key='SP', SA_key='SA')
            ## Re-grid temperature and salinity data
            if not isinstance(profile_filters.regrid_TS, type(None)):
                # Figure out which salinity and temperature variable to re-grid
                T_var = profile_filters.regrid_TS[0]
                S_var = profile_filters.regrid_TS[2]
                # Make coarse grids for temp and salt based on given arguments
                d_temp = profile_filters.regrid_TS[1]
                d_salt = profile_filters.regrid_TS[3]
                t_arr  = np.array(df[T_var])
                s_arr  = np.array(df[S_var])
                t_min, t_max = min(t_arr), max(t_arr)
                s_min, s_max = min(s_arr), max(s_arr)
                t_grid = np.arange(t_min-d_temp, t_max+d_temp, d_temp)
                s_grid = np.arange(s_min-d_salt, s_max+d_salt, d_salt)
                # Make re-gridded arrays of temp and salt
                temp_rg = t_grid[abs(t_arr[None, :] - t_grid[:, None]).argmin(axis=0)]
                salt_rg = s_grid[abs(s_arr[None, :] - s_grid[:, None]).argmin(axis=0)]
                # Overwrite original temp and salt values with those regridded arrays
                df[T_var] = temp_rg
                df[S_var] = salt_rg
                # Note the regridding in the notes column
                df['notes'] = df['notes'] + r'$\Delta t_{rg}=$'+str(d_temp)+r', $\Delta s_{rg}=$'+str(d_salt)
            #
            # Check whether or not to take first differences
            first_dfs = pp.first_dfs
            if any(first_dfs):
                # Get list of profiles
                pfs = np.unique(np.array(df['prof_no']))
                # Take first differences in x variables
                if first_dfs[0]:
                    # Loop across each profile
                    new_dfs = []
                    for pf in pfs:
                        # Find the data for just that profile
                        data_pf = df[df['prof_no'] == pf]
                        # Loop across all x variables
                        for var in pp.x_vars:
                            # Replace the plot variable names
                            dvar = 'd_'+var
                            # Make the given variable into first differences
                            data_pf[dvar] = data_pf[var].diff()
                        new_dfs.append(data_pf)
                    # Re-form the dataframe
                    df = pd.concat(new_dfs)
                    # Replace the plot variable names
                    for i in range(len(pp.x_vars)):
                        dvar = 'd_'+pp.x_vars[i]
                        pp.x_vars[i] = dvar
                        # Remove rows with null values (one per profile because of diff())
                        df = df[df[dvar].notnull()]
                    #
                # Take finite differences in y variables
                if len(first_dfs)==2 and first_dfs[1]:
                    # Loop across each profile
                    new_dfs = []
                    for pf in pfs:
                        # Find the data for just that profile
                        data_pf = df[df['prof_no'] == pf]
                        # Loop across all y variables
                        for var in pp.y_vars:
                            # Replace the plot variable names
                            dvar = 'd_'+var
                            # Make the given variable into first differences
                            data_pf[dvar] = data_pf[var].diff()
                        new_dfs.append(data_pf)
                    # Re-form the dataframe
                    df = pd.concat(new_dfs)
                    # Repalce the plot varaible names
                    for i in range(len(pp.y_vars)):
                        dvar = 'd_'+pp.y_vars[i]
                        pp.y_vars[i] = dvar
                        # Remove rows with null values (one per profile because of diff())
                        df = df[df[dvar].notnull()]
                    #
                # 
            #
            # Check whether or not to take finite differences
            finit_dfs = pp.finit_dfs
            if any(finit_dfs):
                # Get list of profiles
                pfs = np.unique(np.array(df['prof_no']))
                # Take finite differences in x variables
                if finit_dfs[0]:
                    # Find y variable
                    y_var = pp.y_vars[0]
                    # Loop across each profile
                    new_dfs = []
                    for pf in pfs:
                        # Find the data for just that profile
                        data_pf = df[df['prof_no'] == pf]
                        # Loop across all x variables
                        for var in pp.x_vars:
                            # Replace the plot variable names
                            dvar = 'd_'+var
                            # Make the given variable into finite difference
                            data_pf[dvar] = data_pf[var].diff() / data_pf[y_var].diff()
                        new_dfs.append(data_pf)
                    # Re-form the dataframe
                    df = pd.concat(new_dfs)
                    # Replace the plot variable names
                    for i in range(len(pp.x_vars)):
                        dvar = 'd_'+pp.x_vars[i]
                        pp.x_vars[i] = dvar
                        # Remove rows with null values (one per profile because of diff())
                        df = df[df[dvar].notnull()]
                    #
                # Take finite differences in y variables
                if len(finit_dfs)==2 and finit_dfs[1]:
                    # Find x variable
                    x_var = pp.x_vars[0]
                    # Loop across each profile
                    new_dfs = []
                    for pf in pfs:
                        # Find the data for just that profile
                        data_pf = df[df['prof_no'] == pf]
                        # Loop across all y variables
                        for var in pp.y_vars:
                            # Replace the plot variable names
                            dvar = 'd_'+var
                            # Make the given variable into finite difference
                            data_pf[dvar] = data_pf[var].diff() / data_pf[x_var].diff()
                        new_dfs.append(data_pf)
                    # Re-form the dataframe
                    df = pd.concat(new_dfs)
                    # Repalce the plot varaible names
                    for i in range(len(pp.y_vars)):
                        dvar = 'd_'+pp.y_vars[i]
                        pp.y_vars[i] = dvar
                        # Remove rows with null values (one per profile because of diff())
                        df = df[df[dvar].notnull()]
                    #
                # 
            #
            # Append the filtered dataframe to the list
            output_dfs.append(df)
        #
    elif plot_scale == 'by_pf':
        for ds in arr_of_ds:
            # Convert to a pandas data frame
            df = ds[vars_to_keep].to_dataframe()
            # Add expedition and instrument columns
            df['source'] = ds.Expedition
            df['instrmt'] = ds.Instrument
            ## Filters on each profile separately
            df = filter_profile_ranges(df, profile_filters, 'press', 'depth', iT_key='iT', CT_key='CT', PT_key='PT', SP_key='SP', SA_key='SA')
            # Drop dimensions, if needed
            if 'Vertical' in df.index.names:
                # Drop duplicates along the `Time` dimension
                #   (need to make the index `Time` a column first)
                df.reset_index(level=['Time'], inplace=True)
                df.drop_duplicates(inplace=True, subset=['Time'])
            # Add a notes column
            df['notes'] = ''
            # Filter to just one entry per profile
            #   Turns out, if `vars_to_keep` doesn't have any multi-dimensional
            #       vars like 'temp' or 'salt', the resulting dataframe only has
            #       one row per profile automatically. Neat
            # Remove missing data (actually, don't because it will get rid of all
            #   data if one var isn't filled, like dt_end often isn't)
            # for var in vars_to_keep:
            #     df = df[df[var].notnull()]
            # print(df)
            output_dfs.append(df)
        #
    #
    return output_dfs

################################################################################

def take_m_avg(df, m_avg_win, vars_available):
    """
    Returns the same pandas dataframe, but with the filters provided applied to
    the data within

    df                  A pandas dataframe
    m_avg_win           The value of the moving average window in dbar
    vars_available      A list of variables available in the dataframe
    """
    print('\tIn take_m_avg(), m_avg_win:',m_avg_win)
    # Use the pandas `rolling` function to get the moving average
    #   center=True makes the first and last window/2 of the profiles are masked
    #   win_type='boxcar' uses a rectangular window shape
    #   on='press' means it will take `press` as the index column
    #   .mean() takes the average of the rolling
    #   multiplying m_avg_win by 4 because the data is in 0.25 dbar increments
    df1 = df.rolling(window=int(m_avg_win*4), center=True, win_type='boxcar', on='press').mean()
    # Put the moving average profiles for temperature, salinity, and density into the dataset
    for var in vars_available:
        if var in ['iT','ma_iT','la_iT']:
            df['ma_iT'] = df1['iT']
            if 'la_iT' in vars_available:
                df['la_iT'] = df['iT'] - df['ma_iT']
        if var in ['CT','ma_CT','la_CT']:
            df['ma_CT'] = df1['CT']
            if 'la_CT' in vars_available:
                df['la_CT'] = df['CT'] - df['ma_CT']
        if var in ['PT','ma_PT','la_PT']:
            df['ma_PT'] = df1['PT']
            if 'la_PT' in vars_available:
                df['la_PT'] = df['PT'] - df['ma_PT']
        if var in ['SP','ma_SP','la_SP']:
            df['ma_SP'] = df1['SP']
            if 'la_SP' in vars_available:
                df['la_SP'] = df['SP'] - df['ma_SP']
        if var in ['SA','ma_SA','la_SA']:
            df['ma_SA'] = df1['SA']
            if 'la_SA' in vars_available:
                df['la_SA'] = df['SA'] - df['ma_SA']
        if var in ['sigma','ma_sigma','la_sigma']:
            df['ma_sigma'] = df1['sigma']
            if 'la_sigma' in vars_available:
                df['la_sigma'] = df['sigma'] - df['ma_sigma']
            #
        #
    #
    return df

################################################################################

def filter_profile_ranges(df, profile_filters, p_key, d_key, iT_key=None, CT_key=None, PT_key=None, SP_key=None, SA_key=None):
    """
    Returns the same pandas dataframe, but with the filters provided applied to
    the data within

    df                  A pandas dataframe
    profile_filters     A custom Profile_Filters object that contains the filters to apply
    p_key               A string of the pressure variable to filter
    d_key               A string of the depth variable to filter
    iT_key              A string of the in-situ temperature variable to filter
    CT_key              A string of the conservative temperature variable to filter
    PT_key              A string of the potential temperature variable to filter
    SP_key              A string of the practical salinity variable to filter
    SA_key              A string of the absolute salinity variable to filter
    """
    #   Filter to a certain pressure range
    if not isinstance(profile_filters.p_range, type(None)):
        # Get endpoints of pressure range
        p_max = max(profile_filters.p_range)
        p_min = min(profile_filters.p_range)
        # Filter the data frame to the specified pressure range
        df = df[(df[p_key] < p_max) & (df[p_key] > p_min)]
    #   Filter to a certain depth range
    if not isinstance(profile_filters.d_range, type(None)):
        # Get endpoints of depth range
        d_max = max(profile_filters.d_range)
        d_min = min(profile_filters.d_range)
        # Filter the data frame to the specified depth range
        df = df[(df[d_key] < d_max) & (df[d_key] > d_min)]
    #   Filter to a certain in-situ temperature range
    if not isinstance(profile_filters.iT_range, type(None)):
        # Get endpoints of in-situ temperature range
        t_max = max(profile_filters.iT_range)
        t_min = min(profile_filters.iT_range)
        # Filter the data frame to the specified in-situ temperature range
        df = df[(df[iT_key] < t_max) & (df[iT_key] > t_min)]
    #   Filter to a certain conservative temperature range
    if not isinstance(profile_filters.CT_range, type(None)):
        # Get endpoints of conservative temperature range
        t_max = max(profile_filters.CT_range)
        t_min = min(profile_filters.CT_range)
        # Filter the data frame to the specified conservative temperature range
        df = df[(df[CT_key] < t_max) & (df[CT_key] > t_min)]
    #   Filter to a certain potential temperature range
    if not isinstance(profile_filters.PT_range, type(None)):
        # Get endpoints of potential temperature range
        t_max = max(profile_filters.PT_range)
        t_min = min(profile_filters.PT_range)
        # Filter the data frame to the specified potential temperature range
        df = df[(df[PT_key] < t_max) & (df[PT_key] > t_min)]
    #   Filter to a certain practical salinity range
    if not isinstance(profile_filters.SP_range, type(None)):
        # Get endpoints of practical salinity range
        s_max = max(profile_filters.SP_range)
        s_min = min(profile_filters.SP_range)
        # Filter the data frame to the specified practical salinity range
        df = df[(df[SP_key] < s_max) & (df[SP_key] > s_min)]
    #   Filter to a certain absolute salinity range
    if not isinstance(profile_filters.SA_range, type(None)):
        # Get endpoints of absolute salinity range
        s_max = max(profile_filters.SA_range)
        s_min = min(profile_filters.SA_range)
        # Filter the data frame to the specified absolute salinity range
        df = df[(df[SA_key] < s_max) & (df[SA_key] > s_min)]
    return df

################################################################################

def calc_extra_vars(ds, vars_to_keep):
    """
    Takes in an xarray object and a list of variables and, if there are extra
    variables to calculate, it will add those to the xarray

    ds                  An xarray from the arr_of_ds of a custom Data_Set object
    vars_to_keep        A list of variables to keep for the analysis
    """
    # Check for variables to calculate
    if 'aiT' in vars_to_keep:
        ds['aiT'] = ds['alpha'] * ds['iT']
    if 'aCT' in vars_to_keep:
        ds['aCT'] = ds['alpha'] * ds['CT']
    if 'aPT' in vars_to_keep:
        ds['aPT'] = ds['alpha_PT'] * ds['PT']
    if 'BSP' in vars_to_keep:
        ds['BSP'] = ds['beta'] * ds['SP']
    if 'BSt' in vars_to_keep:
        ds['BSt'] = ds['beta_PT'] * ds['SP']
    if 'BSA' in vars_to_keep:
        ds['BSA'] = ds['beta'] * ds['SA']
    # Check for variables with an underscore
    for this_var in vars_to_keep:
        if '_' in this_var and this_var not in ['prof_no','clst_prob']:
            # Split the prefix from the original variable (assumes an underscore split)
            split_var = this_var.split('_', 1)
            prefix = split_var[0]
            var = split_var[1]
            # print('prefix:',prefix,'- var:',var)
            if prefix == 'la':
                # Calculate the local anomaly of this variable
                ds[this_var] = ds[var] - ds['ma_'+var]
            #
            if prefix == 'v1':
                # Calculate the local anomaly of this variable
                ds[this_var] = ds[var] * 1.1
            #
            if prefix == 'v2':
                # Calculate the local anomaly of this variable
                ds[this_var] = (ds[var] - ds['ma_'+var]) * 1.1
            #
        #
    return ds

################################################################################

def get_axis_labels(pp, var_attr_dicts):
    """
    Using the dictionaries of variable attributes, this adds the xlabel and ylabel
    attributes to the plot parameters object passed in (this replaced the
    `get_axis_label` function)

    pp                  A custom Plot_Parameters object
    var_attr_dicts      A list of dictionaries containing the attribute info
                            about all the variables in the netcdf file
    """
    # Get x axis labels
    if not isinstance(pp.x_vars, type(None)):
        if len(pp.x_vars) > 0:
            try:
                pp.xlabels[0] = var_attr_dicts[0][pp.x_vars[0]]['label']
            except:
                pp.xlabels[0] = get_axis_label(pp.x_vars[0], var_attr_dicts)
        else:
            pp.xlabels[0] = None
        if len(pp.x_vars) == 2:
            try:
                pp.xlabels[1] = var_attr_dicts[0][pp.x_vars[1]]['label']
            except:
                pp.xlabels[1] = get_axis_label(pp.x_vars[1], var_attr_dicts)
        else:
            pp.xlabels[1] = None
    else:
        pp.xlabels[0] = None
        pp.xlabels[1] = None
    # Get y axis labels
    if not isinstance(pp.y_vars, type(None)):
        if len(pp.y_vars) > 0:
            try:
                pp.ylabels[0] = var_attr_dicts[0][pp.y_vars[0]]['label']
            except:
                pp.ylabels[0] = get_axis_label(pp.y_vars[0], var_attr_dicts)
        else:
            pp.ylabels[0] = None
        if len(pp.y_vars) == 2:
            try:
                pp.ylabels[1] = var_attr_dicts[0][pp.y_vars[1]]['label']
            except:
                pp.ylabels[1] = get_axis_label(pp.y_vars[1], var_attr_dicts)
        else:
            pp.ylabels[1] = None
    else:
        pp.ylabels[0] = None
        pp.ylabels[1] = None
    # Get colormap label
    if not isinstance(pp.clr_map, type(None)):
        try:
            pp.clabel = var_attr_dicts[0][pp.clr_map]['label']
        except:
            pp.clabel = get_axis_label(pp.clr_map, var_attr_dicts)

    # Check for first differences of variables
    if any(pp.first_dfs):
        if pp.first_dfs[0]:
            if not isinstance(pp.xlabels[0], type(None)):
                pp.xlabels[0] = 'First difference in '+pp.xlabels[0]
            if not isinstance(pp.xlabels[1], type(None)):
                pp.xlabels[1] = 'First difference in '+pp.xlabels[1]
        if len(pp.first_dfs)==2 and pp.first_dfs[1]:
            if not isinstance(pp.ylabels[0], type(None)):
                pp.ylabels[0] = 'First difference in '+pp.ylabels[0]
            if not isinstance(pp.ylabels[1], type(None)):
                pp.ylabels[1] = 'First difference in '+pp.ylabels[1]
    if any(pp.finit_dfs):
        if pp.finit_dfs[0]:
            if not isinstance(pp.xlabels[0], type(None)):
                pp.xlabels[0] = r'Finite difference $\Delta$'+pp.xlabels[0]+'('+pp.ylabels[0]+')'
            if not isinstance(pp.xlabels[1], type(None)):
                pp.xlabels[1] = r'Finite difference $\Delta$'+pp.xlabels[1]+'('+pp.ylabels[0]+')'
        if len(pp.finit_dfs)==2 and pp.finit_dfs[1]:
            if not isinstance(pp.ylabels[0], type(None)):
                pp.ylabels[0] = r'First difference $\Delta$'+pp.ylabels[0]+'('+pp.xlabels[0]+')'
            if not isinstance(pp.ylabels[1], type(None)):
                pp.ylabels[1] = r'First difference $\Delta$'+pp.ylabels[1]+'('+pp.xlabels[0]+')'
    if pp.plot_type == 'map':
        pp.plot_scale = 'by_pf'
        pp.x_vars  = ['lon']
        pp.y_vars  = ['lat']
        pp.xlabels = [None, None]
        pp.ylabels = [None, None]
    #
    return pp

################################################################################

def get_axis_label(var_key, var_attr_dicts):
    """
    Takes in a variable name and returns a nicely formatted string for an axis label

    var_key             A string of the variable name
    var_attr_dicts      A list of dictionaries of variable attributes
    """
    # Check for certain modifications to variables,
    #   check longer strings first to avoid mismatching
    # Check for profile cluster average variables
    if 'pca_' in var_key:
        # Take out the first 4 characters of the string to leave the original variable name
        var_str = var_key[4:]
        return 'PCA '+ var_attr_dicts[0][var_str]['label']
        # return 'Profile cluster average of '+ var_attr_dicts[0][var_str]['label']
    # Check for profile cluster span variables
    if 'pcs_' in var_key:
        # Take out the first 4 characters of the string to leave the original variable name
        var_str = var_key[4:]
        return 'Profile cluster span of '+ var_attr_dicts[0][var_str]['label']
    # Check for cluster mean-centered variables
    elif 'cmc_' in var_key:
        # Take out the first 4 characters of the string to leave the original variable name
        var_str = var_key[4:]
        return 'Cluster mean-centered '+ var_attr_dicts[0][var_str]['label']
    # Check for local anomaly variables
    elif 'la_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[3:]
        return r"$\Theta'$ ($^\circ$C)"
        # return 'Local anomaly of '+ var_attr_dicts[0][var_str]['label']
    # Check for local anomaly variables
    elif 'v1_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[3:]
        return 'Stretched '+ var_attr_dicts[0][var_str]['label']
    # Check for local anomaly variables
    elif 'v2_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[3:]
        return 'Stretched local anomaly of '+ var_attr_dicts[0][var_str]['label']
    # Check for cluster average variables
    elif 'ca_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[3:]
        return 'Cluster average of '+ var_attr_dicts[0][var_str]['label']
    # Check for cluster span variables
    elif 'cs_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[3:]
        return 'Cluster span of '+ var_attr_dicts[0][var_str]['label']
    # Check for cluster min/max variables
    elif 'cmm_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[4:]
        return 'Cluster min/max of '+ var_attr_dicts[0][var_str]['label']
    # Check for normalized inter-cluster range variables
    elif 'nir_' in var_key:
        # Take out the first 3 characters of the string to leave the original variable name
        var_str = var_key[4:]
        return r'Normalized inter-cluster range $IR_{S_P}$'
        # return r'Normalized inter-cluster range $IR$ of '+ var_attr_dicts[0][var_str]['label']
    # Build dictionary of axis labels
    ax_labels = {
                 'hist':r'Occurrences',
                 'aiT':r'$\alpha T$',
                 'aCT':r'$\alpha \Theta$',
                 'aPT':r'$\alpha \theta$',
                 'BSP':r'$\beta S_P$',
                 'BSt':r'$\beta_{PT} S_P$',
                 'BSA':r'$\beta S_A$',
                #  'distance':r'Along-path distance (km)',
                 'm_pts':r'Minimum density threshold $m_{pts}$',
                 'DBCV':'Relative validity measure (DBCV)',
                 'n_clusters':'Number of clusters',
                 'cRL':r'Lateral density ratio $R_L$',
                 'cRl':r'Lateral density ratio $R_L$ with $\theta$'
                }
    if var_key in ax_labels.keys():
        return ax_labels[var_key]
    else:
        return 'None'

################################################################################
# Admin summary functions ######################################################
################################################################################

def txt_summary(groups_to_summarize, filename=None):
    """
    Takes in a list of Analysis_Group objects. Analyzes and outputs summary info
    for each one.

    groups_to_plot      A list of Analysis_Group objects
    filename            The file to write the output to. If none, it will print
                            the output to the console
    """
    # Loop over all Analysis_Group objects
    i = 0
    for a_group in groups_to_summarize:
        i += 1
        # Concatonate all the pandas data frames together
        df = pd.concat(a_group.data_frames)
        # Find the total number of profiles
        n_profs = len(np.unique(np.array(df['prof_no'], dtype=type(''))))
        # Put all of the lines together
        lines = ["Group "+str(i)+":",
                 # "\t"+add_std_title(a_group),
                 print_global_variables(a_group.data_set),
                 "\t\tVariables to keep:",
                 "\t\t"+str(a_group.vars_to_keep),
                 "\t\tNumber of profiles: ",
                 "\t\t\t"+str(n_profs),
                 "\t\tNumber of data points: ",
                 "\t\t\t"+str(len(df))
                 ]
        # Find the ranges of the variables available
        print('working on ranges')
        var_attr_dict = a_group.data_set.var_attr_dicts[0]
        for var in a_group.vars_to_keep:
            if var in var_attr_dict.keys():
                lines.append(report_range(df, var, var_attr_dict[var]))
            #
        #
        # Check for SDA variables
        if 'SDA_n_ml_dc' in a_group.vars_to_keep and 'SDA_n_ml_sf' in a_group.vars_to_keep and 'SDA_n_gl_dc' in a_group.vars_to_keep and 'SDA_n_gl_sf' in a_group.vars_to_keep:# and 'SDA_ml_press' in a_group.vars_to_keep and 'SDA_gl_h' in a_group.vars_to_keep:
            print('working on SDA part')
            lines.append("\n\t\t-- SDA output summary --")
            lines.append(summarize_SDA_output(df, n_profs))
        lines.append("\n")
    #
    if filename != None:
        print('Writing summary file to',filename)
        with open(filename, 'w') as f:
            f.write('\n'.join(lines))
    else:
        for line in lines:
            print(line)

################################################################################

def print_profile_filters(pfs):
    """
    Returns a human readable string of all the profile filters in the given object

    pfs         A Profile_Filters object
    """
    return_string = ''
    if not isinstance(pfs.p_range, type(None)): 
        return_string += ('Pressure range: ['+str(pfs.p_range[0])+', '+str(pfs.p_range[1])+'] ')
    if not isinstance(pfs.d_range, type(None)): 
        return_string += ('Depth range: ['+str(pfs.d_range[0])+', '+str(pfs.d_range[1])+'] ')
    if not isinstance(pfs.iT_range, type(None)): 
        return_string += ('iT range: ['+str(pfs.iT_range[0])+', '+str(pfs.iT_range[1])+'] ')
    if not isinstance(pfs.CT_range, type(None)): 
        return_string += ('CT range: ['+str(pfs.CT_range[0])+', '+str(pfs.CT_range[1])+'] ')
    if not isinstance(pfs.PT_range, type(None)): 
        return_string += ('PT range: ['+str(pfs.PT_range[0])+', '+str(pfs.PT_range[1])+'] ')
    if not isinstance(pfs.SP_range, type(None)): 
        return_string += ('SP range: ['+str(pfs.SP_range[0])+', '+str(pfs.SP_range[1])+'] ')
    if not isinstance(pfs.SA_range, type(None)): 
        return_string += ('SA range: ['+str(pfs.SA_range[0])+', '+str(pfs.SA_range[1])+'] ')
    if not isinstance(pfs.lon_range, type(None)): 
        return_string += ('Longitude range: ['+str(pfs.lon_range[0])+', '+str(pfs.lon_range[1])+'] ')
    if not isinstance(pfs.lat_range, type(None)): 
        return_string += ('Latitude range: ['+str(pfs.lat_range[0])+', '+str(pfs.lat_range[1])+'] ')
    if pfs.lt_pCT_max: 
        return_string += ('Pressures less than p(CT_max) ')
    if pfs.subsample: 
        return_string += ('Subsampled ')
    if not isinstance(pfs.regrid_TS, type(None)): 
        return_string += ('Regrid: d'+str(pfs.regrid_TS[0])+'='+str(pfs.regrid_TS[1])+', d'+str(pfs.regrid_TS[2])+'='+str(pfs.regrid_TS[3])+' ')
    if not isinstance(pfs.m_avg_win, type(None)): 
        return_string += ('Moving average window: ['+str(pfs.m_avg_win)+' ')
    return return_string

################################################################################

def report_range(df, key, name_and_units):
    """
    Reports the range of the key value given for the given dataframe

    df          A pandas data frame
    key         A string of the column to search in the data frame
    name_and_units      A dictionary containing:
        name_str    A string of the long version of the column data type
        units       A string of the units of the data type
    """
    name_str  = name_and_units['long_name']
    try:
        units = name_and_units['units']
    except:
        # The above will fail if no units are given, which currently happens for
        #   the time variables in the GDS netcdfs because I needed to sort them
        print('Missing units for:',key)
        units = 'YYYY-MM-DD HH:MM:SS'
    # Find this range, making sure to skip null (nan) values
    this_arr = df[df[key].notnull()][key]
    if len(this_arr) > 1:
        this_min = min(this_arr)
        this_max = max(this_arr)
        if key != 'dt_start' and key != 'dt_end':
            # This line formats all the numbers to have just 2 digits after the decimal point
            return "\t\t"+key+': '+name_str+" \n\t\t\trange: %.2f"%(this_max-this_min)+" "+units+"\n\t\t\t%.2f"%this_min+" to %.2f"%this_max+" "+units
            # This line outputs the numbers with all available digits after the decimal point
            # return "\t\t"+name_str+" range: "+str(this_max-this_min)+" "+units+"\n\t\t\t"+str(this_min)+" to "+str(this_max)+" "+units
        else:
            return "\t\t"+key+': '+name_str+" \n\t\t\trange: \n\t\t\t"+str(this_min)+" to "+str(this_max)+" "+units
    else:
        return "\t\t"+key+': '+name_str+" \n\t\t\trange: N/A \n\t\t\tN/A"

################################################################################

def print_global_variables(Dataset):
    """
    Prints out the global variables of the data sources in a Data_Set object

    Dataset         A custom Data_Set object
    """
    lines = ''
    for ds in Dataset.arr_of_ds:
        for attr in ds.attrs:
            print('\t'+str(attr)+': '+str(ds.attrs[attr]))
            lines = lines+'\t'+str(attr)+': '+str(ds.attrs[attr])+'\n'
            # if attr in ['Creation date', 'Last modified', 'Last modification', 'Sub-sample scheme']:
            #     print('\t'+attr+': '+ds.attrs[attr])
            # else:
            #     lines = lines+'\t'+attr+': '+ds.attrs[attr]+'\n'
    return lines

################################################################################

def find_max_distance(groups_to_analyze):
    """
    Reports the maximum distance between any two profiles in the given dataframe

    groups_to_analyze   A list of Analysis_Group objects
                        Each Analysis_Group contains a list of dataframes
    """
    for ag in groups_to_analyze:
        for df in ag.data_frames:
            # Drop duplicates in lat/lon to have just one row per profile
            df.drop_duplicates(subset=['lon','lat'], keep='first', inplace=True)
            # Get source and instrument for this df
            this_source = df['source'].values[0]
            this_instrmt = df['instrmt'].values[0]
            # Get arrays of values
            prof_nos = df['prof_no'].values
            lon_vals = df['lon'].values
            lat_vals = df['lat'].values
            # Find the number of profiles
            n_pfs = len(prof_nos)
            # Placeholders to store the largest span between any two profiles
            max_span = 0
            ms_i = None
            ms_i_latlon = None
            ms_j = None
            ms_j_latlon = None
            # Nested loop over the profiles, only the upper triangle
            #   Don't need to double-count, comparing profiles twice
            for i in range(n_pfs):
                for j in range(i+1,n_pfs):
                    i_lat_lon = (lat_vals[i], lon_vals[i])
                    j_lat_lon = (lat_vals[j], lon_vals[j])
                    this_span = geodesic(i_lat_lon, j_lat_lon).km
                    if this_span > max_span:
                        max_span = this_span
                        ms_i = prof_nos[i] 
                        ms_i_latlon = i_lat_lon
                        ms_j = prof_nos[j]
                        ms_j_latlon = j_lat_lon
                    # print('i:',i,'j:',j,'i_pf:',prof_nos[i],'j_pf:',prof_nos[j],'this_span:',this_span,'max_span:',max_span)
                #
            #
            print('For',this_source,this_instrmt,'max_span:',max_span,'km between profiles',ms_i,'at',ms_i_latlon,'and',ms_j,'at',ms_j_latlon)

################################################################################
# Admin plotting functions #####################################################
################################################################################

def make_figure(groups_to_plot, filename=None, use_same_x_axis=None, use_same_y_axis=None, row_col_list=None):
    """
    Takes in a list of Analysis_Group objects, one for each subplot. Determines
    the needed arrangement of subplots, then passes one Analysis_Group object to
    each axis for plotting

    groups_to_plot  A list of Analysis_Group objects, one for each subplot
                    Each Analysis_Group contains the info to create each subplot
    filename        The filename in which to save the figure
                        only accepts .png or .pickle filenames
    use_same_x_axis     True/False whether to force subplots to share x axis
                        ranges if they have the same variable
    use_same_y_axis     True/False whether to force subplots to share y axis
                        ranges if they have the same variable
    row_col_list        [rows, cols, f_ratio, f_size], if none given, will use 
                        the defaults specified below in n_row_col_dict
    """
    # Define number of rows and columns based on number of subplots
    #   key: number of subplots, value: (rows, cols, f_ratio, f_size)
    n_row_col_dict = {'1':[1,1, 0.8, 1.25], '2':[1,2, 0.5, 1.25], '2.5':[2,1, 0.8, 1.25],
                      '3':[3,1, 0.3, 1.70], '4':[2,2, 0.7, 2.00],
                      '5':[2,3, 0.5, 2.00], '6':[2,3, 0.5, 2.00],
                      '7':[2,4, 0.4, 2.00], '8':[2,4, 0.4, 2.00],
                      '9':[3,3, 0.8, 2.00]}
    # Figure out what layout of subplots to make
    n_subplots = len(groups_to_plot)
    x_keys = []
    y_keys = []
    for group in groups_to_plot:
        pp = group.plt_params
        u_x_vars = np.unique(pp.x_vars)
        u_x_vars = np.delete(u_x_vars, np.where(u_x_vars == None))
        u_y_vars = np.unique(pp.y_vars)
        u_y_vars = np.delete(u_y_vars, np.where(u_y_vars == None))
        if len(u_x_vars) > 0:
            for var in u_x_vars:
                x_keys.append(var)
        if len(u_y_vars) > 0:
            for var in u_y_vars:
                y_keys.append(var)
            #
        #
    #
    if len(x_keys) > 1 and isinstance(use_same_x_axis, type(None)):
        # If all the x_vars are the same, share the x axis between subplots
        if all(x == x_keys[0] for x in x_keys):
            use_same_x_axis = True
            print('\t- Set share_x_axis to True')
        else:
            use_same_x_axis = False
        #
    else:
        use_same_x_axis = False
    if len(y_keys) > 1 and isinstance(use_same_y_axis, type(None)):
        # If all the y_vars are the same, share the y axis between subplots
        if all(y == y_keys[0] for y in y_keys):
            use_same_y_axis = True
            print('\t- Set share_y_axis to True')
        else:
            use_same_y_axis = False
        #
    else:
        use_same_y_axis = False
    tight_layout_h_pad = 1.0 #None
    tight_layout_w_pad = 1.0 #None
    if n_subplots == 1:
        if isinstance(row_col_list, type(None)):
            rows, cols, f_ratio, f_size = n_row_col_dict[str(n_subplots)]
        else:
            rows, cols, f_ratio, f_size = row_col_list
        n_subplots = int(np.floor(n_subplots))
        fig, ax = set_fig_axes([1], [1], fig_ratio=f_ratio, fig_size=f_size)
        xlabel, ylabel, plt_title, ax, invert_y_axis = make_subplot(ax, groups_to_plot[0], fig, 111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        # Invert y-axis if specified
        if invert_y_axis:
            ax.invert_yaxis()
            print('\t- Inverting y-axis')
        # If axes limits given, limit axes
        if not isinstance(pp.ax_lims, type(None)):
            try:
                ax.set_xlim(pp.ax_lims['x_lims'])
                print('\t- Set x_lims to',pp.ax_lims['x_lims'])
            except:
                foo = 2
            try:
                ax.set_ylim(pp.ax_lims['y_lims'])
                print('\t- Set y_lims to',pp.ax_lims['y_lims'])
            except:
                foo = 2
            # 
        # If called for, add a grid to the plot by default
        print('\t- Adding grid lines:',pp.add_grid)
        if pp.add_grid:
            ax.grid(color=std_clr, linestyle='--', alpha=grid_alpha)
        ax.set_title(plt_title)
    elif n_subplots > 1 and n_subplots < 10:
        if isinstance(row_col_list, type(None)):
            rows, cols, f_ratio, f_size = n_row_col_dict[str(n_subplots)]
        else:
            rows, cols, f_ratio, f_size = row_col_list
        n_subplots = int(np.floor(n_subplots))
        fig, axes = set_fig_axes([1]*rows, [1]*cols, fig_ratio=f_ratio, fig_size=f_size, share_x_axis=use_same_x_axis, share_y_axis=use_same_y_axis)
        for i in range(n_subplots):
            print('- Subplot '+string.ascii_lowercase[i])
            if rows > 1 and cols > 1:
                i_ax = (i//cols,i%cols)
            else:
                i_ax = i
            ax_pos = int(str(rows)+str(cols)+str(i+1))
            xlabel, ylabel, plt_title, ax, invert_y_axis = make_subplot(axes[i_ax], groups_to_plot[i], fig, ax_pos)
            if use_same_x_axis:
                # If on the top row
                if i < cols == 0:
                    tight_layout_h_pad = -1
                    ax.set_title(plt_title)
                # If in the bottom row
                elif i >= n_subplots-cols:
                    ax.set_xlabel(xlabel)
            else:
                ax.set_title(plt_title)
                ax.set_xlabel(xlabel)
            ax.set_title(plt_title)
            if use_same_y_axis:
                # If in the far left column
                if i%cols == 0:
                    tight_layout_w_pad = -1
                    ax.set_ylabel(ylabel)
                # Invert y-axis if specified
                if i == 0 and invert_y_axis:
                    # ax.invert_yaxis()
                    print('\t- Inverting y-axis')
            else:
                ax.set_ylabel(ylabel)
                # Invert y-axis if specified
                if invert_y_axis:
                    ax.invert_yaxis()
                    print('\t- Inverting y-axis-')
            # If axes limits given, limit axes
            this_ax_pp = groups_to_plot[i].plt_params
            if not isinstance(this_ax_pp.ax_lims, type(None)):
                try:
                    ax.set_xlim(this_ax_pp.ax_lims['x_lims'])
                    print('\t- Set x_lims to',this_ax_pp.ax_lims['x_lims'])
                except:
                    foo = 2
                try:
                    ax.set_ylim(this_ax_pp.ax_lims['y_lims'])
                    print('\t- Set y_lims to',this_ax_pp.ax_lims['y_lims'])
                except:
                    foo = 2
                #
            # If called for, add a grid to the plot by default
            print('\t- Adding grid lines:',this_ax_pp.add_grid)
            if this_ax_pp.add_grid:
                ax.grid(color=std_clr, linestyle='--', alpha=grid_alpha)
            # Label subplots a, b, c, ...
            ax.text(-0.1, -0.1, '('+string.ascii_lowercase[i]+')', transform=ax.transAxes, size=mpl.rcParams['axes.labelsize'], fontweight='bold')
        # Turn off unused axes
        if n_subplots < (rows*cols):
            for i in range(rows*cols-1, n_subplots-1, -1):
                axes[i//cols,i%cols].set_axis_off()
    else:
        print('Too many subplots')
        exit(0)
    #
    plt.tight_layout(pad=0.4, h_pad=tight_layout_h_pad, w_pad=tight_layout_w_pad)
    #
    if filename != None:
        print('- Saving figure to outputs/'+filename)
        if '.png' in filename:
            plt.savefig('outputs/'+filename, dpi=1000, transparent=True)
        elif '.pickle' in filename:
            pl.dump(fig, open('outputs/'+filename, 'wb'))
        else:
            print('File extension not recognized in',filename)
    else:
        print('- Displaying figure')
        plt.show()

################################################################################
# Formatting plotting functions ################################################
################################################################################

def set_fig_axes(heights, widths, fig_ratio=0.5, fig_size=1, share_x_axis=None, share_y_axis=None, prjctn=None):
    """
    Creates fig and axes objects based on desired heights and widths of subplots
    Ex: if widths=[1,5], there will be 2 columns, the 1st 1/5 the width of the 2nd

    heights      array of integers for subplot height ratios, len=rows
    widths       array of integers for subplot width  ratios, len=cols
    fig_ratio    ratio of height to width of overall figure
    fig_size     size scale factor, 1 changes nothing, 2 makes it very big
    share_x_axis bool whether the subplots should share their x axes
    share_y_axis bool whether the subplots should share their y axes
    prjctn       projection type for the subplots
    """
    # Set aspect ratio of overall figure
    w, h = mpl.figure.figaspect(fig_ratio)
    # Find rows and columns of subplots
    rows = len(heights)
    cols = len(widths)
    # This dictionary makes each subplot have the desired ratios
    # The length of heights will be nrows and likewise len(widths)=ncols
    plot_ratios = {'height_ratios': heights,
                   'width_ratios': widths}
    # Determine whether to share x or y axes
    if share_x_axis == None and share_y_axis == None:
        if rows == 1 and cols != 1: # if only one row, share y axis
            share_x_axis = False
            share_y_axis = True
        elif rows != 1 and cols == 1: # if only one column, share x axis
            share_x_axis = True
            share_y_axis = False
        else:                       # otherwise, forget about it
            share_x_axis = False
            share_y_axis = False
    elif share_y_axis == None:
        share_y_axis = False
        print('\tSet share_y_axis to', share_y_axis)
    elif share_x_axis == None:
        share_x_axis = False
        print('\tSet share_x_axis to', share_x_axis)
    # Set ratios by passing dictionary as 'gridspec_kw', and share y axis
    fig, axes = plt.subplots(figsize=(w*fig_size,h*fig_size), nrows=rows, ncols=cols, gridspec_kw=plot_ratios, sharex=share_x_axis, sharey=share_y_axis, subplot_kw=dict(projection=prjctn))
    # Set ticklabel format for all axes
    if (rows+cols)>2:
        for ax in axes.flatten():
            ax.ticklabel_format(style='sci', scilimits=sci_lims, useMathText=True)
    else:
        axes.ticklabel_format(style='sci', scilimits=sci_lims, useMathText=True)
    # Turn off excess ticks if the axes are shared
    if share_x_axis and cols == 1:
        print('\t- Turning off excess x axes')
        for ax in axes.flatten()[0:-2]:
            plt.setp(ax.get_xticklabels(), visible=False)
            # plt.setp(ax.get_xlabel(), visible=False)
    return fig, axes

################################################################################

def add_std_title(a_group):
    """
    Adds in standard information to this subplot's title, as appropriate

    a_group         A Analysis_Group object containing the info to create a subplot
    """
    # Check for an override title
    if not isinstance(a_group.plot_title, type(None)):
        return a_group.plot_title
    else:
        plt_title = ''
    # Determine a prefix, if any, for the plot title
    pfs = a_group.profile_filters
    if pfs.subsample:
        # Get a string of the subsampling scheme
        ss_scheme = np.unique(np.array(a_group.data_frames[0]['ss_scheme']))[0]
        plt_title = 'Subsampled w/ ' + ss_scheme + ':'
    sources_dict = a_group.data_set.sources_dict
    # Get a list of the sources in the plot
    sources = list(sources_dict.keys())
    # Get just the unique sources
    sources = np.unique(sources)
    if len(sources) > 3:
        plt_title += 'Many sources'
    else:
        # Note all the data sources in the title
        plt_title += sources[0].replace('_',' ')
        for i in range(1,len(sources)):
            plt_title += ', '
            plt_title += sources[i].replace('_',' ')
    #
    return plt_title

################################################################################

def add_std_legend(ax, data, x_key):
    """
    Adds in standard information to this subplot's legend, as appropriate

    ax          The axis on which to add the legend
    data        A pandas data frame which is plotted on this axis
    x_key       The string of the column name for the x data from the dataframe
    """
    # Add legend to report the total number of points and notes on the data
    n_pts_patch  = mpl.patches.Patch(color='none', label=str(len(data[x_key]))+' points')
    notes_string = ''.join(data.notes.unique())
    # Only add the notes_string if it contains something
    if len(notes_string) > 1:
        notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
        ax.legend(handles=[n_pts_patch, notes_patch])
    else:
        ax.legend(handles=[n_pts_patch])

################################################################################

def get_var_color(var):
    """
    Takes in the variable name and for the variable to colormap

    var         A string of the variable name to get the color for
    """
    # Check whether taking first or finite differences
    if 'd_' in var:
        var = var[2:]
    # Build dictionary of available colors
    vars = {
            'aiT':'salmon',
            'aCT':'salmon',
            'aPT':'salmon',
            'BSP':'darkturquoise',
            'BSt':'darkturquoise',
            'BSA':'darkturquoise'
            }
    if var in ['entry', 'prof_no', 'dt_start', 'dt_end', 'lon', 'lat', 'press', 'depth', 'og_press', 'og_depth']:
        return std_clr
    elif var in ['iT', 'CT', 'PT', 'alpha']:
        return 'lightcoral'
    elif var in ['SP', 'SA', 'beta']:
        return 'cornflowerblue'
    elif var in ['sigma']:
        return 'mediumpurple'
    elif var in ['ma_iT', 'ma_CT', 'ma_PT', 'la_iT', 'la_CT', 'la_PT']:
        return 'tab:red'
    elif var in ['ma_SP', 'ma_SA', 'la_SP', 'la_SA']:
        return 'tab:blue'
    elif var in ['ma_sigma', 'la_sigma']:
        return 'tab:purple'
    if var in vars.keys():
        return vars[var]
    else:
        return None

################################################################################

def get_color_map(cmap_var):
    """
    Takes in the variable name and for the variable to colormap

    cmap_var    A string of the variable name for the colormap
    """
    # Build dictionary of axis labels
    cmaps = {'entry':'plasma',
             'prof_no':'plasma',
             'dt_start':'viridis',
             'dt_end':'viridis',
             'lon':'YlOrRd',
             'lat':'PuBuGn',
             'R_rho':'ocean',
             'density_hist':'inferno'
             }
    if cmap_var in ['press', 'depth']:
        return 'cividis'
    elif cmap_var in ['iT', 'CT', 'PT', 'alpha', 'aiT', 'aCT', 'aPT', 'ma_iT', 'ma_CT', 'ma_PT', 'la_iT', 'la_CT', 'la_PT']:
        return 'Reds'
    elif cmap_var in ['SP', 'SA', 'beta', 'BSP', 'BSt', 'BSA', 'ma_SP', 'ma_SA', 'la_SP', 'la_SA']:
        return 'Blues'
    elif cmap_var in ['sigma', 'ma_sigma', 'la_sigma']:
        return 'Purples'
    elif cmap_var in ['BL_yn', 'up_cast', 'ss_mask']:
        cmap = mpl.colors.ListedColormap(['green'])
        cmap.set_bad(color='red')
        return cmap
    if cmap_var in cmaps.keys():
        return cmaps[cmap_var]
    else:
        return None

################################################################################

def format_datetime_axes(x_key, y_key, ax, tw_x_key=None, tw_ax_y=None, tw_y_key=None, tw_ax_x=None):
    """
    Formats any datetime axes to show actual dates, as appropriate

    x_key       The string of the name for the x data on the main axis
    y_key       The string of the name for the y data on the main axis
    ax          The main axis on which to format
    tw_x_key    The string of the name for the x data on the twin axis
    tw_ax_y     The twin y axis on which to format
    tw_y_key    The string of the name for the y data on the twin axis
    tw_ax_x     The twin x axis on which to format
    """
    loc = mpl.dates.AutoDateLocator()
    if x_key in ['dt_start', 'dt_end']:
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))
    if y_key in ['dt_start', 'dt_end']:
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))
    if tw_x_key in ['dt_start', 'dt_end']:
        tw_ax_y.xaxis.set_major_locator(loc)
        tw_ax_y.xaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))
    if tw_y_key in ['dt_start', 'dt_end']:
        tw_ax_x.yaxis.set_major_locator(loc)
        tw_ax_x.yaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))

################################################################################

def format_sci_notation(x, ndp=2):
    """
    Formats a number into scientific notation

    x       The number to format
    ndp     The number of decimal places
    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    # Check to see whether it's outside the scientific notation exponent limits
    if int(e) < min(sci_lims) or int(e) > max(sci_lims):
        return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))
    else:
        return r'{x:0.{ndp:d}f}'.format(x=x, ndp=ndp)

################################################################################

def add_h_scale_bar(ax, ax_lims, unit="", tw_clr=False):
    """
    Adds a horizontal scale bar to a plot

    ax          The axis on which to add the legend
    ax_lims     The limits on the horizontal axis range
    unit        A string of the units of the horizontal axis
    tw_clr      The color of the twin axis, if necessary
    """
    h_bar_clr = std_clr
    # Find distance between ticks
    ax_ticks = np.array(ax.get_xticks())
    # Use the ticks to decide on the bar length
    bar_len = min(np.diff(ax_ticks))/2
    # Find the bounds of the plot axis
    if isinstance(ax_lims, type(None)):
        x_bnds = ax.get_xbound()
        y_bnds = ax.get_ybound()
    else:
        if 'x_lims' in ax_lims.keys():
            x_bnds = ax_lims['x_lims']
        else:
            x_bnds = ax.get_xbound()
        if 'y_lims' in ax_lims.keys():
            y_bnds = ax_lims['y_lims']
        else:
            y_bnds = ax.get_ybound()
    # print('x_bnds:',x_bnds,'y_bnds:',y_bnds)
    x_span = abs(max(x_bnds)-min(x_bnds))
    y_span = abs(max(y_bnds)-min(y_bnds))
    # Adjust bar length if too small
    # if x_span/bar_len > 17:
    #     bar_len = bar_len*2
    # Define end points for the scale bar
    x_pad = bar_len/3.8
    sb_x_min = min(x_bnds) + x_pad 
    sb_x_max = sb_x_min + bar_len
    if tw_clr==False:
        sb_y_val = max(y_bnds) - y_span/8
    else:
        sb_y_val = max(y_bnds) - y_span/4
        h_bar_clr = tw_clr
    # Plot scale bar
    ax.plot([sb_x_min,sb_x_max], [sb_y_val,sb_y_val], color=h_bar_clr)
    # Plot bar ends
    ax.scatter([sb_x_min,sb_x_max], [sb_y_val,sb_y_val], color=h_bar_clr, marker='|')
    # Annotate the bar length
    ax.annotate(r'%.3f'%(bar_len)+unit, xy=(sb_x_min+(sb_x_max-sb_x_min)/2, sb_y_val+y_span/16), ha='center', xycoords='data', color=h_bar_clr, zorder=12)
    # Turn off the x tick labels
    ax.set_xticklabels([])

################################################################################
# Main plotting function #######################################################
################################################################################

def make_subplot(ax, a_group, fig, ax_pos):
    """
    Takes in an Analysis_Group object which has the data and plotting parameters
    to produce a subplot. Returns the x and y labels and the subplot title

    ax              The axis on which to make the plot
    a_group         A Analysis_Group object containing the info to create this subplot
    fig             The figure in which ax is contained
    ax_pos          A tuple of the ax (rows, cols, linear number of this subplot)
    """
    # Get relevant parameters for the plot
    pp = a_group.plt_params
    plot_type = pp.plot_type
    clr_map   = pp.clr_map
    # Check the extra arguments
    add_isos = False
    if not isinstance(pp.extra_args, type(None)):
        extra_args = pp.extra_args
        if 'isopycnals' in extra_args.keys():
            isopycnals = extra_args['isopycnals']
            if not isinstance(isopycnals, type(None)) and not isopycnals is False:
                add_isos = True
        if 'place_isos' in extra_args.keys():
            place_isos = extra_args['place_isos']
        else:
            place_isos = False
    else:
        extra_args = False
    # Concatonate all the pandas data frames together
    df = pd.concat(a_group.data_frames)
    # Set variable as to whether to invert the y axis
    invert_y_axis = False
    # Decide on a plot type
    if plot_type == 'xy':
        ## Make a standard x vs. y scatter plot
        # Set the main x and y data keys
        x_key = pp.x_vars[0]
        y_key = pp.y_vars[0]
        # Check for histogram
        if x_key == 'hist' or y_key == 'hist':
            plot_hist = True
        else:
            plot_hist = False
        # Check for twin x and y data keys
        try:
            tw_x_key = pp.x_vars[1]
            tw_ax_y  = ax.twiny()
        except:
            tw_x_key = None
            tw_ax_y  = None
        try:
            tw_y_key = pp.y_vars[1]
            tw_ax_x  = ax.twinx()
        except:
            tw_y_key = None
            tw_ax_x  = None
        # Check whether plotting a parameter sweep of clustering
        for var in [x_key, y_key, tw_x_key, tw_y_key]:
            if var in clstr_ps_vars:
                # Add a standard title
                plt_title = add_std_title(a_group)
                # Plot the parameter sweep
                xlabel, ylabel = plot_clstr_param_sweep(ax, tw_ax_x, a_group, plt_title)
                return xlabel, ylabel, plt_title, ax, invert_y_axis
        # Determine the color mapping to be used
        if clr_map in a_group.vars_to_keep and clr_map != 'cluster':
            if pp.plot_scale == 'by_pf':
                if clr_map not in pf_vars:
                    print('Cannot use',clr_map,'with plot scale by_pf')
                    exit(0)
            # Concatonate all the pandas data frames together
            df = pd.concat(a_group.data_frames)
            # Format the dates if necessary
            if x_key in ['dt_start', 'dt_end']:
                df[x_key] = mpl.dates.date2num(df[x_key])
            if y_key in ['dt_start', 'dt_end']:
                df[y_key] = mpl.dates.date2num(df[y_key])
            if clr_map == 'dt_start' or clr_map == 'dt_end':
                cmap_data = mpl.dates.date2num(df[clr_map])
            else:
                cmap_data = df[clr_map]
            # Get the colormap
            this_cmap = get_color_map(clr_map)
            #
            if plot_hist:
                return plot_histogram(x_key, y_key, ax, a_group, pp, clr_map=clr_map, legend=pp.legend)
            # If there aren't that many points, make the markers bigger
            if len(df[x_key]) < 1000:
                m_size = map_mrk_size
            else: 
                m_size = mrk_size
            heatmap = ax.scatter(df[x_key], df[y_key], c=cmap_data, cmap=this_cmap, s=m_size, marker=std_marker, zorder=5)
            # Invert y-axis if specified
            if y_key in y_invert_vars:
                invert_y_axis = True
            # Plot on twin axes, if specified
            if not isinstance(tw_x_key, type(None)):
                tw_clr = get_var_color(tw_x_key)
                if tw_clr == std_clr:
                    tw_clr = alt_std_clr
                # Add backing x to distinguish from main axis
                tw_ax_y.scatter(df[tw_x_key], df[y_key], color=tw_clr, s=m_size*10, marker='x', zorder=3)
                tw_ax_y.scatter(df[tw_x_key], df[y_key], c=cmap_data, cmap=this_cmap, s=m_size, marker=std_marker, zorder=4)
                tw_ax_y.set_xlabel(pp.xlabels[1])
                # Change color of the axis label on the twin axis
                tw_ax_y.xaxis.label.set_color(tw_clr)
                # Change color of the ticks on the twin axis
                tw_ax_y.tick_params(axis='x', colors=tw_clr)
                # Create the colorbar
                cbar = plt.colorbar(heatmap, ax=tw_ax_y)
            elif not isinstance(tw_y_key, type(None)):
                tw_clr = get_var_color(tw_y_key)
                if tw_clr == std_clr:
                    tw_clr = alt_std_clr
                # Add backing x to distinguish from main axis
                tw_ax_x.scatter(df[x_key], df[tw_y_key], color=tw_clr, s=m_size*10, marker='x', zorder=3)
                tw_ax_x.scatter(df[x_key], df[tw_y_key], c=cmap_data, cmap=this_cmap, s=m_size, marker=std_marker, zorder=4)
                # Invert y-axis if specified
                if tw_y_key in y_invert_vars:
                    tw_ax_x.invert_yaxis()
                tw_ax_x.set_ylabel(pp.ylabels[1])
                # Change color of the axis label on the twin axis
                tw_ax_x.yaxis.label.set_color(tw_clr)
                # Change color of the ticks on the twin axis
                tw_ax_x.tick_params(axis='y', colors=tw_clr)
                # Create the colorbar
                cbar = plt.colorbar(heatmap, ax=tw_ax_x)
            else:
                # Create the colorbar
                cbar = plt.colorbar(heatmap, ax=ax)
            # Format the colorbar ticks, if necessary
            if clr_map == 'dt_start' or clr_map == 'dt_end':
                loc = mpl.dates.AutoDateLocator()
                cbar.ax.yaxis.set_major_locator(loc)
                cbar.ax.yaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))
            cbar.set_label(pp.clabel)
            # Add a standard legend
            if pp.legend:
                add_std_legend(ax, df, x_key)
            # Format the axes for datetimes, if necessary
            format_datetime_axes(x_key, y_key, ax, tw_x_key, tw_ax_y, tw_y_key, tw_ax_x)
            # Check whether to plot isopycnals
            if add_isos:
                add_isopycnals(ax, df, x_key, y_key, p_ref=isopycnals, place_isos=place_isos, tw_x_key=tw_x_key, tw_ax_y=tw_ax_y, tw_y_key=tw_y_key, tw_ax_x=tw_ax_x)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
        elif clr_map == 'clr_all_same':
            # Concatonate all the pandas data frames together
            df = pd.concat(a_group.data_frames)
            # Check for cluster-based variables
            if x_key in clstr_vars or y_key in clstr_vars:
                m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt = get_cluster_args(pp)
                df, rel_val = HDBSCAN_(a_group, df, cl_x_var, cl_y_var, m_pts, min_samp=min_s, extra_cl_vars=[x_key,y_key, 'nir_SP'])
            else:
                # Check whether to plot slopes
                try:
                    plot_slopes = pp.extra_args['plot_slopes']
                except:
                    plot_slopes = False
            # Check for histogram
            if plot_hist:
                return plot_histogram(x_key, y_key, ax, a_group, pp, clr_map, legend=pp.legend, df=df, txk=tw_x_key, tay=tw_ax_y, tyk=tw_y_key, tax=tw_ax_x)
            # Format the dates if necessary
            if x_key in ['dt_start', 'dt_end']:
                df[x_key] = mpl.dates.date2num(df[x_key])
            if y_key in ['dt_start', 'dt_end']:
                df[y_key] = mpl.dates.date2num(df[y_key])
            # Drop duplicates
            df.drop_duplicates(subset=[x_key, y_key], keep='first', inplace=True)
            # If there aren't that many points, make the markers bigger
            if len(df[x_key]) < 1000:
                m_size = map_mrk_size
            else: 
                m_size = mrk_size
            # Plot every point the same color, size, and marker
            ax.scatter(df[x_key], df[y_key], color=std_clr, s=m_size, marker=std_marker, alpha=mrk_alpha, zorder=5)
            if plot_slopes:
                # Mark outliers
                mark_outliers(ax, df, x_key, y_key)
                # Get data without outliers
                x_data = np.array(df[df['out_'+x_key]==False][x_key].values, dtype=np.float64)
                y_data = np.array(df[df['out_'+x_key]==False][y_key].values, dtype=np.float64)
                # Find the slope of the ordinary least-squares of the points for this cluster
                # m, c = np.linalg.lstsq(np.array([x_data, np.ones(len(x_data))]).T, y_data, rcond=None)[0]
                # Find the slope of the total least-squares of the points for this cluster
                m, c, sd_m, sd_c = orthoregress(x_data, y_data)
                # Find mean and standard deviation of x and y data
                x_mean = df.loc[:,x_key].mean()
                y_mean = df.loc[:,y_key].mean()
                x_stdv = df.loc[:,x_key].std()
                y_stdv = df.loc[:,y_key].std()
                # Plot the least-squares fit line for this cluster through the centroid
                ax.axline((x_mean, y_mean), slope=m, color=alt_std_clr, zorder=3)
                # Add annotation to say what the slope is
                ax.annotate(r'%.2f$\pm$%.2f'%(m,sd_m), xy=(x_mean+x_stdv/4,y_mean+y_stdv/10), xycoords='data', color=alt_std_clr, weight='bold', zorder=12)
            # Invert y-axis if specified
            if y_key in y_invert_vars:
                invert_y_axis = True
            # Plot on twin axes, if specified
            if not isinstance(tw_x_key, type(None)):
                tw_clr = get_var_color(tw_x_key)
                if tw_clr == std_clr:
                    tw_clr = alt_std_clr
                tw_ax_y.scatter(df[tw_x_key], df[y_key], color=tw_clr, s=m_size, marker=std_marker, alpha=mrk_alpha)
                tw_ax_y.set_xlabel(pp.xlabels[1])
                # Change color of the ticks on the twin axis
                tw_ax_y.tick_params(axis='x', colors=tw_clr)
            elif not isinstance(tw_y_key, type(None)):
                tw_clr = get_var_color(tw_y_key)
                if tw_clr == std_clr:
                    tw_clr = alt_std_clr
                tw_ax_x.scatter(df[x_key], df[tw_y_key], color=tw_clr, s=m_size, marker=std_marker, alpha=mrk_alpha)
                # Invert y-axis if specified
                if tw_y_key in y_invert_vars:
                    tw_ax_x.invert_yaxis()
                tw_ax_x.set_ylabel(pp.ylabels[1])
                # Change color of the ticks on the twin axis
                tw_ax_x.tick_params(axis='y', colors=tw_clr)
            # Add a standard legend
            if pp.legend:
                add_std_legend(ax, df, x_key)
            # Format the axes for datetimes, if necessary
            format_datetime_axes(x_key, y_key, ax, tw_x_key, tw_ax_y, tw_y_key, tw_ax_x)
            # Check whether to plot isopycnals
            if add_isos:
                add_isopycnals(ax, df, x_key, y_key, p_ref=isopycnals, place_isos=place_isos, tw_x_key=tw_x_key, tw_ax_y=tw_ax_y, tw_y_key=tw_y_key, tw_ax_x=tw_ax_x)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
        elif clr_map == 'clr_by_source':
            if plot_hist:
                return plot_histogram(x_key, y_key, ax, a_group, pp, clr_map, legend=pp.legend)
            # Find the list of sources
            sources_list = []
            for df in a_group.data_frames:
                # Get unique sources
                these_sources = np.unique(df['source'])
                for s in these_sources:
                    sources_list.append(s)
            # The pandas version of 'unique()' preserves the original order
            sources_list = pd.unique(pd.Series(sources_list))
            # If there aren't that many points, make the markers bigger
            if len(df[x_key]) < 1000:
                m_size = map_mrk_size
            else: 
                m_size = mrk_size
            i = 0
            lgnd_hndls = []
            for source in sources_list:
                # Decide on the color, don't go off the end of the array
                my_clr = distinct_clrs[i%len(distinct_clrs)]
                these_dfs = []
                for df in a_group.data_frames:
                    these_dfs.append(df[df['source'] == source])
                this_df = pd.concat(these_dfs)
                # Format the dates if necessary
                if x_key in ['dt_start', 'dt_end']:
                    this_df[x_key] = mpl.dates.date2num(this_df[x_key])
                if y_key in ['dt_start', 'dt_end']:
                    this_df[y_key] = mpl.dates.date2num(this_df[y_key])
                # Plot every point from this df the same color, size, and marker
                ax.scatter(this_df[x_key], this_df[y_key], color=my_clr, s=m_size, marker=std_marker, alpha=mrk_alpha, zorder=5)
                i += 1
                # Add legend to report the total number of points for this instrmt
                lgnd_label = source+': '+str(len(this_df[x_key]))+' points'
                lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label))
            # Invert y-axis if specified
            if y_key in y_invert_vars:
                invert_y_axis = True
            notes_string = ''.join(this_df.notes.unique())
            # Only add the notes_string if it contains something
            if len(notes_string) > 1:
                notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
                lgnd_hndls.append(notes_patch)
            # Format the axes for datetimes, if necessary
            format_datetime_axes(x_key, y_key, ax)
            # Add legend with custom handles
            if pp.legend:
                lgnd = ax.legend(handles=lgnd_hndls)
            # Check whether to plot isopycnals
            if add_isos:
                add_isopycnals(ax, df, x_key, y_key, p_ref=isopycnals, place_isos=place_isos, tw_x_key=tw_x_key, tw_ax_y=tw_ax_y, tw_y_key=tw_y_key, tw_ax_x=tw_ax_x)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
        elif clr_map == 'clr_by_instrmt':
            if plot_hist:
                return plot_histogram(x_key, y_key, ax, a_group, pp, clr_map, legend=pp.legend)
            # If there aren't that many points, make the markers bigger
            if len(df[x_key]) < 1000:
                m_size = map_mrk_size
            else: 
                m_size = mrk_size
            i = 0
            lgnd_hndls = []
            # Loop through each data frame, the same as looping through instrmts
            for df in a_group.data_frames:
                # Format the dates if necessary
                if x_key in ['dt_start', 'dt_end']:
                    df[x_key] = mpl.dates.date2num(df[x_key])
                if y_key in ['dt_start', 'dt_end']:
                    df[y_key] = mpl.dates.date2num(df[y_key])
                df['source-instrmt'] = df['source']+' '+df['instrmt']
                # Get instrument name
                s_instrmt = np.unique(df['source-instrmt'])[0]
                # Decide on the color, don't go off the end of the array
                my_clr = distinct_clrs[i%len(distinct_clrs)]
                # Plot every point from this df the same color, size, and marker
                ax.scatter(df[x_key], df[y_key], color=my_clr, s=m_size, marker=std_marker, alpha=mrk_alpha, zorder=5)
                i += 1
                # Add legend to report the total number of points for this instrmt
                lgnd_label = s_instrmt+': '+str(len(df[x_key]))+' points'
                lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label))
                notes_string = ''.join(df.notes.unique())
            # Invert y-axis if specified
            if y_key in y_invert_vars:
                invert_y_axis = True
            # Only add the notes_string if it contains something
            if len(notes_string) > 1:
                notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
                lgnd_hndls.append(notes_patch)
            # Add legend with custom handles
            if pp.legend:
                lgnd = ax.legend(handles=lgnd_hndls)
            # Format the axes for datetimes, if necessary
            format_datetime_axes(x_key, y_key, ax)
            # Check whether to plot isopycnals
            if add_isos:
                add_isopycnals(ax, df, x_key, y_key, p_ref=isopycnals, place_isos=place_isos, tw_x_key=tw_x_key, tw_ax_y=tw_ax_y, tw_y_key=tw_y_key, tw_ax_x=tw_ax_x)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
        elif clr_map == 'density_hist':
            if plot_hist:
                print('Cannot use density_hist colormap with the `hist` variable')
                exit(0)
            # Concatonate all the pandas data frames together
            df = pd.concat(a_group.data_frames)
            # Check for cluster-based variables
            if x_key in clstr_vars or y_key in clstr_vars:
                m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt = get_cluster_args(pp)
                df, rel_val = HDBSCAN_(a_group, df, cl_x_var, cl_y_var, m_pts, min_samp=min_s, extra_cl_vars=[x_key,y_key])
            # Format the dates if necessary
            if x_key in ['dt_start', 'dt_end']:
                df[x_key] = mpl.dates.date2num(df[x_key])
            if y_key in ['dt_start', 'dt_end']:
                df[y_key] = mpl.dates.date2num(df[y_key])
            # Plot a density histogram where each grid box is colored to show how
            #   many points fall within it
            # Get the density histogram plotting parameters from extra args
            d_hist_dict = pp.extra_args
            if isinstance(d_hist_dict, type(None)):
                # If none are given, use the defaults here
                clr_min = 0
                clr_max = 20
                clr_ext = 'max'        # adds arrow indicating values go past the bounds
                #                       #   use 'min', 'max', or 'both'
                xy_bins = 250
            else:
                # Pull out the values from the dictionary provided
                clr_min = d_hist_dict['clr_min']
                clr_max = d_hist_dict['clr_max']
                clr_ext = d_hist_dict['clr_ext']
                xy_bins = d_hist_dict['xy_bins']
            # Make the 2D histogram, the number of bins really changes the outcome
            heatmap = ax.hist2d(df[x_key], df[y_key], bins=xy_bins, cmap=get_color_map(clr_map), vmin=clr_min, vmax=clr_max)
            # Invert y-axis if specified
            if y_key in y_invert_vars:
                invert_y_axis = True
            # `hist2d` returns a tuple, the index 3 of which is the mappable for a colorbar
            cbar = plt.colorbar(heatmap[3], ax=ax, extend=clr_ext)
            cbar.set_label('points per pixel')
            # Add legend to report the total number of points and pixels on the plot
            n_pts_patch  = mpl.patches.Patch(color='none', label=str(len(df[x_key]))+' points')
            pixel_patch  = mpl.patches.Patch(color='none', label=str(xy_bins)+'x'+str(xy_bins)+' pixels')
            notes_string = ''.join(df.notes.unique())
            # Only add the notes_string if it contains something
            if pp.legend:
                if len(notes_string) > 1:
                    notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
                    ax.legend(handles=[n_pts_patch, pixel_patch, notes_patch])
                else:
                    ax.legend(handles=[n_pts_patch, pixel_patch])
            # Format the axes for datetimes, if necessary
            format_datetime_axes(x_key, y_key, ax)
            # Check whether to plot isopycnals
            if add_isos:
                add_isopycnals(ax, df, x_key, y_key, p_ref=isopycnals, place_isos=place_isos, tw_x_key=tw_x_key, tw_ax_y=tw_ax_y, tw_y_key=tw_y_key, tw_ax_x=tw_ax_x)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
        elif clr_map == 'cluster':
            # Add a standard title
            plt_title = add_std_title(a_group)
            #   Legend is handled inside plot_clusters()
            # Concatonate all the pandas data frames together
            df = pd.concat(a_group.data_frames)
            # Check for histogram
            if plot_hist:
                return plot_histogram(x_key, y_key, ax, a_group, pp, clr_map, legend=pp.legend, df=df, txk=tw_x_key, tay=tw_ax_y, tyk=tw_y_key, tax=tw_ax_x)
            # Format the dates if necessary
            if x_key in ['dt_start', 'dt_end']:
                df[x_key] = mpl.dates.date2num(df[x_key])
            if y_key in ['dt_start', 'dt_end']:
                df[y_key] = mpl.dates.date2num(df[y_key])
            m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt = get_cluster_args(pp)
            invert_y_axis = plot_clusters(a_group, ax, pp, df, x_key, y_key, cl_x_var, cl_y_var, clr_map, m_pts, min_samp=min_s, box_and_whisker=b_a_w_plt, plot_slopes=plot_slopes)
            # Format the axes for datetimes, if necessary
            format_datetime_axes(x_key, y_key, ax)
            # Check whether to plot isopycnals
            if add_isos:
                add_isopycnals(ax, df, x_key, y_key, p_ref=isopycnals, place_isos=place_isos, tw_x_key=tw_x_key, tw_ax_y=tw_ax_y, tw_y_key=tw_y_key, tw_ax_x=tw_ax_x)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
            #
        else:
            # Did not provide a valid clr_map
            print('Colormap',clr_map,'not valid')
            exit(0)
        #
    #
    elif plot_type == 'map':
        ## Plot profile locations on a map of the Arctic Ocean
        # Get the map plotting parameters from extra args
        map_dict = pp.extra_args
        # See what was in that dictionary
        try:
            map_extent = map_dict['map_extent']
        except:
            map_extent = None
        #   Set latitude and longitude extents
        #   I found these values by guess-and-check, there really isn't a good way
        #       to know beforehand what you'll actually get
        if map_extent == 'Canada_Basin':
            cent_lon = -135 #-140
            ex_N = 78 #80
            ex_S = 73.5 #69
            ex_E = -145 #-156
            ex_W = -125 #-124
        elif map_extent == 'Western_Arctic':
            cent_lon = -140
            ex_N = 80
            ex_S = 69
            ex_E = -165
            ex_W = -124
        else:
            cent_lon = -120
            ex_N = 90
            ex_S = 70
            ex_E = -180
            ex_W = 180
        # Remove the current axis
        ax.remove()
        # Replace axis with one that can be made into a map
        ax = fig.add_subplot(ax_pos, projection=ccrs.NorthPolarStereo(central_longitude=cent_lon))
        ax.set_extent([ex_E, ex_W, ex_S, ex_N], ccrs.PlateCarree())
        #   Add ocean first, then land. Otherwise the ocean covers the land shapes
        ax.add_feature(cartopy.feature.OCEAN, color=bathy_clrs[0])
        ax.add_feature(cartopy.feature.LAND, color=clr_land, alpha=0.5)
        # Make bathymetry features
        bathy_0200 = cartopy.feature.NaturalEarthFeature(category='physical',name='bathymetry_K_200',scale='10m')
        bathy_1000 = cartopy.feature.NaturalEarthFeature(category='physical',name='bathymetry_J_1000',scale='10m')
        bathy_2000 = cartopy.feature.NaturalEarthFeature(category='physical',name='bathymetry_I_2000',scale='10m')
        bathy_3000 = cartopy.feature.NaturalEarthFeature(category='physical',name='bathymetry_H_3000',scale='10m')
        bathy_4000 = cartopy.feature.NaturalEarthFeature(category='physical',name='bathymetry_G_4000',scale='10m')
        bathy_5000 = cartopy.feature.NaturalEarthFeature(category='physical',name='bathymetry_F_5000',scale='10m')
        # Add bathymetry lines
        ax.add_feature(bathy_0200, facecolor=bathy_clrs[1], zorder=1)
        ax.add_feature(bathy_1000, facecolor=bathy_clrs[2], zorder=2)
        ax.add_feature(bathy_2000, facecolor=bathy_clrs[3], zorder=3)
        ax.add_feature(bathy_3000, facecolor=bathy_clrs[4], zorder=4)
        # ax.add_feature(bathy_4000, facecolor=bathy_clrs[5], zorder=5)
        # ax.add_feature(bathy_5000, facecolor=bathy_clrs[6], zorder=6)
        #   Add gridlines to show longitude and latitude
        gl = ax.gridlines(draw_labels=True, color=clr_lines, alpha=grid_alpha, linestyle='--', zorder=7)
        #       x is actually all labels around the edge
        gl.xlabel_style = {'size':8, 'color':clr_lines, 'zorder':7}
        #       y is actually all labels within map
        gl.ylabel_style = {'size':8, 'color':clr_lines, 'zorder':7}
        #   Plotting the coastlines takes a really long time
        # ax.coastlines()
        # Add bounding box for Canada Basin
        #   Don't plot outside the extent chosen
        if map_extent == 'Canada_Basin':
            # Only the Eastern boundary appears in this extent
            # ax.plot([-130,-130], [73.7, 78.15], color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Eastern boundary
            # Add bathymetry lines
            ax.add_feature(bathy_0200, facecolor='none', edgecolor=std_clr, linestyle='-.', alpha=0.3, zorder=1)
            ax.add_feature(bathy_1000, facecolor='none', edgecolor=std_clr, linestyle='-.', alpha=0.3, zorder=2)
            ax.add_feature(bathy_2000, facecolor='none', edgecolor=std_clr, linestyle='-.', alpha=0.3, zorder=3)
            ax.add_feature(bathy_3000, facecolor='none', edgecolor=std_clr, linestyle='-.', alpha=0.3, zorder=4)
            # ax.add_feature(bathy_4000, facecolor='none', edgecolor=std_clr, linestyle='-.', alpha=0.3, zorder=5)
            # ax.add_feature(bathy_5000, facecolor='none', edgecolor=std_clr, linestyle='-.', alpha=0.3, zorder=6)
        elif map_extent == 'Western_Arctic':
            CB_lons = np.linspace(-130, -155, 50)
            # Northern boundary does not appear
            ax.plot(CB_lons, 72*np.ones(len(CB_lons)), color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Southern boundary
            ax.plot([-130,-130], [72, 80.5], color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Eastern boundary
            ax.plot([-155,-155], [72, 80.5], color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Western boundary
        else:
            CB_lons = np.linspace(-130, -155, 50)
            ax.plot(CB_lons, 84*np.ones(len(CB_lons)), color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Northern boundary
            ax.plot(CB_lons, 72*np.ones(len(CB_lons)), color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Southern boundary
            ax.plot([-130,-130], [72, 84], color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Eastern boundary
            ax.plot([-155,-155], [72, 84], color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Western boundary# Only the Eastern boundary appears in this extent
            ax.plot([-130,-130], [73.7, 78.15], color='red', linewidth=1, linestyle='-', transform=ccrs.Geodetic(), zorder=8) # Eastern boundary
        # Determine the color mapping to be used
        if clr_map in a_group.vars_to_keep:
            # Make sure it isn't a vertical variable
            if clr_map in vertical_vars:
                print('Cannot use',clr_map,'as colormap for a map plot')
                exit(0)
            # Concatonate all the pandas data frames together
            df = pd.concat(a_group.data_frames)
            # If not plotting very many points, increase the marker size
            if len(df) < 10:
                mrk_s = big_map_mrkr
            else:
                mrk_s = map_mrk_size
            # Format the dates if necessary
            if clr_map == 'dt_start' or clr_map == 'dt_end':
                cmap_data = mpl.dates.date2num(df[clr_map])
            else:
                cmap_data = df[clr_map]
            # The color of each point corresponds to the number of the profile it came from
            heatmap = ax.scatter(df['lon'], df['lat'], c=cmap_data, cmap=get_color_map(clr_map), s=mrk_s, marker=map_marker, linewidths=map_ln_wid, transform=ccrs.PlateCarree(), zorder=10)
            # Create the colorbar
            cbar = plt.colorbar(heatmap, ax=ax)
            # Format the colorbar ticks, if necessary
            if clr_map == 'dt_start' or clr_map == 'dt_end':
                loc = mpl.dates.AutoDateLocator()
                cbar.ax.yaxis.set_major_locator(loc)
                cbar.ax.yaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))
            cbar.set_label(a_group.data_set.var_attr_dicts[0][clr_map]['label'])
            # Add a standard legend
            add_std_legend(ax, df, 'lon')
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, False
        if clr_map == 'clr_all_same':
            # Concatonate all the pandas data frames together
            df = pd.concat(a_group.data_frames)
            # If not plotting very many points, increase the marker size
            if len(df) < 10:
                mrk_s = big_map_mrkr
            else:
                mrk_s = map_mrk_size
            # Plot every point the same color, size, and marker
            ax.scatter(df['lon'], df['lat'], color=std_clr, s=mrk_s, marker=map_marker, alpha=mrk_alpha, linewidths=map_ln_wid, transform=ccrs.PlateCarree(), zorder=10)
            # Add a standard legend
            add_std_legend(ax, df, 'lon')
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, False
        elif clr_map == 'clr_by_source':
            # Find the list of sources
            sources_list = []
            for df in a_group.data_frames:
                # Get unique sources
                these_sources = np.unique(df['source'])
                for s in these_sources:
                    sources_list.append(s)
            # The pandas version of 'unique()' preserves the original order
            sources_list = pd.unique(pd.Series(sources_list))
            i = 0
            lgnd_hndls = []
            for source in sources_list:
                # Decide on the color, don't go off the end of the array
                my_clr = distinct_clrs[i%len(distinct_clrs)]
                these_dfs = []
                for df in a_group.data_frames:
                    these_dfs.append(df[df['source'] == source])
                this_df = pd.concat(these_dfs)
                # If not plotting very many points, increase the marker size
                if len(this_df) < 10:
                    mrk_s = big_map_mrkr
                else:
                    mrk_s = map_mrk_size
                # Plot every point from this df the same color, size, and marker
                ax.scatter(this_df['lon'], this_df['lat'], color=my_clr, s=mrk_s, marker=map_marker, alpha=mrk_alpha, linewidths=map_ln_wid, transform=ccrs.PlateCarree(), zorder=10)
                i += 1
                # Add legend to report the total number of points for this instrmt
                lgnd_label = source+': '+str(len(this_df['lon']))+' points'
                lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label))
                notes_string = ''.join(this_df.notes.unique())
            # Only add the notes_string if it contains something
            if len(notes_string) > 1:
                notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
                lgnd_hndls.append(notes_patch)
            # Add legend with custom handles
            lgnd = ax.legend(handles=lgnd_hndls)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, False
        elif clr_map == 'clr_by_instrmt':
            i = 0
            lgnd_hndls = []
            # Loop through each data frame, the same as looping through instrmts
            for df in a_group.data_frames:
                df['source-instrmt'] = df['source']+' '+df['instrmt']
                # If not plotting very many points, increase the marker size
                if len(df) < 10:
                    if map_extent == 'Canada_Basin':
                        mrk_s = big_map_mrkr
                    else:
                        mrk_s = map_mrk_size
                else:
                    if map_extent == 'Canada_Basin':
                        mrk_s = map_mrk_size
                    else:
                        mrk_s = sml_map_mrkr
                # Get instrument name
                s_instrmt = np.unique(df['source-instrmt'])[0]
                # Decide on the color, don't go off the end of the array
                my_clr = distinct_clrs[i%len(distinct_clrs)]
                # Plot every point from this df the same color, size, and marker
                ax.scatter(df['lon'], df['lat'], color=my_clr, s=mrk_s, marker=map_marker, alpha=map_alpha, linewidths=map_ln_wid, transform=ccrs.PlateCarree(), zorder=10)
                # Find first and last profile, based on entry
                entries = df['entry'].values
                first_pf_df = df[df['entry']==min(entries)]
                last_pf_df  = df[df['entry']==max(entries)]
                # Plot first and last profiles on map
                ax.scatter(first_pf_df['lon'], first_pf_df['lat'], color=my_clr, s=mrk_s*5, marker='>', alpha=map_alpha, linewidths=map_ln_wid, transform=ccrs.PlateCarree(), zorder=11)
                ax.scatter(last_pf_df['lon'], last_pf_df['lat'], color=my_clr, s=mrk_s*5, marker='s', alpha=map_alpha, linewidths=map_ln_wid, transform=ccrs.PlateCarree(), zorder=11)
                # Increase i to get next color in the array
                i += 1
                # Add legend to report the total number of points for this instrmt
                lgnd_label = s_instrmt+': '+str(len(df['lon']))+' profiles'
                lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label))
                notes_string = ''.join(df.notes.unique())
            # Only add the notes_string if it contains something
            if len(notes_string) > 1:
                notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
                lgnd_hndls.append(notes_patch)
            # Add legend entries for stopping and starting locations
            lgnd_hndls.append(mpl.lines.Line2D([],[],color=std_clr, label='Start', marker='>', linewidth=0))
            lgnd_hndls.append(mpl.lines.Line2D([],[],color=std_clr, label='Stop', marker='s', linewidth=0))
            # Add legend with custom handles
            # Add a standard legend
            if pp.legend:
                lgnd = ax.legend(handles=lgnd_hndls, loc='lower left')
            # Create the colorbar for the bathymetry
            if map_extent == 'Full_Arctic':
                cbar = plt.colorbar(bathy_smap, ax=ax, extend='min', shrink=0.7, pad=0.15)# fraction=0.05)
                # Fixing ticks
                ticks_loc = cbar.ax.get_yticks().tolist()
                cbar.ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(ticks_loc))
                cbar.ax.set_yticklabels(['','3000','2000','1000','200','0'])
                cbar.ax.tick_params(labelsize=font_size_ticks)
                cbar.set_label('Bathymetry (m)', size=font_size_labels)
            # Add a standard title
            plt_title = add_std_title(a_group)
            return pp.xlabels[0], pp.ylabels[0], plt_title, ax, False
    elif plot_type == 'profiles':
        # Check for clustering
        # Call profile plotting function
        return plot_profiles(ax, a_group, pp)
    else:
        # Did not provide a valid plot type
        print('Plot type',plot_type,'not valid')
        exit(0)
    #
    print('ERROR: Should not have gotten here')
    return a_group.xlabel, a_group.ylabel, plt_title, ax, False, False

################################################################################
# Auxiliary plotting functions #################################################
################################################################################

def add_isopycnals(ax, df, x_key, y_key, p_ref=None, place_isos=False, tw_x_key=None, tw_ax_y=None, tw_y_key=None, tw_ax_x=None):
    """
    Adds lines of constant density anomaly, if applicable

    ax          The main axis on which to format
    df          Pandas dataframe of the data these isopycnals will be plotted under
    x_key       The string of the name for the x data on the main axis
    y_key       The string of the name for the y data on the main axis
    p_ref       Value of pressure to reference the isopycnals to. Default is surface (0 dbar)
    place_isos  How to place the isopycnals, either 'auto' or 'manual'
    tw_x_key    The string of the name for the x data on the twin axis
    tw_ax_y     The twin y axis on which to format
    tw_y_key    The string of the name for the y data on the twin axis
    tw_ax_x     The twin x axis on which to format
    """
    # Check whether to make isopycnals or not
    intersect_arr = set([x_key, y_key, tw_x_key, tw_y_key]).intersection(['iT', 'CT', 'PT', 'SP', 'SA'])
    if len(intersect_arr) < 2:
        intersect_arr = set([x_key, y_key, tw_x_key, tw_y_key]).intersection(['aiT','aCT','aPT','BSP','BSt','BSA'])
        if len(intersect_arr) < 2:
            return
        else:
            ## Plot lines of slope -1 all across the domain
            # Get bounds of axes
            x_bnds = ax.get_xbound()
            y_bnds = ax.get_ybound()
            # print('x_bnds:',x_bnds)
            # print('y_bnds:',y_bnds)
            x_arr = np.linspace(x_bnds[0], x_bnds[1], 50)
            y_arr = [(y_bnds[1]-y_bnds[0])/2]*len(x_arr)
            # Plot the least-squares fit line for this cluster through the centroid
            for i in range(len(x_arr)):
                ax.axline((x_arr[i], y_arr[i]), slope=1, color=std_clr, alpha=0.25, linestyle='--', zorder=1)
            # Set offset on axis, incredibly specific for making a nice plot for the paper
            ax.ticklabel_format(axis='x', style='sci', scilimits=sci_lims, useMathText=True, useOffset=0.0268)
            return
    # Zoom in on data so there's no border space
    x_var_s = df.loc[:,x_key]
    x_var_min = x_var_s.min()
    x_var_max = x_var_s.max()
    ax.set_xlim([x_var_min, x_var_max])
    y_var_s = df.loc[:,y_key]
    y_var_min = y_var_s.min()
    y_var_max = y_var_s.max()
    ax.set_ylim([y_var_min, y_var_max])
    # Get bounds of axes
    x_bnds = ax.get_xbound()
    y_bnds = ax.get_ybound()
    # print('x_bnds:',x_bnds)
    # print('y_bnds:',y_bnds)
    # Figure out which orientation to make the grid
    if x_key in ['iT','CT','PT']:
        iso_x_key = 'iT'
    elif x_key in ['SP','SA']:
        iso_x_key = 'SA'
    if y_key in ['iT','CT','PT']:
        iso_y_key = 'iT'
    elif y_key in ['SP','SA']:
        iso_y_key = 'SA'
    # If no reference pressure is given, take p_ref to be the median pressure
    if p_ref==True:
        press_s = df.loc[:, 'press']
        p_ref = press_s.median()
    print('\t- Adding isopycnals referenced to:',p_ref,'dbar')
    # Get bounds of iso keys that fit within axes bounds
    iso_x_min = df[(df[x_key] == x_var_min)].loc[:,iso_x_key].min()
    iso_x_max = df[(df[x_key] == x_var_max)].loc[:,iso_x_key].max()
    iso_y_min = df[(df[y_key] == y_var_min)].loc[:,iso_y_key].min()
    iso_y_max = df[(df[y_key] == y_var_max)].loc[:,iso_y_key].max()
    # Number of points for the mesh grid in both directions
    n_grid_pts = 100
    # Make meshgrid for calculating
    iso_x_arr = np.arange(iso_x_min, iso_x_max, abs(iso_x_max-iso_x_min)/n_grid_pts)
    iso_y_arr = np.arange(iso_y_min, iso_y_max, abs(iso_y_max-iso_y_min)/n_grid_pts)
    iso_X, iso_Y = np.meshgrid(iso_x_arr, iso_y_arr)
    # Make grid of pressure values
    P = iso_X*0+p_ref
    # Calculate potential density referenced to p_ref
    if iso_x_key == 'SA' and iso_y_key == 'iT':
        Z = gsw.pot_rho_t_exact(iso_X, iso_Y, P, p_ref)-1000
    # Make meshgrid for plotting
    x_arr = np.arange(x_var_min, x_var_max, abs(x_var_max-x_var_min)/n_grid_pts)
    y_arr = np.arange(y_var_min, y_var_max, abs(y_var_max-y_var_min)/n_grid_pts)
    X, Y = np.meshgrid(x_arr, y_arr)
    # Plot the contours
    CS = ax.contour(X, Y, Z, colors=std_clr, alpha=0.4, zorder=1, linestyles='dashed')
    # Place contour inline labels
    if place_isos == 'manual':
        # Place contour labels manualy, interactively
        ax.clabel(CS, manual=True, fontsize=10, fmt='%.1f')
    else:
        # Place contour labels automatically
        ax.clabel(CS, inline=True, fontsize=10, fmt='%.1f')

################################################################################

def plot_histogram(x_key, y_key, ax, a_group, pp, clr_map, legend=True, df=None, txk=None, tay=None, tyk=None, tax=None):
    """
    Takes in an Analysis_Group object which has the data and plotting parameters
    to produce a subplot of individual profiles. Returns the x and y labels and
    the subplot title

    x_key       The string of the name for the x data on the main axis
    y_key       The string of the name for the y data on the main axis
    ax          The axis on which to make the plot
    a_group     A Analysis_Group object containing the info to create this subplot
    pp          The Plot_Parameters object for a_group
    clr_map     A string to determine what color map to use in the plot
    legend      True/False whether to add a legend
    df          A pandas dataframe
    txk         Twin x variable key
    tay         Twin y axis
    tyk         Twin y variable key
    tax         Twin x axis
    """
    # Find the histogram parameters, if given
    if not isinstance(pp.extra_args, type(None)):
        try:
            n_h_bins = pp.extra_args['n_h_bins']
        except:
            n_h_bins = None
        try:
            plt_hist_lines = pp.extra_args['plt_hist_lines']
        except:
            plt_hist_lines = False
        try:
            plt_noise = pp.extra_args['plt_noise']
        except:
            plt_noise = True
    else:
        n_h_bins = None
        plt_hist_lines = False
        plt_noise = True
    # Load in variables
    if x_key == 'hist':
        var_key = y_key
        orientation = 'horizontal'
        x_label = 'Number of points'
        y_label = pp.ylabels[0]
        tw_var_key = tyk
        tw_ax = tax
        tw_label = pp.ylabels[1]
    elif y_key == 'hist':
        var_key = x_key
        orientation = 'vertical'
        x_label = pp.xlabels[0]
        y_label = 'Number of points'
        tw_var_key = txk
        tw_ax = tay
        tw_label = pp.xlabels[1]
    # Histograms don't work with datetime data
    if var_key in ['dt_start', 'dt_end']:
        print('Cannot make histogram with',var_key,'data')
        exit(0)
    # Set variable as to whether to invert the y axis
    invert_y_axis = False
    # Determine the color mapping to be used
    if clr_map == 'clr_all_same':
        # Get histogram parameters
        h_var, res_bins, median, mean, std_dev = get_hist_params(df, var_key, n_h_bins)
        # Plot the histogram
        ax.hist(h_var, bins=res_bins, color=std_clr, orientation=orientation)
        # Invert y-axis if specified
        if y_key in y_invert_vars:
            invert_y_axis = True
        # Check whether to plot lines for mean and standard deviation
        if plt_hist_lines:
            if orientation == 'vertical':
                ax.axvline(mean, color='r')
                ax.axvline(mean-2*std_dev, color='r', linestyle='--')
                ax.axvline(mean+2*std_dev, color='r', linestyle='--')
            elif orientation == 'horizontal':
                ax.axhline(mean, color='r')
                ax.axhline(mean-2*std_dev, color='r', linestyle='--')
                ax.axhline(mean+2*std_dev, color='r', linestyle='--')
        # Add legend to report overall statistics
        n_pts_patch   = mpl.patches.Patch(color=std_clr, label=str(len(h_var))+' points')
        median_patch  = mpl.patches.Patch(color=std_clr, label='Median:  '+'%.4f'%median)
        mean_patch    = mpl.patches.Patch(color=std_clr, label='Mean:    ' + '%.4f'%mean)
        std_dev_patch = mpl.patches.Patch(color=std_clr, label='Std dev: '+'%.4f'%std_dev)
        hndls = [n_pts_patch, median_patch, mean_patch, std_dev_patch]
        notes_string = ''.join(df.notes.unique())
        # Plot twin axis if specified
        if not isinstance(tw_var_key, type(None)):
            tw_clr = get_var_color(tw_var_key)
            if tw_clr == std_clr:
                tw_clr = alt_std_clr
            # Get histogram parameters
            h_var, res_bins, median, mean, std_dev = get_hist_params(df, tw_var_key)
            # Plot the histogram
            tw_ax.hist(h_var, bins=res_bins, color=tw_clr, alpha=mrk_alpha, orientation=orientation)
            if orientation == 'vertical':
                tw_ax.set_xlabel(tw_label)
                tw_ax.tick_params(axis='x', colors=tw_clr)
            elif orientation == 'horizontal':
                tw_ax.set_ylabel(tw_label)
                tw_ax.tick_params(axis='y', colors=tw_clr)
            # Check whether to plot lines for mean and standard deviation
            if plt_hist_lines:
                if orientation == 'vertical':
                    tw_ax.axvline(mean, color='r', alpha=0.5)
                    tw_ax.axvline(mean-2*std_dev, color='r', linestyle='--', alpha=0.5)
                    tw_ax.axvline(mean+2*std_dev, color='r', linestyle='--', alpha=0.5)
                elif orientation == 'horizontal':
                    tw_ax.axhline(mean, color='r', alpha=0.5)
                    tw_ax.axhline(mean-2*std_dev, color='r', linestyle='--', alpha=0.5)
                    tw_ax.axhline(mean+2*std_dev, color='r', linestyle='--', alpha=0.5)
            # Add legend to report overall statistics
            n_pts_patch   = mpl.patches.Patch(color=tw_clr, label=str(len(h_var))+' points', alpha=mrk_alpha)
            median_patch  = mpl.patches.Patch(color=tw_clr, label='Median:  '+'%.4f'%median, alpha=mrk_alpha)
            mean_patch    = mpl.patches.Patch(color=tw_clr, label='Mean:    ' + '%.4f'%mean, alpha=mrk_alpha)
            std_dev_patch = mpl.patches.Patch(color=tw_clr, label='Std dev: '+'%.4f'%std_dev, alpha=mrk_alpha)
            tw_notes_string = ''.join(df.notes.unique())
            # Only add the notes_string if it contains something
            if len(tw_notes_string) > 1:
                tw_notes_patch  = mpl.patches.Patch(color='none', label=notes_string, alpha=mrk_alpha)
                tw_hndls=[n_pts_patch, median_patch, mean_patch, std_dev_patch, tw_notes_patch]
            else:
                tw_hndls=[n_pts_patch, median_patch, mean_patch, std_dev_patch]
        else:
            tw_hndls = []
        # Only add the notes_string if it contains something
        if len(notes_string) > 1:
            notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
            ax.legend(handles=hndls+notes_patch+tw_hndls)
        else:
            ax.legend(handles=hndls+tw_hndls)
        # Add a standard title
        plt_title = add_std_title(a_group)
        return x_label, y_label, plt_title, ax, invert_y_axis
    elif clr_map == 'clr_by_source':
        # Can't make this type of plot for certain variables yet
        if var_key in clstr_vars:
            print('Cannot yet make a histogram plot of',var_key)
        # Find the list of sources
        sources_list = []
        for df in a_group.data_frames:
            # Get unique sources
            these_sources = np.unique(df['source'])
            for s in these_sources:
                sources_list.append(s)
        # The pandas version of 'unique()' preserves the original order
        sources_list = pd.unique(pd.Series(sources_list))
        i = 0
        lgnd_hndls = []
        for source in sources_list:
            # Decide on the color, don't go off the end of the array
            my_clr = distinct_clrs[i%len(distinct_clrs)]
            these_dfs = []
            for df in a_group.data_frames:
                these_dfs.append(df[df['source'] == source])
            this_df = pd.concat(these_dfs)
            # Get histogram parameters
            h_var, res_bins, median, mean, std_dev = get_hist_params(this_df, var_key, n_h_bins)
            # Plot the histogram
            ax.hist(h_var, bins=res_bins, color=my_clr, alpha=mrk_alpha, orientation=orientation)
            # Check whether to plot lines for mean and standard deviation
            if plt_hist_lines:
                if orientation == 'vertical':
                    ax.axvline(mean, color='r')
                    ax.axvline(mean-2*std_dev, color='r', linestyle='--')
                    ax.axvline(mean+2*std_dev, color='r', linestyle='--')
                elif orientation == 'horizontal':
                    ax.axhline(mean, color='r')
                    ax.axhline(mean-2*std_dev, color='r', linestyle='--')
                    ax.axhline(mean+2*std_dev, color='r', linestyle='--')
            i += 1
            # Add legend handle to report the total number of points for this source
            lgnd_label = source+': '+str(len(this_df[var_key]))+' points, Median:'+'%.4f'%median
            lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label, alpha=mrk_alpha))
            # Add legend handle to report overall statistics
            lgnd_label = 'Mean:'+ '%.4f'%mean+', Std dev:'+'%.4f'%std_dev
            lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label, alpha=mrk_alpha))
            notes_string = ''.join(this_df.notes.unique())
        # Only add the notes_string if it contains something
        if len(notes_string) > 1:
            notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
            lgnd_hndls.append(notes_patch)
        # Add legend with custom handles
        lgnd = ax.legend(handles=lgnd_hndls)
        # Invert y-axis if specified
        if y_key in y_invert_vars:
            invert_y_axis = True
        # Add a standard title
        plt_title = add_std_title(a_group)
        return x_label, y_label, plt_title, ax, invert_y_axis
    elif clr_map == 'clr_by_instrmt':
        # Can't make this type of plot for certain variables yet
        if var_key in clstr_vars:
            print('Cannot yet make a histogram plot of',var_key)
        i = 0
        lgnd_hndls = []
        # Loop through each data frame, the same as looping through instrmts
        for df in a_group.data_frames:
            df['source-instrmt'] = df['source']+' '+df['instrmt']
            # Get instrument name
            s_instrmt = np.unique(df['source-instrmt'])[0]
            # Decide on the color, don't go off the end of the array
            my_clr = distinct_clrs[i%len(distinct_clrs)]
            # Get histogram parameters
            h_var, res_bins, median, mean, std_dev = get_hist_params(df, var_key, n_h_bins)
            # Plot the histogram
            ax.hist(h_var, bins=res_bins, color=my_clr, alpha=mrk_alpha, orientation=orientation)
            # Check whether to plot lines for mean and standard deviation
            if plt_hist_lines:
                if orientation == 'vertical':
                    ax.axvline(mean, color='r')
                    ax.axvline(mean-2*std_dev, color='r', linestyle='--')
                    ax.axvline(mean+2*std_dev, color='r', linestyle='--')
                elif orientation == 'horizontal':
                    ax.axhline(mean, color='r')
                    ax.axhline(mean-2*std_dev, color='r', linestyle='--')
                    ax.axhline(mean+2*std_dev, color='r', linestyle='--')
            i += 1
            # Add legend to report the total number of points for this instrmt
            lgnd_label = s_instrmt+': '+str(len(df[var_key]))+' points'
            lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label, alpha=mrk_alpha))
            # Add legend handle to report overall statistics
            lgnd_label = 'Mean:'+ '%.4f'%mean+', Std dev:'+'%.4f'%std_dev
            lgnd_hndls.append(mpl.patches.Patch(color=my_clr, label=lgnd_label, alpha=mrk_alpha))
            notes_string = ''.join(df.notes.unique())
        # Only add the notes_string if it contains something
        if len(notes_string) > 1:
            notes_patch  = mpl.patches.Patch(color='none', label=notes_string)
            lgnd_hndls.append(notes_patch)
        # Add legend with custom handles
        lgnd = ax.legend(handles=lgnd_hndls)
        # Invert y-axis if specified
        if y_key in y_invert_vars:
            invert_y_axis = True
        # Add a standard title
        plt_title = add_std_title(a_group)
        return x_label, y_label, plt_title, ax, invert_y_axis
    elif clr_map == 'cluster':
        # Run clustering algorithm
        m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt = get_cluster_args(pp)
        df, rel_val = HDBSCAN_(a_group, df, cl_x_var, cl_y_var, m_pts, min_samp=min_s, extra_cl_vars=[x_key,y_key])
        # Remove rows where the plot variables are null
        df = df[df[var_key].notnull()]
        # Clusters are labeled starting from 0, so total number of clusters is
        #   the largest label plus 1
        n_clusters  = int(df['cluster'].max()+1)
        # Make blank lists to record values
        pts_per_cluster = []
        clstr_means = []
        clstr_stdvs = []
        # Loop through each cluster
        for i in range(n_clusters):
            # Decide on the color and symbol, don't go off the end of the arrays
            my_clr = distinct_clrs[i%len(distinct_clrs)]
            my_mkr = mpl_mrks[i%len(mpl_mrks)]
            # Get relevant data
            this_clstr_df = df[df.cluster == i]
            h_data = this_clstr_df[var_key]
            h_mean = np.mean(h_data)
            # Get histogram parameters
            h_var, res_bins, median, mean, std_dev = get_hist_params(this_clstr_df, var_key, n_h_bins)
            # Plot the histogram
            n, bins, patches = ax.hist(h_var, bins=res_bins, color=my_clr, alpha=mrk_alpha, orientation=orientation, zorder=5)
            # Find the maximum value for this histogram
            h_max = n.max()
            # Find where that max value occured
            bin_h_max = bins[np.where(n == h_max)][0]
            # Plot a symbol to indicate which cluster is which histogram
            if orientation == 'vertical':
                ax.scatter(bin_h_max, h_max, color=my_clr, s=cent_mrk_size, marker=my_mkr, zorder=1)
            elif orientation == 'horizontal':
                ax.scatter(h_max, bin_h_max, color=my_clr, s=cent_mrk_size, marker=my_mkr, zorder=1)
                # ax.scatter(h_max, bin_h_max, color=my_clr, s=cent_mrk_size, marker=r"${}$".format(str(i)), zorder=10)
            #
        # Noise points are labeled as -1
        # Plot noise points
        df_noise = df[df.cluster==-1]
        if plt_noise:
            # Get histogram parameters
            h_var, res_bins, median, mean, std_dev = get_hist_params(df_noise, var_key, n_h_bins)
            # Plot the noise histogram on a twin axis
            if orientation == 'vertical':
                tw_ax = ax.twinx()
                tw_ax.set_ylabel('Number of noise points')
            elif orientation == 'horizontal':
                tw_ax = ax.twiny()
                tw_ax.set_xlabel('Number of noise points')
            n, bins, patches = tw_ax.hist(h_var, bins=res_bins, color=std_clr, alpha=noise_alpha, orientation=orientation, zorder=1)
        n_noise_pts = len(df_noise)
        # Add legend to report the total number of points and notes on the data
        n_pts_patch   = mpl.patches.Patch(color='none', label=str(len(df[var_key]))+' points')
        m_pts_patch = mpl.patches.Patch(color='none', label=r'$m_{pts}$: '+str(m_pts))
        n_clstr_patch = mpl.lines.Line2D([],[],color=cnt_clr, label=r'$n_{clusters}$: '+str(n_clusters), marker='*', linewidth=0)
        n_noise_patch = mpl.patches.Patch(color=std_clr, label=r'$n_{noise pts}$: '+str(n_noise_pts), alpha=noise_alpha, edgecolor=None)
        rel_val_patch = mpl.patches.Patch(color='none', label='DBCV: %.4f'%(rel_val))
        if legend:
            ax.legend(handles=[n_pts_patch, m_pts_patch, n_clstr_patch, n_noise_patch, rel_val_patch])
        # Invert y-axis if specified
        if y_key in y_invert_vars:
            invert_y_axis = True
        # Add a standard title
        plt_title = add_std_title(a_group)
        return x_label, y_label, plt_title, ax, invert_y_axis
    else:
        # Did not provide a valid clr_map
        print('Colormap',clr_map,'not valid')
        exit(0)

################################################################################

def get_hist_params(df, h_key, n_h_bins=25):
    """
    Returns the needed information to make a histogram

    df          A pandas data frame of the data to plot
    h_key       A string of the column header to plot
    n_h_bins    The number of histogram bins to use
    """
    if isinstance(n_h_bins, type(None)):
        n_h_bins = 25
    # Pull out the variable to plot, removing null values
    try:
        h_var = np.array(df[df[h_key].notnull()][h_key])
    except:
        print('Cannot find',h_key,'in the dataframe. Aborting script')
        exit(0)
    # Find overall statistics
    median  = np.median(h_var)
    mean    = np.mean(h_var)
    std_dev = np.std(h_var)
    # Define the bins to use in the histogram, np.arange(start, stop, step)
    start = mean - 3*std_dev
    stop  = mean + 3*std_dev
    step  = std_dev / n_h_bins
    res_bins = np.arange(start, stop, step)
    return h_var, res_bins, median, mean, std_dev

################################################################################

def plot_profiles(ax, a_group, pp, clr_map=None):
    """
    Takes in an Analysis_Group object which has the data and plotting parameters
    to produce a subplot of individual profiles. Returns the x and y labels and
    the subplot title

    ax              The axis on which to make the plot
    a_group         An Analysis_Group object containing the info to create this subplot
    pp              The Plot_Parameters object for a_group
    clr_map         A string to determine what color map to use in the plot
    """
    # Find extra arguments, if given
    legend = pp.legend
    if not isinstance(pp.extra_args, type(None)):
        try:
            plt_noise = pp.extra_args['plt_noise']
        except:
            plt_noise = True
    else:
        plt_noise = True
    scale = pp.plot_scale
    ax_lims = pp.ax_lims
    if scale == 'by_pf':
        print('Cannot use',pp.plot_type,'with plot scale by_pf')
        exit(0)
    clr_map = pp.clr_map
    ## Make a plot of the given profiles vs. the vertical
    # Set the main x and y data keys
    x_key = pp.x_vars[0]
    y_key = pp.y_vars[0]
    var_clr = get_var_color(x_key)
    var_clr = std_clr
    # Check for histogram
    if x_key == 'hist' or y_key == 'hist':
        print('Cannot plot histograms with profiles plot type')
        exit(0)
    # Check for twin x and y data keys
    try:
        tw_x_key = pp.x_vars[1]
        tw_ax_y  = ax.twiny()
        tw_clr = get_var_color(tw_x_key)
        if tw_clr == std_clr:
            tw_clr = alt_std_clr
    except:
        tw_x_key = None
        tw_ax_y  = None
    try:
        tw_y_key = pp.y_vars[1]
        tw_ax_x  = ax.twinx()
        tw_clr = get_var_color(tw_y_key)
        if tw_clr == std_clr:
            tw_clr = alt_std_clr
    except:
        tw_y_key = None
        tw_ax_x  = None
    # Invert y-axis if specified
    if y_key in y_invert_vars:
        invert_y_axis = True
    else:
        invert_y_axis = False
    if not isinstance(tw_y_key, type(None)):
        if tw_y_key in y_invert_vars:
            invert_tw_y_axis = True
        else:
            invert_tw_y_axis = False
    else:
        invert_tw_y_axis = False
    # Whether to add a scale bar or not
    add_scale_bar = False
    # Get extra args dictionary, if it exists
    try:
        extra_args = pp.extra_args
    except:
        extra_args = False
    # Make a blank list for dataframes of each profile
    profile_dfs = []
    # Concatonate all the pandas data frames together
    df = pd.concat(a_group.data_frames)
    # Check whether to run the clustering algorithm
    #   Need to run the clustering algorithm BEFORE filtering to specified range
    cluster_this = False
    #   Find all the plotting variables
    plot_vars = [x_key,y_key,tw_x_key,tw_y_key,clr_map]
    for var in plot_vars:
        if var in clstr_vars:
            cluster_this = True
        #
    if cluster_this:
        m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt = get_cluster_args(pp)
        df, rel_val = HDBSCAN_(a_group, df, cl_x_var, cl_y_var, m_pts, min_samp=min_s, extra_cl_vars=plot_vars)
    # Filter to specified range if applicable
    if not isinstance(a_group.plt_params.ax_lims, type(None)):
        try:
            y_lims = a_group.plt_params.ax_lims['y_lims']
            # Get endpoints of vertical range, make sure they are positive
            y_lims = [abs(ele) for ele in y_lims]
            y_max = max(y_lims)
            y_min = min(y_lims)
            # Filter the data frame to the specified vertical range
            df = df[(df[y_key] < y_max) & (df[y_key] > y_min)]
        except:
            foo = 2
        #
    # Clean out the null values
    df = df[df[x_key].notnull() & df[y_key].notnull()]
    if tw_x_key:
        df = df[df[tw_x_key].notnull()]
    if tw_y_key:
        df = df[df[tw_y_key].notnull()]
    # Get notes
    notes_string = ''.join(df.notes.unique())
    # Check for extra arguments
    if extra_args:
        # Check whether to narrow down the profiles to plot
        try:
            pfs_to_plot = extra_args['pfs_to_plot']
            # Make a temporary list of dataframes for the profiles to plot
            tmp_df_list = []
            for pf_no in pfs_to_plot:
                tmp_df_list.append(df[df['prof_no'] == pf_no])
            df = pd.concat(tmp_df_list)
        except:
            extra_args = None
        try:
            plt_noise = extra_args['plt_noise']
        except:
            plt_noise = True
    # Find the unique profiles for this instrmt
    pfs_in_this_df = np.unique(np.array(df['prof_no']))
    print('\t- Profiles to plot:',pfs_in_this_df)
    # Make sure you're not trying to plot too many profiles
    # if len(pfs_in_this_df) > 15:
    #     print('You are trying to plot',len(pfs_in_this_df),'profiles')
    #     print('That is too many. Try to plot less than 15')
    #     exit(0)
    # else:
    #     print('Plotting',len(pfs_in_this_df),'profiles')
    # Loop through each profile to create a list of dataframes
    for pf_no in pfs_in_this_df:
        # Get just the part of the dataframe for this profile
        pf_df = df[df['prof_no']==pf_no]
        profile_dfs.append(pf_df)
    #
    # Check to see whether any profiles were actually loaded
    n_pfs = len(profile_dfs)
    if n_pfs < 1:
        print('No profiles loaded')
        # Add a standard title
        plt_title = add_std_title(a_group)
        return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis
    # 
    # Plot each profile
    for i in range(n_pfs):
        # Pull the dataframe for this profile
        pf_df = profile_dfs[i]
        # Make a label for this profile
        pf_label = pf_df['source'][0]+pf_df['instrmt'][0]+'-'+str(int(pf_df['prof_no'][0]))
        # Decide on marker and line styles, don't go off the end of the array
        mkr     = mpl_mrks[i%len(mpl_mrks)]
        l_style = l_styles[i%len(l_styles)]
        # Get array to plot and find bounds
        if i == 0:
            # Pull data to plot
            xvar = pf_df[x_key]
            # Find upper and lower bounds of first profile, for reference points
            xvar_low  = min(xvar)
            xvar_high = max(xvar)
            # Find span of the profile
            xv_span = abs(xvar_high - xvar_low)
            # Define far left bound for plot
            left_bound = xvar_low
            # Define pads for x bound (applied at the end)
            x_pad = xv_span/15
            # Find upper and lower bounds of first profile for twin axis
            if tw_x_key:
                tvar = pf_df[tw_x_key]
                twin_low  = min(tvar)
                twin_high = max(tvar)
                tw_span   = abs(twin_high - twin_low)
                # Find normalized profile difference for spacing
                norm_xv = (xvar-xvar_low)/xv_span
                norm_tw = (tvar-twin_low)/tw_span
                norm_pf_diff = min(norm_tw-norm_xv)
                # Shift things over
                right_bound = xvar_high - norm_pf_diff*xv_span
                tw_left_bound = twin_low + norm_pf_diff*tw_span
                # Find index of largest x value
                # xvar_max_idx = np.argmax(xvar)
                # Find value of twin profile at that index
                # tw_xvar_max = tvar[xvar_max_idx]
                tw_x_pad = tw_span/15
            else:
                tvar = None
                right_bound = xvar_high
            #
        # Adjust the starting points of each subsequent profile
        elif i > 0:
            add_scale_bar = True
            # Keep the old profile for reference
            old_xvar = xvar
            old_xvar_low = xvar_low
            old_xvar_high = xvar_high
            old_xv_span = xv_span
            # Find array of data for this profile
            xvar = pf_df[x_key] - min(pf_df[x_key]) + old_xvar_low
            # Find new upper and lower bounds of profile
            xvar_high = max(xvar)
            xvar_low  = min(xvar)
            # Find span of this profile
            xv_span = abs(xvar_high - xvar_low)
            if tw_x_key:
                xvar = xvar + old_xv_span*0.8
                xvar_high = max(xvar)
                xvar_low = min(xvar)
                tvar = pf_df[tw_x_key] - min(pf_df[tw_x_key]) + twin_low + tw_span*0.8
                # 
                twin_low  = min(tvar)
                twin_high = max(tvar)
                tw_span   = abs(twin_high - twin_low)
                # Find normalized profile difference for spacing
                norm_xv = (xvar-xvar_low)/xv_span
                norm_tw = (tvar-twin_low)/tw_span
                norm_pf_diff = min(norm_tw-norm_xv)
                # Shift things over
                right_bound = xvar_high - norm_pf_diff*xv_span
                twin_high = max(tvar)
                # Find index of largest x value
                # xvar_max_idx = np.argmax(xvar)
                # Find value of twin profile at that index
                # tw_xvar_max = tvar[xvar_max_idx]
            else:
                tvar = None
                # Shift things over
                xvar = xvar + old_xv_span*0.45
                xvar_low = min(xvar)
                xvar_high = max(xvar)
                right_bound = xvar_high
            #
        # Determine the color mapping to be used
        if clr_map in a_group.vars_to_keep and clr_map != 'cluster':
            # Format the dates if necessary
            if clr_map == 'dt_start' or clr_map == 'dt_end':
                cmap_data = mpl.dates.date2num(pf_df[clr_map])
            else:
                cmap_data = pf_df[clr_map]
            # Plot a background line for each profile
            ax.plot(xvar, pf_df[y_key], color=var_clr, linestyle=l_style, label=pf_label, zorder=1)
            # Get the colormap
            this_cmap = get_color_map(clr_map)
            # Plot the points as a heatmap
            heatmap = ax.scatter(xvar, pf_df[y_key], c=cmap_data, cmap=this_cmap, s=pf_mrk_size, marker=mkr)
            # Plot on twin axes, if specified
            if not isinstance(tw_x_key, type(None)):
                tw_ax_y.plot(tvar, pf_df[y_key], color=tw_clr, linestyle=l_style, zorder=1)
                tw_ax_y.scatter(tvar, pf_df[y_key], c=cmap_data, cmap=this_cmap, s=pf_mrk_size, marker=mkr)
            #
            # Create the colorbar, but only on the first profile
            if i == 0:
                cbar = plt.colorbar(heatmap, ax=ax)
                # Format the colorbar ticks, if necessary
                if clr_map == 'dt_start' or clr_map == 'dt_end':
                    loc = mpl.dates.AutoDateLocator()
                    cbar.ax.yaxis.set_major_locator(loc)
                    cbar.ax.yaxis.set_major_formatter(mpl.dates.ConciseDateFormatter(loc))
                cbar.set_label(pp.clabel)
        if clr_map == 'clr_all_same':
            # Plot a background line for each profile
            ax.plot(xvar, pf_df[y_key], color=var_clr, linestyle=l_style, label=pf_label, zorder=1)
            # Plot every point the same color, size, and marker
            # ax.scatter(xvar, pf_df[y_key], color=var_clr, s=pf_mrk_size, marker=mkr, alpha=mrk_alpha)
            # Plot on twin axes, if specified
            if not isinstance(tw_x_key, type(None)):
                tw_ax_y.plot(tvar, pf_df[y_key], color=tw_clr, linestyle=l_style, zorder=1)
                # tw_ax_y.scatter(tvar, pf_df[y_key], color=tw_clr, s=pf_mrk_size, marker=mkr, alpha=mrk_alpha)
            #
        if clr_map == 'cluster':
            # Plot a background line for each profile
            ax.plot(xvar, pf_df[y_key], color=var_clr, linestyle=l_style, alpha=0.5, label=pf_label, zorder=1)
            # Make a dataframe with adjusted xvar and tvar
            df_clstrs = pd.DataFrame({x_key:xvar, tw_x_key:tvar, y_key:pf_df[y_key], 'cluster':pf_df['cluster'], 'clst_prob':pf_df['clst_prob']})
            # Get a list of unique cluster numbers, but delete the noise point label "-1"
            cluster_numbers = np.unique(np.array(df_clstrs['cluster'].values, dtype=int))
            cluster_numbers = np.delete(cluster_numbers, np.where(cluster_numbers == -1))
            # print('\tcluster_numbers:',cluster_numbers)
            # Plot noise points first
            if plt_noise:
                ax.scatter(df_clstrs[df_clstrs.cluster==-1][x_key], df_clstrs[df_clstrs.cluster==-1][y_key], color=noise_clr, s=pf_mrk_size, marker=std_marker, alpha=noise_alpha, zorder=2)
                #ax.scatter(df_clstrs[df_clstrs.cluster==-1][x_key], df_clstrs[df_clstrs.cluster==-1][y_key], color=std_clr, s=pf_mrk_size, marker=std_marker, alpha=pf_alpha, zorder=1)
            # Plot on twin axes, if specified
            if not isinstance(tw_x_key, type(None)):
                tw_ax_y.plot(tvar, pf_df[y_key], color=tw_clr, linestyle=l_style, label=pf_label, zorder=1)
                if plt_noise:
                    tw_ax_y.scatter(df_clstrs[df_clstrs.cluster==-1][tw_x_key], df_clstrs[df_clstrs.cluster==-1][y_key], color=std_clr, s=pf_mrk_size, marker=std_marker, alpha=noise_alpha, zorder=2)
            # Loop through each cluster
            for i in cluster_numbers:
                # Decide on the color and symbol, don't go off the end of the arrays
                my_clr = distinct_clrs[i%len(distinct_clrs)]
                my_mkr = mpl_mrks[i%len(mpl_mrks)]
                # print('\t\tcluster:',i,'my_clr:',my_clr,'my_mkr:',my_mkr)
                # Get relevant data
                x_data = df_clstrs[df_clstrs.cluster == i][x_key]
                y_data = df_clstrs[df_clstrs.cluster == i][y_key]
                alphas = df_clstrs[df_clstrs.cluster == i]['clst_prob']
                # Plot the points for this cluster with the specified color, marker, and alpha value
                # ax.scatter(x_data, y_data, color=my_clr, s=pf_mrk_size, marker=my_mkr, alpha=alphas, zorder=5)
                ax.scatter(x_data, y_data, color=my_clr, s=pf_mrk_size, marker=my_mkr, alpha=pf_alpha, zorder=5)
                # Plot on twin axes, if specified
                if not isinstance(tw_x_key, type(None)):
                    t_data = df_clstrs[df_clstrs.cluster == i][tw_x_key]
                    tw_ax_y.scatter(t_data, y_data, color=my_clr, s=pf_mrk_size, marker=my_mkr, alpha=pf_alpha, zorder=5)
        #
    #
    # Plot on twin axes, if specified
    if not isinstance(tw_x_key, type(None)):
        # Adjust bounds on axes
        tw_ax_y.set_xlim([tw_left_bound-tw_x_pad, twin_high+tw_x_pad])
        ax.set_xlim([left_bound-x_pad, right_bound+x_pad])
        # Add label to twin axis
        tw_ax_y.set_xlabel(pp.xlabels[1])
        # Change color of the axis label on the twin axis
        tw_ax_y.xaxis.label.set_color(tw_clr)
        # Change color of the ticks on the twin axis
        tw_ax_y.tick_params(axis='x', colors=tw_clr)
        # Add a grid
        tw_ax_y.grid(color=tw_clr, linestyle='--', alpha=grid_alpha+0.2, axis='x')
        if invert_tw_y_axis:
            tw_ax_y.invert_yaxis()
            print('\t- Inverting twin y axis')
    else:
        ax.set_xlim([left_bound-x_pad, right_bound+x_pad])
    # Check whether to add a scale bar
    if add_scale_bar:
        add_h_scale_bar(ax, ax_lims, unit=' g/kg')
        if tw_x_key:
            add_h_scale_bar(tw_ax_y, ax_lims, unit=r' $^\circ$C', tw_clr=tw_clr)
    # Add legend
    if legend:
        lgnd = ax.legend()
        # Only add the notes_string if it contains something
        if len(notes_string) > 1:
            handles, labels = ax.get_legend_handles_labels()
            handles.append(mpl.patches.Patch(color='none'))
            labels.append(notes_string)
            lgnd = ax.legend(handles=handles, labels=labels)
        # Need to change the marker size for each label in the legend individually
        for hndl in lgnd.legendHandles:
            hndl._sizes = [lgnd_mrk_size]
        #
    # Add a standard title
    plt_title = add_std_title(a_group)
    return pp.xlabels[0], pp.ylabels[0], plt_title, ax, invert_y_axis

################################################################################

def get_cluster_args(pp):
    """
    Finds cluster-realted variables in the extra_args dictionary

    pp              The Plot_Parameters object for a_group
    """
    # Get the dictionary stored in extra_args
    cluster_plt_dict = pp.extra_args
    # print('cluster_plt_dict:',cluster_plt_dict)
    # Get minimum cluster size parameter, if it was given
    try:
        m_pts = int(cluster_plt_dict['m_pts'])
    except:
        m_pts = None
    # Get minimum samples parameter, if it was given
    try:
        min_s = cluster_plt_dict['min_samp']
        if not isinstance(min_s, type(None)):
            min_s = int(min_s)
    except:
        min_s = None
    # Get x and y parameters, if given
    try:
        cl_x_var = cluster_plt_dict['cl_x_var']
        cl_y_var = cluster_plt_dict['cl_y_var']
    except:
        cl_x_var = None
        cl_y_var = None
    # Get the boolean for whether to plot trendline slopes or not
    try:
        plot_slopes = cluster_plt_dict['plot_slopes']
    except:
        plot_slopes = False
    # Check whether or not the box and whisker plots were called for
    try:
        b_a_w_plt = cluster_plt_dict['b_a_w_plt']
    except:
        b_a_w_plt = False
    return m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt

################################################################################

def HDBSCAN_(run_group, df, x_key, y_key, m_pts, min_samp=None, extra_cl_vars=[None]):
    """
    Runs the HDBSCAN algorithm on the set of data specified. Returns a pandas
    dataframe with columns for x_key, y_key, 'cluster', and 'clst_prob' and a
    rough measure of the DBCV score from `relative_validity_`

    run_group   The Analysis_Group object to run HDBSCAN on
    df          A pandas data frame with x_key and y_key as equal length columns
    x_key       String of the name of the column to use on the x-axis
    y_key       String of the name of the column to use on the y-axis
    m_pts      An integer, the minimum number of points for a cluster
    min_samp    An integer, number of points in neighborhood for a core point
    extra_cl_vars   A list of extra variables to potentially calculate
    """
    # print('-- in HDBSCAN')
    # print('-- df columns:',df.columns.values.tolist())
    # If run_group == None, then run the algorithm again
    if isinstance(run_group, type(None)):
        re_run = True
        print('-- run_group is None, re_run:',re_run)
    # Check to make sure there is actually a 'cluster' column in the dataframe
    elif not 'cluster' in df.columns.values.tolist():
        re_run = True
        print('-- `cluster` not in df columns, re_run:',re_run)
    # Check to make sure the global clustering attributes make sense
    else:
        # Make a copy of the global clustering attribute dictionary
        gcattr_dict = {'Last clustered':[],
                       'Clustering x-axis':[],
                       'Clustering y-axis':[],
                       'Clustering m_pts':[],
                       'Clustering filters':[],
                       'Clustering DBCV':[]}
        # Get the global clustering attributes from each dataset
        for ds in run_group.data_set.arr_of_ds:
            for attr in gcattr_dict.keys():
                gcattr_dict[attr].append(ds.attrs[attr])
        # If any attribute has multiple different values, need to re-run HDBSCAN
        re_run = False
        for attr in gcattr_dict.keys():
            gcattr_dict[attr] = list(set(gcattr_dict[attr]))
            if len(gcattr_dict[attr]) > 1:
                re_run = True
                print('-- found multiple distinct attrs in list, re_run:',re_run)
        if gcattr_dict['Last clustered'][0] == 'Never':
            re_run = True
            print('-- `Last clustered` attr is `Never`, re_run:',re_run)
    print('\t- Re-run HDBSCAN:',re_run)
    if re_run:
        print('\t- Running HDBSCAN')
        print('\t\tClustering x-axis:',x_key)
        print('\t\tClustering y-axis:',y_key)
        print('\t\tClustering m_pts: ',m_pts)
        # Set the parameters of the HDBSCAN algorithm
        #   Note: must set gen_min_span_tree=True or you can't get `relative_validity_`
        hdbscan_1 = hdbscan.HDBSCAN(gen_min_span_tree=True, min_cluster_size=m_pts, min_samples=min_samp, cluster_selection_method='leaf')
        # Run the HDBSCAN algorithm
        hdbscan_1.fit_predict(df[[x_key,y_key]])
        # Add the cluster labels and probabilities to the dataframe
        df['cluster']   = hdbscan_1.labels_
        df['clst_prob'] = hdbscan_1.probabilities_
        rel_val = hdbscan_1.relative_validity_
        # Determine whether there are any new variables to calculate
        new_cl_vars = list(set(extra_cl_vars) & set(clstr_vars))
        # Don't need to calculate `cluster` so remove it if its there
        if 'cluster' in new_cl_vars:
            new_cl_vars.remove('cluster')
        if len(new_cl_vars) > 0:
            print('\t\tCalculating extra clustering variables')
            df = calc_extra_cl_vars(df, new_cl_vars)
        # Write this dataframe and DBCV value back to the analysis group
        if not isinstance(run_group, type(None)):
            run_group.data_frames = [df]
            run_group.data_set.arr_of_ds[0].attrs['Clustering DBCV'] = rel_val
        return df, rel_val
    else:
        # Use the clustering results that are already in the dataframe
        print('\t- Using clustering results from file')
        print('\t\tClustering x-axis:',gcattr_dict['Clustering x-axis'])
        print('\t\tClustering y-axis:',gcattr_dict['Clustering y-axis'])
        print('\t\tClustering m_pts: ',gcattr_dict['Clustering m_pts'])
        print('\t\tClustering filters:  ',gcattr_dict['Clustering filters'])
        print('\t\tClustering DBCV:  ',gcattr_dict['Clustering DBCV'])
        # Determine whether there are any new variables to calculate
        new_cl_vars = list(set(extra_cl_vars) & set(clstr_vars))
        # Don't need to calculate `cluster` so remove it if its there
        if 'cluster' in new_cl_vars:
            new_cl_vars.remove('cluster')
        if len(new_cl_vars) > 0:
            print('\t\tCalculating extra clustering variables')
            df = calc_extra_cl_vars(df, new_cl_vars)
        return df, gcattr_dict['Clustering DBCV'][0]

################################################################################

def calc_extra_cl_vars(df, new_cl_vars):
    """
    Takes in an already-clustered pandas data frame and a list of variables and,
    if there are extra variables to calculate, it will add those to the data frame

    df                  A pandas data frame of the data to plot
    new_cl_vars         A list of clustering-related variables to calculate
    """
    # Find the number of clusters
    n_clusters = int(max(df['cluster']+1))
    # Check for variables to calculate
    for this_var in new_cl_vars:
        # Split the prefix from the original variable (assumes an underscore split)
        split_var = this_var.split('_', 1)
        prefix = split_var[0]
        try:
            var = split_var[1]
        except:
            var = None
        # print('prefix:',prefix,'- var:',var)
        # Make a new blank column in the data frame for this variable
        df[this_var] = None
        # Calculate the new values based on the prefix
        if prefix == 'pca':
            # Calculate the per profile cluster average version of the variable
            #   Reduces the number of points to just one per cluster per profile
            # Loop over each profile
            n_pfs = int(max(df['entry'].values))
            for pf in range(n_pfs+1):
                # Find the data from just this profile
                df_this_pf = df[df['entry']==pf]
                # Get a list of clusters in this profile
                clstr_ids = np.unique(np.array(df_this_pf['cluster'].values))
                # Remove the noise points
                clstr_ids = clstr_ids[clstr_ids != -1]
                # Loop over every cluster in this profile
                for i in clstr_ids:
                    # Find the data from this cluster
                    df_this_pf_this_cluster = df_this_pf[df_this_pf['cluster']==i]
                    # Find the mean of this var for this cluster for this profile
                    pf_clstr_mean = np.mean(df_this_pf_this_cluster[var].values)
                    # Put that value back into the original dataframe
                    #   Need to make a mask first for some reason
                    this_pf_this_cluster = (df['entry']==pf) & (df['cluster']==i)
                    df.loc[this_pf_this_cluster, this_var] = pf_clstr_mean
        elif prefix == 'pcs':
            # Calculate the per profile cluster span version of the variable
            #   Reduces the number of points to just one per cluster per profile
            # Loop over each profile
            n_pfs = int(max(df['entry'].values))
            for pf in range(n_pfs+1):
                # Find the data from just this profile
                df_this_pf = df[df['entry']==pf]
                # Get a list of clusters in this profile
                clstr_ids = np.unique(np.array(df_this_pf['cluster'].values))
                # Remove the noise points
                clstr_ids = clstr_ids[clstr_ids != -1]
                # Loop over every cluster in this profile
                for i in clstr_ids:
                    # Find the data from this cluster
                    df_this_pf_this_cluster = df_this_pf[df_this_pf['cluster']==i]
                    # Find the span of this var for this cluster for this profile
                    pf_clstr_span = max(df_this_pf_this_cluster[var].values) - min(df_this_pf_this_cluster[var].values)
                    # Replace any values of zero with `None`
                    if pf_clstr_span == 0:
                        pf_clstr_span = None
                    # Put that value back into the original dataframe
                    #   Need to make a mask first for some reason
                    this_pf_this_cluster = (df['entry']==pf) & (df['cluster']==i)
                    df.loc[this_pf_this_cluster, this_var] = pf_clstr_span
        elif prefix == 'cmc':
            # Calculate the cluster mean-centered version of the variable
            #   Should not change the number of points to display
            # Loop over each cluster
            for i in range(n_clusters):
                # Find the data from this cluster
                df_this_cluster = df[df['cluster']==i].copy()
                # Find the mean of this var for this cluster
                clstr_mean = np.mean(df_this_cluster[var].values)
                # Calculate normalized values
                df_this_cluster[this_var] = df_this_cluster[var] - clstr_mean
                # Put those values back into the original dataframe
                df.loc[df['cluster']==i, this_var] = df_this_cluster[this_var]
        elif prefix == 'ca':
            # Calculate the cluster average version of the variable
            #   Reduces the number of points to just one per cluster
            # Loop over each cluster
            # print('cluster,'+var)
            for i in range(n_clusters):
                # Find the data from this cluster
                df_this_cluster = df[df['cluster']==i].copy()
                # Find the mean of this var for this cluster
                clstr_mean = np.mean(df_this_cluster[var].values)
                # Put those values back into the original dataframe
                df.loc[df['cluster']==i, this_var] = clstr_mean
                # print(str(i)+','+str(clstr_mean))
        elif prefix == 'cs':
            # Calculate the cluster span version of the variable
            #   Reduces the number of points to just one per cluster
            # Loop over each cluster
            for i in range(n_clusters):
                # Find the data from this cluster
                df_this_cluster = df[df['cluster']==i].copy()
                # Find the mean of this var for this cluster
                clstr_span = max(df_this_cluster[var].values) - min(df_this_cluster[var].values)
                # Put those values back into the original dataframe
                df.loc[df['cluster']==i, this_var] = clstr_span
        elif prefix == 'cmm':
            # Find the min/max of each cluster for the variable
            #   Reduces the number of points to just one per cluster
            # Loop over each cluster
            for i in range(n_clusters):
                # Find the data from this cluster
                df_this_cluster = df[df['cluster']==i].copy()
                # Find the min/max of this var for this cluster
                clstr_min = min(df_this_cluster[var].values)
                clstr_max = max(df_this_cluster[var].values)
                # Put those values back into the original dataframe
                #   Won't actually use the data in `this_var` so I'll make it obvious it's to be ignored
                df.loc[df['cluster']==i, this_var] = -999
                df.loc[df['cluster']==i, 'cmin_'+var] = clstr_min
                df.loc[df['cluster']==i, 'cmax_'+var] = clstr_max
            #
        elif prefix == 'nir':
            # Find the normalized inter-cluster range for the variable
            #   Reduces the number of points to just one per cluster, minus one
            #   because it depends on the difference between adjacent clusters
            # Create dataframe to store means and standard deviations of each cluster
            clstr_df = pd.DataFrame(data=np.arange(n_clusters), columns=['cluster'])
            clstr_df['clstr_mean'] = None
            clstr_df['clstr_rnge'] = None
            print('cluster,std_'+var)
            # Loop over each cluster
            for i in range(n_clusters):
                # Find the data from this cluster
                df_this_cluster = df[df['cluster']==i].copy()
                # Find the mean and standard deviation of this var for this cluster
                clstr_mean = np.mean(df_this_cluster[var].values)
                clstr_min = min(df_this_cluster[var].values)
                clstr_max = max(df_this_cluster[var].values)
                clstr_std = np.std(df_this_cluster[var].values)
                # print(i,',',clstr_std)
                clstr_rnge = abs(clstr_max - clstr_min)
                # Put those values into the cluster dataframe
                clstr_df.loc[clstr_df['cluster']==i, 'clstr_mean'] = clstr_mean
                clstr_df.loc[clstr_df['cluster']==i, 'clstr_rnge'] = clstr_rnge
            # Sort the cluster dataframe by the mean values
            sorted_clstr_df = clstr_df.sort_values(by='clstr_mean')
            # Find normalized inter-cluster range for the first cluster in sorted order
            clstr_id_here  = sorted_clstr_df['cluster'].values[0]
            clstr_id_below = sorted_clstr_df['cluster'].values[1]
            clstr_mean_here = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_here, 'clstr_mean'].values[0]
            clstr_mean_below = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_below, 'clstr_mean'].values[0]
            this_diff = abs(clstr_mean_below - clstr_mean_here)
            this_rnge = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_here, 'clstr_rnge'].values[0]
            this_nir = this_rnge / this_diff
            # print('cluster,rnge,diff,nir,diff_below')
            # print(clstr_id_here,',',this_rnge,',',this_diff,',',this_nir,',',this_diff)
            # Put that value back into the original dataframe
            df.loc[df['cluster']==clstr_id_here, this_var] = this_nir
            # Loop over each middle cluster, finding the normalized inter-cluster ranges
            for i in range(1,n_clusters-1):
                # Get cluster id's of this cluster, plus those above and below
                clstr_id_above = sorted_clstr_df['cluster'].values[i-1]
                clstr_id_here  = sorted_clstr_df['cluster'].values[i]
                clstr_id_below = sorted_clstr_df['cluster'].values[i+1]
                # Get the mean value for each of the three clusters for this variable
                clstr_mean_above = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_above, 'clstr_mean'].values[0]
                clstr_mean_here = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_here, 'clstr_mean'].values[0]
                clstr_mean_below = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_below, 'clstr_mean'].values[0]
                # Find the distances from this cluster to the clusters above and below
                diff_above = abs(clstr_mean_above - clstr_mean_here)
                diff_below = abs(clstr_mean_below - clstr_mean_here)
                # Find minimum of the distances above and below
                this_diff = min(diff_above, diff_below)
                # this_diff = np.mean([diff_above, diff_below])
                # Calculate the normalized inter-cluster range for this cluster
                this_rnge = sorted_clstr_df.loc[sorted_clstr_df['cluster']==i, 'clstr_rnge'].values[0]
                this_nir = this_rnge / this_diff
                # print(clstr_id_here,',',this_rnge,',',this_diff,',',this_nir,',',diff_below)
                # Put that value back into the original dataframe
                df.loc[df['cluster']==clstr_id_here, this_var] = this_nir
            # Find normalized inter-cluster range for the last cluster
            clstr_id_above  = sorted_clstr_df['cluster'].values[-2]
            clstr_id_here = sorted_clstr_df['cluster'].values[-1]
            clstr_mean_here = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_here, 'clstr_mean'].values[0]
            clstr_mean_above = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_above, 'clstr_mean'].values[0]
            this_diff = abs(clstr_mean_above - clstr_mean_here)
            this_nir = sorted_clstr_df.loc[sorted_clstr_df['cluster']==clstr_id_here, 'clstr_rnge'].values[0] / this_diff
            # print(clstr_id_here,',',this_rnge,',',this_diff,',',this_nir)
            # Put that value back into the original dataframe
            df.loc[df['cluster']==clstr_id_here, this_var] = this_nir
            #
        if this_var == 'cRL':
            # Find the lateral density ratio R_L for each cluster
            #   Reduces the number of points to just one per cluster
            # Loop over each cluster
            for i in range(n_clusters):
                # Find the data from this cluster
                df_this_cluster = df[df['cluster']==i].copy()
                # Find the variables needed
                alphas = df_this_cluster['alpha'].values
                temps  = df_this_cluster['CT'].values
                betas  = df_this_cluster['beta'].values
                salts  = df_this_cluster['SP'].values
                # Calculate variables needed
                aTs = alphas * temps
                BSs = betas * salts
                # Find the slope of this cluster in aT-BS space
                # m, c = np.linalg.lstsq(np.array([BSs, np.ones(len(BSs))]).T, aTs, rcond=None)[0]
                # Find the slope of the total least-squares of the points for this cluster
                m, c, sd_m, sd_c = orthoregress(BSs, aTs)
                # The lateral density ratio is the inverse of the slope
                this_cRL = 1/m
                # print('cluster:',i,'cRL:',this_cRL)
                # Put those values back into the original dataframe
                df.loc[df['cluster']==i, this_var] = this_cRL
            #
        #
    #
    # Remove rows where the plot variables are null
    for this_var in new_cl_vars:
        if 'nir_' not in this_var:
            df = df[df[this_var].notnull()]
    return df

################################################################################

def find_outliers(df, var_keys, threshold=2):
    """
    Finds any outliers in the dataframe with respect to the x and y keys

    df              A pandas data frame
    var_keys        A list of strings of the names of the columns on which to 
                        find outliers
    threshold       The threshold zscore for which to consider an outlier
    """
    for v_key in var_keys:
        # Get the values of the variable for this key
        v_values = np.array(df[v_key].values, dtype=np.float64)
        # print(v_key,v_values)
        # Find the zscores 
        v_zscores = stats.zscore(v_values)
        # Put those values back into the dataframe
        df['zs_'+ v_key] = v_zscores
        # Make a column of True/False whether that row is an outlier
        df['out_'+v_key] = (df['zs_'+v_key] > threshold) | (df['zs_'+v_key] < -threshold)
    # print(df)
    # exit(0)
    return df

################################################################################

def mark_outliers(ax, df, x_key, y_key, find_all=False, threshold=2, mrk_clr='r', mk_size=mrk_size):
    """
    Finds and marks outliers in the dataframe with respect to the x and y keys
    on the provided axis

    ax              The axis on which to plot the outlier markers
    df              A pandas data frame
    x_key           A string of the variable from df to use on the x axis
    y_key           A string of the variable from df to use on the y axis
    find_all        True/False as whether to find outliers in both cRL and nir_SP
    threshold       The threshold zscore for which to consider an outlier
    mrk_clr         The color in which to mark the outliers
    mk_size         The size of the marker to use to mark the outliers
    """
    print('\t- Marking outliers')
    # Find outliers
    df = find_outliers(df, [x_key, y_key], threshold)
    # If finding outliers in nir or R_L, find outliers in the other as well
    if x_key == 'cRL' and find_all == True:
        df = find_outliers(df, ['nir_SP', y_key], threshold)
        df.loc[df['out_nir_SP']==True, 'out_'+x_key] = True
    #
    # for i in range(len(df)):
    #     this_row = df.loc[df['cluster']==i]
    #     print(this_row['cluster'].values[0], this_row[x_key].values[0], this_row['out_'+x_key].values[0], this_row[y_key].values[0], this_row['out_'+y_key].values[0])#.sort_values('cluster'))
    # Get data with outliers
    x_data = np.array(df[df['out_'+x_key]==True][x_key].values, dtype=np.float64)
    y_data = np.array(df[df['out_'+x_key]==True][y_key].values, dtype=np.float64)
    print('found outliers:')
    print(x_data)
    print(y_data)
    # Mark outliers
    ax.scatter(x_data, y_data, edgecolors=mrk_clr, s=mk_size*5, marker='o', facecolors='none', zorder=2)
    # Run it again
    if mrk_clr == 'r':
        mark_outliers(ax, df.loc[df['out_'+x_key]==False], x_key, y_key, find_all, threshold, mrk_clr='b', mk_size=mk_size)

################################################################################

def plot_clusters(a_group, ax, pp, df, x_key, y_key, cl_x_var, cl_y_var, clr_map, m_pts, min_samp=None, box_and_whisker=True, plot_slopes=False):
    """
    Plots the clusters found by HDBSCAN on the x-y plane

    ax              The axis on which to plot
    pp
    df              A pandas data frame output from HDBSCAN_
    x_key           String of the name of the column to use on the x-axis
    y_key           String of the name of the column to use on the y-axis
    cl_x_var        String of the name of the column used on x-axis of clustering
    cl_y_var        String of the name of the column used on y-axis of clustering
    clr_map         String of the name of the colormap to use (ex: 'clusters')
    m_pts          An integer, the minimum number of points for a cluster
    min_samp        An integer, number of points in neighborhood for a core point
    box_and_whisker True/False whether to include the box and whisker plot
    plot_slopes     True/False whether to plot lines of least-squares slopes for
                        each cluster
    """
    # Needed to make changes to a copy of the global variable
    global mrk_size
    m_size = mrk_size
    # Find extra arguments, if given
    legend = pp.legend
    if not isinstance(pp.extra_args, type(None)):
        try:
            plt_noise = pp.extra_args['plt_noise']
        except:
            plt_noise = True
    else:
        plt_noise = True
    # Decide whether to plot the centroid or not
    if x_key in pf_vars or y_key in pf_vars:
        plot_centroid = False
    elif y_key == 'CT':
        plot_centroid = False
    else:
        plot_centroid = True
    # Run the HDBSCAN algorithm on the provided dataframe
    df, rel_val = HDBSCAN_(a_group, df, cl_x_var, cl_y_var, m_pts, min_samp=min_samp, extra_cl_vars=[x_key,y_key])
    # Clusters are labeled starting from 0, so total number of clusters is
    #   the largest label plus 1
    n_clusters = int(df['cluster'].max()+1)
    # Remove rows where the plot variables are null
    for var in [x_key, y_key]:
        df = df[df[var].notnull()]
    # Noise points are labeled as -1
    # Plot noise points first
    df_noise = df[df.cluster==-1]
    if plt_noise:
        ax.scatter(df_noise[x_key], df_noise[y_key], color=std_clr, s=m_size, marker=std_marker, alpha=noise_alpha, zorder=1)
    n_noise_pts = len(df_noise)
    # Check whether to just plot one point per cluster
    if x_key in ca_vars:
        # Drop duplicates to have just one row per cluster
        df.drop_duplicates(subset=[x_key], keep='first', inplace=True)
        plot_centroid = False
        mrk_size = cent_mrk_size
        plot_slopes = False
    if y_key in ca_vars:
        # Drop duplicates to have just one row per cluster
        df.drop_duplicates(subset=[y_key], keep='first', inplace=True)
        plot_centroid = False
        m_size = cent_mrk_size
        plot_slopes = False
    if 'cRL' in [x_key, y_key]:
        # Drop duplicates to have just one row per cluster
        df.drop_duplicates(subset=['cRL'], keep='first', inplace=True)
        plot_centroid = False
        m_size = cent_mrk_size
        plot_slopes = False
    if 'cRl' in [x_key, y_key]:
        # Drop duplicates to have just one row per cluster
        df.drop_duplicates(subset=['cRl'], keep='first', inplace=True)
        plot_centroid = False
        m_size = cent_mrk_size
        plot_slopes = False
    if x_key in cmm_vars:
        # Split the prefix from the original variable (assumes an underscore split)
        split_var = x_key.split('_', 1)
        var = split_var[1]
        # Drop duplicates to have just one row per cluster
        df.drop_duplicates(subset=['cmin_'+var,'cmax_'+var], keep='first', inplace=True)
        # Calculate
        df['cmm_mid'] = (df['cmax_'+var] + df['cmin_'+var]) / 2
        df['bar_len'] = (df['cmax_'+var] - df['cmin_'+var]) / 2
        x_key = 'cmm_mid'
        plot_centroid = False
        plot_slopes = False
    if y_key in cmm_vars:
        # Split the prefix from the original variable (assumes an underscore split)
        split_var = y_key.split('_', 1)
        var = split_var[1]
        # Drop duplicates to have just one row per cluster
        df.drop_duplicates(subset=['cmin_'+var,'cmax_'+var], keep='first', inplace=True)
        # Calculate
        df['cmm_mid'] = (df['cmax_'+var] + df['cmin_'+var]) / 2
        df['bar_len'] = (df['cmax_'+var] - df['cmin_'+var]) / 2
        y_key = 'cmm_mid'
        plot_centroid = False
        plot_slopes = False
    # Look for outliers
    if 'nir' in x_key or 'cRL' in x_key or 'ca' in x_key:
        mrk_outliers = True
    else:
        mrk_outliers = False
    # Set variable as to whether to invert the y axis
    invert_y_axis = False
    # Which colormap?
    if clr_map == 'cluster':
        # Make blank lists to record values
        pts_per_cluster = []
        clstr_means = []
        clstr_stdvs = []
        # Loop through each cluster
        for i in range(n_clusters):
            # Decide on the color and symbol, don't go off the end of the arrays
            my_clr = distinct_clrs[i%len(distinct_clrs)]
            my_mkr = mpl_mrks[i%len(mpl_mrks)]
            # Find the data from this cluster
            df_this_cluster = df[df['cluster']==i]
            # print(df_this_cluster)
            # Get relevant data
            x_data = df_this_cluster[x_key] 
            x_mean = np.mean(x_data)
            x_stdv = np.std(x_data)
            y_data = df_this_cluster[y_key] 
            y_mean = np.mean(y_data)
            y_stdv = np.std(y_data)
            alphas = df_this_cluster['clst_prob'] #df[df.cluster == i]['clst_prob']
            # Plot the points for this cluster with the specified color, marker, and alpha value
            #   Having an issue with actually using the alphas from above without a TypeError
            if x_key == 'cmm_mid':
                xerrs = df[df.cluster == i]['bar_len']
                ax.errorbar(x_data, y_data, xerr=xerrs, color=my_clr, capsize=l_cap_size)
            elif y_key == 'cmm_mid':
                yerrs = df[df.cluster == i]['bar_len']
                ax.errorbar(x_data, y_data, yerr=yerrs, color=my_clr, capsize=l_cap_size)
            elif 'nir' in x_key or 'cRL' in x_key or 'cRl' in x_key:
                ax.scatter(x_data, y_data, color=my_clr, s=m_size, marker=my_mkr, alpha=1, zorder=5)
            else:
                ax.scatter(x_data, y_data, color=my_clr, s=m_size, marker=my_mkr, alpha=mrk_alpha, zorder=5)
            # Plot the centroid of this cluster
            if plot_centroid:
                # This will plot a marker at the centroid
                ax.scatter(x_mean, y_mean, color=std_clr, s=cent_mrk_size*1.1, marker='o', zorder=9)
                ax.scatter(x_mean, y_mean, color=cnt_clr, s=cent_mrk_size, marker=my_mkr, zorder=10)
                # This will plot the cluster number at the centroid
                # ax.scatter(x_mean, y_mean, color=cnt_clr, s=cent_mrk_size, marker=r"${}$".format(str(i)), zorder=10)
            if plot_slopes:
                # Find the slope of the ordinary least-squares of the points for this cluster
                # m, c = np.linalg.lstsq(np.array([x_data, np.ones(len(x_data))]).T, y_data, rcond=None)[0]
                # Find the slope of the total least-squares of the points for this cluster
                m, c, sd_m, sd_c = orthoregress(x_data, y_data)
                # Plot the least-squares fit line for this cluster through the centroid
                ax.axline((x_mean, y_mean), slope=m, color=my_clr, zorder=3)
                # Add annotation to say what the slope is
                ax.annotate('%.2f'%(1/m), xy=(x_mean+x_stdv/4,y_mean+y_stdv/10), xycoords='data', color=std_clr, weight='bold', zorder=12)
            #
            # Record the number of points in this cluster
            pts_per_cluster.append(len(x_data))
            # Record mean and standard deviation to calculate normalized inter-cluster range
            clstr_means.append(y_mean)
            clstr_stdvs.append(y_stdv)
        # Mark outliers, if specified
        if mrk_outliers:
            mark_outliers(ax, df, x_key, y_key, mk_size=m_size, mrk_clr='red')
        if 'nir_' in x_key:
            # Get bounds of axes
            x_bnds = ax.get_xbound()
            y_bnds = ax.get_ybound()
            x_span = abs(x_bnds[1] - x_bnds[0])
            y_span = abs(y_bnds[1] - y_bnds[0])
            # Add line at R_L = -1
            ax.axvline(1, color=std_clr, alpha=0.5, linestyle='-')
            ax.annotate(r'$IR_{S_P}=1$', xy=(1+x_span/15,y_bnds[0]+y_span/5), xycoords='data', color=std_clr, weight='bold', alpha=0.5, zorder=12)
        # Add some lines if plotting cRL
        if 'cRL' in [x_key, y_key]:
            # Get bounds of axes
            x_bnds = ax.get_xbound()
            y_bnds = ax.get_ybound()
            x_span = abs(x_bnds[1] - x_bnds[0])
            y_span = abs(y_bnds[1] - y_bnds[0])
            if x_key == 'cRL':
                # Plot polynomial fit line
                if True:
                    # Get the data without the outliers
                    x_data = df[df['out_'+x_key] == False][x_key].astype('float')
                    y_data = df[df['out_'+x_key] == False][y_key].astype('float')
                    # Use polyfit 
                    deg = 2
                    z = np.polyfit(y_data, x_data, deg)
                    x_fit = np.poly1d(z)
                    # Get the means and standard deviations
                    x_mean = np.mean(x_data)
                    x_stdv = np.std(x_data)
                    y_mean = np.mean(y_data)
                    y_stdv = np.std(y_data)
                    # Find R^2 (coefficient of determination) value which is 1 
                    #   minus the ratio of the sum of squares of residuals over 
                    #   the total sum of squares
                    x_fit_data = x_fit(y_data)
                    ss_res = ((x_data-x_fit_data)**2).sum()
                    ss_tot = ((x_data-x_mean)**2).sum()
                    R2 = 1 - ss_res/ss_tot
                    print('\t- Sum of square residuals:',ss_res)
                    print('\t- Total sum of squares:',ss_tot)
                    print('\t- R^2:',R2)
                    print('\t- Median cRL:',np.median(x_data))
                    print('\t- Mean cRL:',np.mean(x_data))
                    print('\t- Standard deviation of cRL:',np.std(x_data))
                    # Make an array of values to plot log fit line 
                    y_arr = np.linspace(y_bnds[0], y_bnds[1], 50)
                    # Plot the fit
                    ax.plot(x_fit(y_arr), y_arr, color=alt_std_clr)
                    # Add annotation to say what the line is
                    line_label = "$R_L = " + format_sci_notation(z[0]) + " p^2 + " + format_sci_notation(z[1])+ "p " + format_sci_notation(z[2])+"$"
                    # Find a good place for the annotation
                    if x_bnds[1] < 0:
                        ann_x_placement = x_span/10+min(x_bnds)
                    else:
                        ann_x_placement = x_mean+0.5*x_stdv
                    ax.annotate(line_label, xy=(ann_x_placement,y_span/2.7+min(y_bnds)), xycoords='data', color=alt_std_clr, weight='bold', zorder=12)
                    # For Figure 5, add the line of the opposite ITP
                    print('format_sci_notation(z[2]):',format_sci_notation(z[2]))
                    if format_sci_notation(z[2]) == "-39.85":
                        print('\t- Adding additional curve for ITP3')
                        ax.plot(-3.38e-4*(y_arr)**2 + 0.24*y_arr - 43.29, y_arr, color=std_clr, alpha=0.5, linestyle='--', zorder=1)
                    elif format_sci_notation(z[2]) == "-43.29":
                        print('\t- Adding additional curve for ITP2')
                        ax.plot(-3.81e-4*(y_arr)**2 + 0.24*y_arr - 39.85, y_arr, color=std_clr, alpha=0.5, linestyle='--', zorder=1)
                    # Limit the y axis
                    ax.set_ylim((y_bnds[1],y_bnds[0]))
                # Plot exponential fit line
                if False:
                    # Get the data without the outliers
                    x_data = df[df['out_'+x_key] == False][x_key].astype('float')
                    y_data = df[df['out_'+x_key] == False][y_key].astype('float')
                    # Get the means and standard deviations
                    x_mean = np.mean(x_data)
                    x_stdv = np.std(x_data)
                    y_mean = np.mean(y_data)
                    y_stdv = np.std(y_data)
                    # Make an array of values to plot log fit line (needs to be all negative values)
                    x_arr = np.linspace(x_bnds[0], 0, 25)
                    # Find the slope of the total least-squares of the points in semilogx
                    m, c, sd_m, sd_c = orthoregress(np.log(-1*x_data), y_data)
                    print('sd_m:',sd_m,'sd_c:',sd_c)
                    # Plot the semilogx least-squares fit on the non-semilogx axes
                    ax.plot(x_arr, m*np.log(-1*x_arr)+c, color=alt_std_clr)
                    # Add annotation to say what the line is
                    # line_label = "$%.2f\\ln(-R_L)+%.2f $"%(m,c)
                    line_label = "$R_L = -\exp( %.2f - p/ %.2f )$"%(c/abs(m),abs(m))
                    ax.annotate(line_label, xy=(x_mean+2*x_stdv,y_mean+y_stdv), xycoords='data', color=alt_std_clr, weight='bold', zorder=12)
                    # Limit the y axis
                    ax.set_ylim((y_bnds[1],y_bnds[0]))
                #
            #
        # Add legend to report the total number of points and notes on the data
        n_pts_patch   = mpl.patches.Patch(color='none', label=str(len(df[x_key]))+' points')
        m_pts_patch = mpl.patches.Patch(color='none', label=r'$m_{pts}$: '+str(m_pts))
        n_clstr_patch = mpl.lines.Line2D([],[],color=cnt_clr, label=r'$n_{clusters}$: '+str(n_clusters), marker='*', linewidth=0)
        n_noise_patch = mpl.patches.Patch(color=std_clr, label=r'$n_{noise pts}$: '+str(n_noise_pts), alpha=noise_alpha, edgecolor=None)
        rel_val_patch = mpl.patches.Patch(color='none', label='DBCV: %.4f'%(rel_val))
        if legend:
            ax.legend(handles=[n_pts_patch, m_pts_patch, n_clstr_patch, n_noise_patch, rel_val_patch])
        # Add inset axis for box-and-whisker plot
        if box_and_whisker:
            inset_ax = inset_axes(ax, width="30%", height="10%", loc=4)
            inset_ax.boxplot(pts_per_cluster, vert=0, sym='X')
            #   Adjust the look of the inset axis
            inset_ax.xaxis.set_ticks_position('top')
            inset_ax.set_xticks([int(min(pts_per_cluster)),int(np.median(pts_per_cluster)),int(max(pts_per_cluster))])
            inset_ax.spines['left'].set_visible(False)
            inset_ax.get_yaxis().set_visible(False)
            inset_ax.spines['right'].set_visible(False)
            inset_ax.spines['bottom'].set_visible(False)
            inset_ax.set_xlabel('Points per cluster')
            inset_ax.xaxis.set_label_position('top')
        #
        # Invert y-axis if specified
        if y_key in y_invert_vars:
            invert_y_axis = True
    #
    return invert_y_axis

################################################################################

def plot_clstr_param_sweep(ax, tw_ax_x, a_group, plt_title=None):
    """
    Plots the number of clusters found by HDBSCAN vs. the number of profiles
    included in the data set

    ax              The axis on which to plot
    tw_ax_x         The twin x axis on which to plot
    a_group         An Analysis_Group object containing the info to create this subplot
    plt_title       A string to use as the title of this subplot
    """
    ## Get relevant parameters for the plot
    pp = a_group.plt_params
    # Concatonate all the pandas data frames together
    df = pd.concat(a_group.data_frames)
    # Get the dictionary stored in extra_args
    cluster_plt_dict = pp.extra_args
    # Set the main x and y data keys
    x_key = pp.x_vars[0]
    y_key = pp.y_vars[0]
    # Check for twin axis data key
    try:
        tw_y_key = pp.y_vars[1]
    except:
        tw_y_key = None
    # Get cluster arguments
    m_pts, min_s, cl_x_var, cl_y_var, plot_slopes, b_a_w_plt = get_cluster_args(pp)
    # Get the parameter sweep tuple, expect the form [start,stop,step]
    try:
        cl_ps_tuple = cluster_plt_dict['cl_ps_tuple']
    except:
        print('Cannot find `cl_ps_tuple` in `extra_args`, aborting script')
        exit(0)
    # Build the array for the x_var axis
    x_var_array = np.arange(cl_ps_tuple[0], cl_ps_tuple[1], cl_ps_tuple[2])
    # Check for a z variable
    try:
        z_key = cluster_plt_dict['z_var']
        z_list = cluster_plt_dict['z_list']
    except:
        z_key = None
        z_list = [0]
    # Open a text file to record values from the parameter sweep
    sweep_txt_file = 'outputs/ps_x_'+x_key+'_z_'+str(z_key)+'.txt'
    f = open(sweep_txt_file,'w')
    f.write('Parameter Sweep for '+plt_title+'\n')
    f.write(datetime.now().strftime("%I:%M%p on %B %d, %Y"))
    f.close()
    # If limiting the number of pfs, find total number of pfs in the given df
    #   In the multi-index of df, level 0 is 'Time'
    if x_key == 'n_pfs':
        pf_nos = np.unique(np.array(df['prof_no'].values))
        number_of_pfs = len(pf_nos)
        print('\tNumber of profiles:',number_of_pfs)
        x_var_array = x_var_array[x_var_array <= number_of_pfs]
    if z_key == 'n_pfs':
        pf_nos = np.unique(np.array(df['prof_no'].values))
        number_of_pfs = len(pf_nos)
        print('\tNumber of profiles:',number_of_pfs)
        z_list = np.array(z_list)
        z_list = z_list[z_list <= number_of_pfs]
    #
    print('\tPlotting these x values of',x_key,':',x_var_array)
    f = open(sweep_txt_file,'a')
    f.write('\nPlotting these x values of '+x_key+':\n')
    f.write(str(x_var_array))
    f.close()
    if z_key:
        print('\tPlotting these z values of',z_key,':',z_list)
        f = open(sweep_txt_file,'a')
        f.write('\nPlotting these z values of '+z_key+':\n')
        f.write(str(z_list))
        f.close()
    for i in range(len(z_list)):
        y_var_array = []
        tw_y_var_array = []
        for x in x_var_array:
            # Set initial values for some variables
            zlabel = None
            this_df = df.copy()
            # Set parameters based on variables selected
            #   NOTE: need to run `ell_size` BEFORE `n_pfs`
            if x_key == 'ell_size':
                # Need to apply moving average window to original data, before
                #   the data filters were applied, so make a new Analysis_Group
                a_group.profile_filters.m_avg_win = x
                new_a_group = Analysis_Group(a_group.data_set, a_group.profile_filters, a_group.plt_params)
                this_df = pd.concat(new_a_group.data_frames)
                xlabel = r'$\ell$ (dbar)'
            if z_key == 'ell_size':
                # Need to apply moving average window to original data, before
                #   the data filters were applied, so make a new Analysis_Group
                a_group.profile_filters.m_avg_win = z_list[i]
                new_a_group = Analysis_Group(a_group.data_set, a_group.profile_filters, a_group.plt_params)
                this_df = pd.concat(new_a_group.data_frames)
                zlabel = r'$\ell=$'+str(z_list[i])+' dbar'
            if x_key == 'n_pfs':
                this_df = this_df[this_df['prof_no'] <= pf_nos[x-1]].copy()
                xlabel = 'Number of profiles included'
            if z_key == 'n_pfs':
                this_df = this_df[this_df['prof_no'] <= pf_nos[z_list[i]-1]].copy()
                zlabel = r'$n_{profiles}=$'+str(z_list[i])
            if x_key == 'm_pts':
                # min cluster size must be an integer
                m_pts = int(x)
                xlabel = r'$m_{pts}$'
            elif x_key == 'min_samps':
                # min samples must be an integer, or None
                if not isinstance(x, type(None)):
                    min_s = int(x)
                else:
                    min_s = x
                xlabel = 'Minimum samples'
            if z_key == 'm_pts':
                # min cluster size must be an integer
                m_pts = int(z_list[i])
                zlabel = r'$m_{pts}=$: '+str(m_pts)
            elif z_key == 'min_samps':
                # min samps must be an integer, or None
                if not isinstance(z_list[i], type(None)):
                    min_s = int(z_list[i])
                else:
                    min_s = z_list[i]
                zlabel = 'Minimum samples: '+str(min_s)
            # Run the HDBSCAN algorithm on the provided dataframe
            new_df, rel_val = HDBSCAN_(None, this_df, cl_x_var, cl_y_var, m_pts, min_samp=min_s)
            # Record outputs to plot
            if y_key == 'DBCV':
                # relative_validity_ is a rough measure of DBCV
                y_var_array.append(rel_val)
                ylabel = 'DBCV'
            elif y_key == 'n_clusters':
                # Clusters are labeled starting from 0, so total number of clusters is
                #   the largest label plus 1
                y_var_array.append(new_df['cluster'].max()+1)
                ylabel = 'Number of clusters'
            if tw_y_key:
                if tw_y_key == 'DBCV':
                    # relative_validity_ is a rough measure of DBCV
                    tw_y_var_array.append(rel_val)
                    tw_ylabel = 'DBCV'
                elif tw_y_key == 'n_clusters':
                    # Clusters are labeled starting from 0, so total number of clusters is
                    #   the largest label plus 1
                    tw_y_var_array.append(new_df['cluster'].max()+1)
                    tw_ylabel = 'Number of clusters'
                #
            #
            f = open(sweep_txt_file,'a')
            if z_key:
                f.write('\n m_pts: '+str(m_pts)+' '+str(x_key)+': '+str(x)+' '+str(z_key)+': '+str(z_list[i])+' n_clstrs: '+str(new_df['cluster'].max()+1)+' DBCV: '+str(rel_val))
            else:
                f.write('\n m_pts: '+str(m_pts)+' '+str(x_key)+': '+str(x)+' n_clstrs: '+str(new_df['cluster'].max()+1)+' DBCV: '+str(rel_val))
            f.close()
        ax.plot(x_var_array, y_var_array, color=std_clr, linestyle=l_styles[i], label=zlabel)
        # Add gridlines
        ax.grid(color=std_clr, linestyle='--', alpha=grid_alpha)
        if tw_y_key:
            if z_key:
                tw_ax_x.plot(x_var_array, tw_y_var_array, color=alt_std_clr, linestyle=l_styles[i])
            else: # If not plotting multiple z values, make the twin axis line a different linestyle
                tw_ax_x.plot(x_var_array, tw_y_var_array, color=alt_std_clr, linestyle=l_styles[i+1])
            tw_ax_x.set_ylabel(tw_ylabel)
            # Change color of the axis label on the twin axis
            tw_ax_x.yaxis.label.set_color(alt_std_clr)
            # Change color of the ticks on the twin axis
            tw_ax_x.tick_params(axis='y', colors=alt_std_clr)
            # Add gridlines
            # tw_ax_x.grid(color=alt_std_clr, linestyle='--', alpha=grid_alpha+0.3, axis='y')
        f = open(sweep_txt_file,'a')
        f.write('\n')
        f.close()
    if z_key:
        ax.legend()
    return xlabel, ylabel

################################################################################

################################################################################
