"""
Author: Mikhail Schee
Created: 2022-08-18

This script is set up to make figures of Arctic Ocean profile data that has been
formatted into netcdfs by the `make_netcdf` function.

"""

# For custom analysis functions
import analysis_helper_functions as ahf

# Common filters
staircase_range = [200,300]

################################################################################
# Make dictionaries for what data to load in and analyze
################################################################################

# All profiles from certain ITPs
ITP2_all = {'ITP_2':'all'}

# Just specific profiles
ITP2_pfs = {'ITP_2':[183, 185, 187]}

################################################################################
# Create data filtering objects
print('- Creating data filtering objects')
################################################################################

dfs0 = ahf.Data_Filters()

################################################################################
# Create data sets by combining filters and the data to load in
print('- Creating data sets')
################################################################################

ds_ITP2_all = ahf.Data_Set(ITP2_all, dfs0)
ds_ITP2_pfs = ahf.Data_Set(ITP2_pfs, dfs0)

################################################################################
# Create profile filtering objects
print('- Creating profile filtering objects')
################################################################################

pfs_f0 = ahf.Profile_Filters()
pfs_f1 = ahf.Profile_Filters(p_range=staircase_range)

################################################################################
# Create plotting parameter objects

print('- Creating plotting parameter objects')
################################################################################

### xy plots
pp_xy_default = ahf.Plot_Parameters()
pp_test = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['DBCV','n_clusters'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':200, 'cl_ps_tuple':[10,300,50], 'z_var':'min_cs', 'z_list':[100,150,200]})
pp_test0 = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_iT'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':200})
# pp_test = ahf.Plot_Parameters(x_vars=['iT'], y_vars=['press'], plot_type='profiles', clr_map='clr_all_same')

################################################################################
# Create analysis group objects
print('- Creating analysis group objects')
################################################################################

# Analysis Groups
my_group0 = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_test)
my_group1 = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_f1, pp_test)

################################################################################
# Declare figures or summaries to output
print('- Creating outputs')
################################################################################

ahf.make_figure([my_group0, my_group1])
