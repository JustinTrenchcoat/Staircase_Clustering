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

### Filters for reproducing plots from Timmermans et al. 2008
# Timmermans 2008 Figure 4 depth range
T2008_fig4_y_lims = {'y_lims':[-260,-220]}
# Timmermans 2008 Figure 4 shows profile 185
T2008_fig4_pfs = [183, 185, 187]
# Filters used in Timmermans 2008 T-S and aT-BS plots
T2008_p_range = [180,300]
T2008_S_range = [34.4,34.6]
T2008_fig5a_ax_lims = {'x_lims':[34.05,34.75], 'y_lims':[-1.3,0.5]}
T2008_fig6a_ax_lims = {'x_lims':[0.027,0.027045], 'y_lims':[-13e-6,3e-6]}
# The actual limits are above, but need to adjust the x lims for some reason
T2008_fig6a_ax_lims = {'x_lims':[0.026835,0.026880], 'y_lims':[-13e-6,3e-6]}

################################################################################
# Make dictionaries for what data to load in and analyze
################################################################################

# All profiles from certain ITPs
ITP2_all = {'ITP_2':'all'}

# Just specific profiles
ITP2_pfs = {'ITP_2':T2008_fig4_pfs}

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
pfs_f1 = ahf.Profile_Filters(p_range=staircase_range, m_avg_win=50)

################################################################################
# Create plotting parameter objects

print('- Creating plotting parameter objects')
################################################################################

### Test plots
pp_xy_default = ahf.Plot_Parameters()
pp_test0 = ahf.Plot_Parameters(x_vars=['BSP'], y_vars=['aiT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':90}, ax_lims=T2008_fig6a_ax_lims)
# pp_test = ahf.Plot_Parameters(x_vars=['iT', 'la_iT'], y_vars=['press'], plot_type='profiles', clr_map='clr_all_same')
# pp_test = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_iT'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':200})

### Figures for paper
## Parameter sweeps
pp_ps_min_cs = ahf.Plot_Parameters(x_vars=['min_cs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':90, 'cl_ps_tuple':[20,210,10], 'z_var':'maw_size', 'z_list':[50,70,140]})
pp_ps_l_maw  = ahf.Plot_Parameters(x_vars=['maw_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':90, 'cl_ps_tuple':[10,210,10], 'z_var':'min_cs', 'z_list':[90,120,240]})
pp_ps_n_pfs  = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_iT', 'min_cs':90, 'cl_ps_tuple':[10,410,10], 'z_var':'min_cs', 'z_list':[90,120,240]})

## Reproducing Timmermans et al. 2008 Figure 5a, but with cluster coloring
pp_T2008_fig5a = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, ax_lims=T2008_fig5a_ax_lims)
## Reproducing Timmermans et al. 2008 Figure 6a, but with cluster coloring
pp_T2008_fig6a = ahf.Plot_Parameters(x_vars=['BSP'], y_vars=['aCT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, ax_lims=T2008_fig6a_ax_lims)

################################################################################
# Create analysis group objects
print('- Creating analysis group objects')
################################################################################

## Test Analysis Groups
my_group0 = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_test0)
# my_group1 = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_test)
# my_group0 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_f0, pp_test)
# my_group1 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_f1, pp_test)

### Figures for paper
## Parameter sweeps
# group_ps_min_cs = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_ps_min_cs)
# group_ps_l_maw  = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_ps_l_maw)
# group_ps_n_pfs  = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_ps_n_pfs)
## Reproducing figures from Timmermans et al. 2008
group_T2008_fig5a = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_T2008_fig5a)
group_T2008_fig6a = ahf.Analysis_Group(ds_ITP2_all, pfs_f1, pp_T2008_fig6a)

################################################################################
# Declare figures or summaries to output
print('- Creating outputs')
################################################################################

# ahf.make_figure([my_group0])
# ahf.make_figure([my_group0, my_group1])

### Figures for paper
## Parameter sweeps
# ahf.make_figure([group_ps_min_cs, group_ps_l_maw, group_ps_n_pfs])
## Reproducing figures from Timmermans et al. 2008
ahf.make_figure([group_T2008_fig5a, group_T2008_fig6a])
