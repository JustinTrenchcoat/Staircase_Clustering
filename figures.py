"""
Author: Mikhail Schee
Created: 2022-08-18

This script is set up to make figures of Arctic Ocean profile data that has been
formatted into netcdfs by the `make_netcdf` function.

Note: You unfortunately cannot pickle plots with parasite axes. So, the pickling
functionality won't work if you try to make a plot with the box and whisker plot
in an inset.

"""

# For custom analysis functions
import analysis_helper_functions as ahf

### Filters for reproducing plots from Timmermans et al. 2008
ITP2_p_range = [185,300]
ITP2_S_range = [34.05,34.75]
T2008_m_pts = 170
# Timmermans 2008 Figure 4 depth range
T2008_fig4_y_lims = {'y_lims':[260,220]}
# Timmermans 2008 Figure 4 shows profile 185
T2008_fig4_pfs = [183, 185, 187]
# Filters used in Timmermans 2008 T-S and aT-BS plots
T2008_p_range = [180,300]
T2008_S_range = [34.4,34.6]
T2008_fig5a_ax_lims = {'x_lims':[34.05,34.75], 'y_lims':[-1.3,0.5]}
T2008_fig6a_ax_lims = {'x_lims':[0.027002,0.027042], 'y_lims':[-13e-6,3e-6]}
# The actual limits are above, but need to adjust the x lims for some reason
T2008_fig6a_ax_lims = {'x_lims':[0.026838,0.026878], 'y_lims':[-13e-6,3e-6]}

# A list of many profiles to plot from ITP2
start_pf = 85
import numpy as np
n_pfs_to_plot = 15
ITP2_some_pfs = list(np.arange(start_pf, start_pf+(n_pfs_to_plot*2), 2))
# For showing multiple layers grouped into one cluster
ITP2_some_pfs_0 = [87, 89, 95, 97, 99, 101, 103, 105, 109, 111]
ITP2_some_pfs_ax_lims_0 = {'y_lims':[240,220]}
# For showing one layer split into multiple clusters
ITP2_some_pfs_1 = [67, 69, 73, 75, 81, 83, 91, 93, 97, 99]
ITP2_some_pfs_ax_lims_1 = {'y_lims':[295,270]}

### Filters for reproducing plots from Lu et al. 2022
Lu2022_p_range = [200,355]
Lu2022_T_range = [-1.0,0.9]
Lu2022_S_range = [34.21,34.82]
Lu2022_m_pts = 580

################################################################################
# Make dictionaries for what data to load in and analyze
################################################################################

# All profiles from all ITPs in this study
all_ITPs = {'ITP_2':'all','ITP_3':'all'}

# All profiles from certain ITPs
ITP2_all  = {'ITP_2':'all'}
ITP3_all  = {'ITP_3':'all'}

# Just specific profiles
ITP2_pfs  = {'ITP_2':[185]}#T2008_fig4_pfs}

ITP2_pfs1  = {'ITP_2':ITP2_some_pfs}

################################################################################
# Create data filtering objects
print('- Creating data filtering objects')
################################################################################

dfs0 = ahf.Data_Filters()
dfs_S2019 = ahf.Data_Filters(date_range=['2007/09/01','2007/12/31'])
dfs_test = ahf.Data_Filters(date_range=['2005/08/25 00:00:00','2005/10/31 00:00:00'])
# dfs_test = ahf.Data_Filters(date_range=['2004/09/01 00:00:00','2004/09/15 00:00:00'])

################################################################################
# Create data sets by combining filters and the data to load in
print('- Creating data sets')
################################################################################

ds_all_ITPs = ahf.Data_Set(all_ITPs, dfs0)

ds_ITP2_all = ahf.Data_Set(ITP2_all, dfs0)
ds_ITP2_pfs = ahf.Data_Set(ITP2_pfs, dfs0)

ds_ITP2_some_pfs = ahf.Data_Set(ITP2_pfs1, dfs0)

ds_ITP3_all = ahf.Data_Set(ITP3_all, dfs0)

################################################################################
# Create profile filtering objects
print('- Creating profile filtering objects')
################################################################################

pfs_0 = ahf.Profile_Filters()
pfs_1 = ahf.Profile_Filters(p_range=T2008_p_range)

test_p_range = [200,400]
pfs_maw_10  = ahf.Profile_Filters(p_range=test_p_range, m_avg_win=10)
pfs_maw_50  = ahf.Profile_Filters(p_range=test_p_range, m_avg_win=50)
pfs_maw_100 = ahf.Profile_Filters(p_range=test_p_range)
pfs_maw_150 = ahf.Profile_Filters(p_range=test_p_range, m_avg_win=150)

pfs_T2008  = ahf.Profile_Filters(SP_range=ITP2_S_range)
pfs_Lu2022 = ahf.Profile_Filters(SP_range=Lu2022_S_range)

################################################################################
# Create plotting parameter objects

print('- Creating plotting parameter objects')
################################################################################

### Test plots
pp_xy_default = ahf.Plot_Parameters()
pp_test0 = ahf.Plot_Parameters(x_vars=['alpha'], y_vars=['beta'], clr_map='clr_all_same', legend=True, extra_args={'b_a_w_plt':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':100})
pp_og_ma_la_pf = ahf.Plot_Parameters(x_vars=['CT', 'ma_CT'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':[185], 'plt_noise':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':210}, ax_lims=T2008_fig4_y_lims)
pp_test1 = ahf.Plot_Parameters(x_vars=['cRL'], y_vars=['ca_press'], clr_map='clr_all_same', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts, 'plot_slopes':True}, legend=True)
# pp_test1 = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['hist'], clr_map='cluster', legend=True, extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':90})
# pp_test1 = ahf.Plot_Parameters(x_vars=['pcs_press'], y_vars=['pca_press'], clr_map='density_hist', legend=True, extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':90, 'clr_min':0, 'clr_max':15, 'clr_ext':'max', 'xy_bins':250})
# pp_test0 = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['lat'], clr_map='clr_all_same')
pp_map = ahf.Plot_Parameters(plot_type='map', clr_map='dt_start')

pp_pfs0 = ahf.Plot_Parameters(x_vars=['CT','ma_CT'], y_vars=['depth'], plot_type='profiles')
pp_pfs1 = ahf.Plot_Parameters(x_vars=['CT','ma_CT'], y_vars=['depth'], first_dfs=[True], plot_type='profiles')
pp_pfs2 = ahf.Plot_Parameters(x_vars=['CT','ma_CT'], y_vars=['depth'], finit_dfs=[True], plot_type='profiles')
pp_pfs3 = ahf.Plot_Parameters(x_vars=['sigma','ma_sigma'], y_vars=['depth'], finit_dfs=[True], plot_type='profiles')

## Test Figures
# The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_Lu2022_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts}, legend=True)

## Testing impact of maw_size
pp_CT_SP  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='clr_all_same', legend=False)
pp_CT_pfs = ahf.Plot_Parameters(x_vars=['CT'], y_vars=['press'], plot_type='profiles', clr_map='clr_all_same', extra_args={'pfs_to_plot':T2008_fig4_pfs})
pp_ma_CT_SP  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['ma_CT'], clr_map='clr_all_same', legend=False)
pp_ma_CT_pfs = ahf.Plot_Parameters(x_vars=['ma_CT'], y_vars=['press'], plot_type='profiles', clr_map='clr_all_same', extra_args={'pfs_to_plot':T2008_fig4_pfs})
#
pp_CT_and_ma_CT_pfs = ahf.Plot_Parameters(x_vars=['CT','ma_CT'], y_vars=['press'], plot_type='profiles', clr_map='clr_all_same', extra_args={'pfs_to_plot':[185]})

## Finding average profile
pp_avg_pf = ahf.Plot_Parameters(x_vars=['CT','SP'], y_vars=['press'])

### Figures for paper
## Map of ITP drifts
pp_ITP_map = ahf.Plot_Parameters(plot_type='map', clr_map='clr_by_instrmt')
pp_ITP_map_full_Arctic = ahf.Plot_Parameters(plot_type='map', clr_map='clr_by_instrmt', extra_args={'map_extent':'Full_Arctic'}, legend=False)

## Parameter sweeps
# ITP 2, Timmermans et al. 2008
pp_ITP2_ps_m_pts = ahf.Plot_Parameters(x_vars=['m_pts'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts, 'cl_ps_tuple':[20,455,10]})#, 'z_var':'maw_size', 'z_list':[60,100,150]})
pp_ITP2_ps_l_maw  = ahf.Plot_Parameters(x_vars=['maw_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts, 'cl_ps_tuple':[10,271,10]})#, 'z_var':'m_pts', 'z_list':[140,170,200]})
pp_ITP2_ps_n_pfs  = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts, 'cl_ps_tuple':[10,210,5], 'z_var':'m_pts', 'z_list':[70,100,120]})
# ITP 3
pp_ITP3_ps_m_pts = ahf.Plot_Parameters(x_vars=['m_pts'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts, 'cl_ps_tuple':[800,1001,10], 'z_var':'maw_size', 'z_list':[100]})
pp_ITP3_ps_l_maw  = ahf.Plot_Parameters(x_vars=['maw_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts, 'cl_ps_tuple':[10,210,10]})
pp_ITP3_ps_n_pfs  = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts, 'cl_ps_tuple':[50,770,50]})

## Evaluating clusterings with the overlap ratio and lateral density ratio
# Choose the value of m_pts depending on which ITP to evaluate the clustering of
eval_m_pts = T2008_m_pts
eval_m_pts = Lu2022_m_pts
pp_press_nir   = ahf.Plot_Parameters(x_vars=['nir_press'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':eval_m_pts, 'b_a_w_plt':False, 'plt_noise':False}, legend=True)
pp_temp_nir   = ahf.Plot_Parameters(x_vars=['nir_CT'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':eval_m_pts, 'b_a_w_plt':False, 'plt_noise':False}, legend=True)
pp_salt_nir   = ahf.Plot_Parameters(x_vars=['nir_SP'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':eval_m_pts, 'b_a_w_plt':False, 'plt_noise':False}, legend=False)
pp_salt_R_L = ahf.Plot_Parameters(x_vars=['cRL'], y_vars=['ca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':eval_m_pts, 'plot_slopes':True}, legend=False)

## Reproducing figures from Timmermans et al. 2008
## The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_T2008_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, ax_lims={'x_lims':ITP2_S_range}, legend=False)
# pp_T2008_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts},  legend=True)
# pp_T2008_clstr2 = ahf.Plot_Parameters(x_vars=['v1_SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'v1_SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=True)
## Reproducing Timmermans et al. 2008 Figure 4, with cluster coloring and 2 extra profiles
pp_T2008_fig4  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':T2008_fig4_pfs, 'plt_noise':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=True, ax_lims=T2008_fig4_y_lims)
## Reproducing Timmermans et al. 2008 Figure 5a, but with cluster coloring
pp_T2008_fig5a = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':False, 'isopycnals':0, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, ax_lims=T2008_fig5a_ax_lims, legend=False)
## Reproducing Timmermans et al. 2008 Figure 6a, but with cluster coloring
pp_T2008_fig6a = ahf.Plot_Parameters(x_vars=['BSP'], y_vars=['aCT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'isopycnals':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, ax_lims=T2008_fig6a_ax_lims, legend=False)

## Tracking clusters across a subset of profiles
pp_ITP2_some_pfs_0  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':ITP2_some_pfs_0, 'plt_noise':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=False, ax_lims=ITP2_some_pfs_ax_lims_0)
pp_ITP2_some_pfs_1  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':ITP2_some_pfs_1, 'plt_noise':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=False, ax_lims=ITP2_some_pfs_ax_lims_1)
pp_ITP2_T2008_fig4 = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='clr_all_same', extra_args={'pfs_to_plot':T2008_fig4_pfs, 'plt_noise':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=True, ax_lims=T2008_fig4_y_lims)

## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
pp_Lu2022_fig3a = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts}, legend=False)
pp_Lu2022_fig3b = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts}, legend=False)
pp_Lu2022_fig3c = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_SP'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts}, legend=False)

################################################################################
# Create analysis group objects
print('- Creating analysis group objects')
################################################################################

## Test Analysis Groups
# my_group0 = ahf.Analysis_Group(ds_ITP13_all, pfs_ITP13, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_test0)
# my_group1 = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_test1)
# my_group0 = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP3_all, pfs_0, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP1_all, pfs_0, pp_map)

# Finding the typical gradients of example profile
# pf_test0 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_1, pp_pfs0, plot_title='')
# pf_test1 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_1, pp_pfs1, plot_title='')
# pf_test2 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_1, pp_pfs2, plot_title='')
# pf_test3 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_1, pp_pfs3, plot_title='')

## Test figures
## Viewing example profiles with original, moving average, and local anomaly 
# group_og_ma_la_pf = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_og_ma_la_pf)

## Testing moving average window size
# group_maw_010 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_10, pp_test0, plot_title=r'ITP2 $\ell=10$ dbar')
# group_maw_050 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_50, pp_test0, plot_title=r'ITP2 $\ell=50$ dbar')
# group_maw_100 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100,pp_test0,plot_title=r'ITP2 $\ell=100$ dbar')
# group_maw_150 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_150,pp_test0,plot_title=r'ITP2 $\ell=150$ dbar')

## Testing impact of maw_size
# group_CT_SP  = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100, pp_CT_SP)
# group_CT_pfs = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100, pp_CT_pfs)
# group_ma_CT_SP  = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100, pp_ma_CT_SP)
# group_ma_CT_pfs = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100, pp_ma_CT_pfs)
# # 
# group_CT_and_ma_CT_pfs__10 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_10, pp_CT_and_ma_CT_pfs, plot_title=r'ITP2 $\ell=10$ dbar')
# group_CT_and_ma_CT_pfs__50 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_50, pp_CT_and_ma_CT_pfs, plot_title=r'ITP2 $\ell=50$ dbar')
# group_CT_and_ma_CT_pfs_100 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100, pp_CT_and_ma_CT_pfs, plot_title=r'ITP2 $\ell=100$ dbar')
# group_CT_and_ma_CT_pfs_150 = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_150, pp_CT_and_ma_CT_pfs, plot_title=r'ITP2 $\ell=150$ dbar')

## Finding average profile
# group_avg_pf = ahf.Analysis_Group(ds_ITP2_all, pfs_maw_100, pp_avg_pf)

### Figures for paper

## Map of ITP drifts
# group_ITP_map = ahf.Analysis_Group(ds_all_ITPs, pfs_0, pp_ITP_map, plot_title='Profile locations')
# group_ITP_map_full_Arctic = ahf.Analysis_Group(ds_all_ITPs, pfs_0, pp_ITP_map_full_Arctic, plot_title='Profile locations')

## Parameter sweeps
## ITP2, Timmermans et al. 2008
# group_ps_l_maw = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_ITP2_ps_l_maw)
# group_ps_m_pts = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_ITP2_ps_m_pts)
# group_ps_n_pfs = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_ITP2_ps_n_pfs)
# ITP3, Lu et al. 2022
# group_ps_l_maw = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_ITP3_ps_l_maw)
# group_ps_m_pts = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_ITP3_ps_m_pts)
# group_ps_n_pfs = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_ITP3_ps_n_pfs)

## Evaluating clusterings with the overlap ratio and lateral density ratio
# For the reproduction of Timmermans et al. 2008
# group_clstr_ST = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_clstr_ST)
# group_salt_hist= ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_salt_hist)
# group_press_nir = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_press_nir)
# group_temp_nir = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_temp_nir)
# group_salt_nir = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_salt_nir)
# group_salt_R_L = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_salt_R_L)
# For the reproduction of Lu et al. 2022
# group_clstr_ST = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_clstr_ST)
# group_salt_hist= ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_salt_hist)
# group_press_nir = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_press_nir)
# group_temp_nir = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_temp_nir)
# group_salt_nir = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_salt_nir)
# group_salt_R_L = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_salt_R_L)

## Reproducing figures from Timmermans et al. 2008
# group_T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_clstr, plot_title='')
# group_T2008_clstr2 = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_clstr2, plot_title='stretched')
# group_T2008_fig4  = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_fig4, plot_title='')
group_T2008_fig5a = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_fig5a, plot_title='')
# group_T2008_fig6a = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_T2008_fig6a, plot_title='')

## Tracking clusters across a subset of profiles
# group_ITP2_some_pfs_0 = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_ITP2_some_pfs_0)
# group_ITP2_some_pfs_1 = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_ITP2_some_pfs_1)
# group_ITP2_T2008_pfs = ahf.Analysis_Group(ds_ITP2_all, pfs_T2008, pp_ITP2_T2008_fig4)

## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
# group_Lu2022_fig3a = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3a)
# group_Lu2022_fig3b = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3b)
# group_Lu2022_fig3c = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3c)
#
# group_Lu2022_clstr = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_clstr)

################################################################################
# Declare figures or summaries to output
print('- Creating outputs')
################################################################################

# ahf.make_figure([group_T2008_fig4])#, filename='test.pickle')
# ahf.make_figure([group_salt_nir])
# ahf.make_figure([group_ps_n_pfs])
# ahf.make_figure([my_group0])#, filename='test.pickle')
# ahf.make_figure([my_group0, my_group1])
# ahf.make_figure([group_T2008_clstr])

# ahf.make_figure([pf_test0, pf_test1, pf_test2, pf_test3])

## Test Figures
# ahf.make_figure([group_press_hist, group_press_nir, group_sigma_hist, group_sigma_nir, group_temp_hist, group_temp_nir, group_salt_hist, group_salt_nir], filename='ITP2_nir_vs_press_all_var.pickle')
# ahf.make_figure([group_clstr_ST, group_salt_nir, group_salt_R_L])#, filename='ITP2_nir_vs_press_all_var.pickle')
# ahf.make_figure([group_salt_R_L, group_salt_R_l])#, filename='ITP2_nir_vs_press_all_var.pickle')
# ahf.make_figure([group_T2008_clstr, group_salt_hist, group_salt_nir, group_salt_R_L])#, filename='ITP2_nir_vs_press_all_var.pickle')
# ahf.make_figure([group_maw_001, group_maw_005, group_maw_010, group_maw_050, group_maw_100, group_maw_200], filename='ITP2_maw_size_tests.png')

## Testing impact of maw_size
# ahf.make_figure([group_maw_010, group_maw_050, group_maw_100, group_maw_150])
# ahf.make_figure([group_CT_SP, group_ma_CT_SP, group_CT_pfs, group_ma_CT_pfs])
# ahf.make_figure([group_CT_and_ma_CT_pfs__10, group_CT_and_ma_CT_pfs__50, group_CT_and_ma_CT_pfs_100, group_CT_and_ma_CT_pfs_150])
# ahf.make_figure([group_CT_and_ma_CT_pfs_100])

## Finding average profile
# ahf.find_avg_profile(group_avg_pf, 'press', 0.25)

### Figures for paper

## Map of ITP drifts
# ahf.make_figure([group_ITP_map], use_same_x_axis=False, use_same_y_axis=False)#, filename='ITP_map.pickle')
# ahf.make_figure([group_ITP_map_full_Arctic], use_same_x_axis=False, use_same_y_axis=False)
# ahf.make_figure([group_ITP_map_full_Arctic, group_ITP_map], use_same_x_axis=False, use_same_y_axis=False)#, filename='ITP_map.pickle')
# ahf.find_max_distance([group_ITP_map])

## Parameter sweeps
# ahf.make_figure([group_ps_l_maw, group_ps_m_pts], use_same_y_axis=False)#, filename='ITP2_param_sweep.pickle')
# ahf.make_figure([group_ps_m_pts], filename='ITP2_test_sweep.pickle')

## Evaluating clusterings with the overlap ratio and lateral density ratio
# ahf.make_figure([group_salt_R_L])
# ahf.make_figure([group_salt_nir, group_salt_R_L], filename='nir_R_L_vs_SP.pickle')
# ahf.make_figure([group_salt_nir])#, filename='nir_vs_press.pickle')
# ahf.make_figure([group_press_nir, group_temp_nir], filename='press_and_temp_nir_ITP2.pickle')

## Reproducing figures from Timmermans et al. 2008
# ahf.make_figure([group_T2008_fig5a, group_T2008_fig4, group_T2008_clstr, group_T2008_fig6a], filename='T2008.png')
# ahf.make_figure([group_T2008_clstr])
# ahf.make_figure([group_T2008_clstr, group_T2008_clstr2], use_same_y_axis=False, use_same_x_axis=False)
# ahf.make_figure([group_T2008_fig4])
ahf.make_figure([group_T2008_fig5a])
# ahf.make_figure([group_T2008_fig6a])

## Tracking clusters across a subset of profiles
# ahf.make_figure([group_ITP2_some_pfs_0], use_same_x_axis=False, use_same_y_axis=False)
# ahf.make_figure([group_ITP2_some_pfs_1], use_same_x_axis=False, use_same_y_axis=False)
# ahf.make_figure([group_ITP2_T2008_pfs], use_same_x_axis=False, use_same_y_axis=False)

## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
# ahf.make_figure([group_Lu2022_fig3a, group_Lu2022_fig3b, group_Lu2022_fig3c], filename='Lu2022_f3.pickle')
# ahf.make_figure([group_Lu2022_clstr])
