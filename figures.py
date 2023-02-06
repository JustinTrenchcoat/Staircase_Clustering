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

# Depth filters
ITP2_p_range = [185,300]
ITP3_p_range = [180,300]
ITP13_p_range = [195,300]

### Filters for reproducing plots from Timmermans et al. 2008
# Timmermans 2008 Figure 4 depth range
T2008_fig4_y_lims = {'y_lims':[260,220]}
# Timmermans 2008 Figure 4 shows profile 185
T2008_fig4_pfs = [183, 185, 187]
# Filters used in Timmermans 2008 T-S and aT-BS plots
T2008_p_range = [180,300]
T2008_S_range = [34.4,34.6]
T2008_fig5a_ax_lims = {'x_lims':[34.05,34.75], 'y_lims':[-1.3,0.5]}
T2008_fig6a_ax_lims = {'x_lims':[0.027,0.027045], 'y_lims':[-13e-6,3e-6]}
# The actual limits are above, but need to adjust the x lims for some reason
T2008_fig6a_ax_lims = {'x_lims':[0.026835,0.026880], 'y_lims':[-13e-6,3e-6]}

# Filters used in Bebieva et al. 2019
B2019_p_range = [220,400]

# Shibley et al. 2019 Figure 6a I think shows profile 655
S2019_fig6a_pfs = [654, 655, 656]
# Filters to reproduce Shibley et al. 2019 Figure 6b
S2019_p_range = [180,300]

### Filters for reproducing prlots from Lu et al. 2022
Lu2022_p_range = [170,410]
Lu2022_t_range = [-1.5,1]
Lu2022_s_range = [34.2,34.9]

################################################################################
# Make dictionaries for what data to load in and analyze
################################################################################

# All profiles from all ITPs in this study
all_ITPs = {'ITP_1':'all','ITP_2':'all','ITP_3':'all','ITP_13':'all'}

# All profiles from certain ITPs
ITP1_all  = {'ITP_1':'all'}
ITP2_all  = {'ITP_2':'all'}
ITP3_all  = {'ITP_3':'all'}
ITP13_all = {'ITP_13':'all'}

# Just specific profiles
ITP2_pfs  = {'ITP_2':T2008_fig4_pfs}
ITP13_pfs = {'ITP_13':S2019_fig6a_pfs}

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

ds_ITP1_all = ahf.Data_Set(ITP1_all, dfs0)

ds_ITP2_all = ahf.Data_Set(ITP2_all, dfs0)
ds_ITP2_pfs = ahf.Data_Set(ITP2_pfs, dfs0)

ds_ITP3_all = ahf.Data_Set(ITP3_all, dfs0)
# ds_ITP3_all = ahf.Data_Set(ITP3_all, dfs_test)

ds_ITP13_all = ahf.Data_Set(ITP13_all, dfs0)#dfs_S2019)
ds_ITP13_pfs = ahf.Data_Set(ITP13_pfs, dfs0)

################################################################################
# Create profile filtering objects
print('- Creating profile filtering objects')
################################################################################

pfs_f0 = ahf.Profile_Filters()
pfs_ITP2  = ahf.Profile_Filters(p_range=ITP2_p_range)#, m_avg_win=50)
pfs_ITP3  = ahf.Profile_Filters(p_range=ITP3_p_range)
pfs_Lu2022= ahf.Profile_Filters(p_range=Lu2022_p_range, CT_range=Lu2022_t_range, SP_range=Lu2022_s_range)
pfs_ITP13 = ahf.Profile_Filters(p_range=ITP13_p_range, SP_range=[34.1,34.8])
pfs_B2019 = ahf.Profile_Filters(p_range=B2019_p_range)

################################################################################
# Create plotting parameter objects

print('- Creating plotting parameter objects')
################################################################################

### Test plots
pp_xy_default = ahf.Plot_Parameters()
pp_test0 = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['cmc_press'], clr_map='cluster', legend=True, extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90})
# pp_test1 = ahf.Plot_Parameters(x_vars=['ca_SP'], y_vars=['cmm_press'], clr_map='cluster', legend=True, extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90})
# pp_test1 = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['hist'], clr_map='cluster', legend=True, extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90})
# pp_test1 = ahf.Plot_Parameters(x_vars=['pcs_press'], y_vars=['pca_press'], clr_map='density_hist', legend=True, extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'clr_min':0, 'clr_max':15, 'clr_ext':'max', 'xy_bins':250})
# pp_test0 = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['lat'], clr_map='clr_all_same')
pp_map = ahf.Plot_Parameters(plot_type='map', clr_map='prof_no')
pp_pfs = ahf.Plot_Parameters(x_vars=['CT'], y_vars=['press'], plot_type='profiles')

## Test Figures
#
pp_cmm_SP = ahf.Plot_Parameters(x_vars=['cmm_SP'], y_vars=['ca_press'], clr_map='cluster', legend=True, extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90})
pp_cmm_press = ahf.Plot_Parameters(x_vars=['ca_SP'], y_vars=['cmm_press'], clr_map='cluster', legend=True, extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90})
#
pp_press_cor  = ahf.Plot_Parameters(x_vars=['cor_press'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_press_hist = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_sigma_cor  = ahf.Plot_Parameters(x_vars=['cor_sigma'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_sigma_hist = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['sigma'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_temp_cor   = ahf.Plot_Parameters(x_vars=['cor_CT'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_temp_hist  = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['CT'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_salt_cor   = ahf.Plot_Parameters(x_vars=['cor_SP'], y_vars=['ca_press'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})
pp_salt_hist  = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['SP'], clr_map='cluster', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'b_a_w_plt':False, 'plt_noise':False})



### Figures for paper
## Map of ITP drifts
pp_ITP_map = ahf.Plot_Parameters(plot_type='map', clr_map='clr_by_instrmt')

## Parameter sweeps
# ITP 2, Timmermans et al. 2008
pp_ITP2_ps_min_cs = ahf.Plot_Parameters(x_vars=['min_cs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'cl_ps_tuple':[20,450,10], 'z_var':'maw_size', 'z_list':[50,90,110]})
pp_ITP2_ps_l_maw  = ahf.Plot_Parameters(x_vars=['maw_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'cl_ps_tuple':[10,310,10], 'z_var':'min_cs', 'z_list':[50,90,140]})
pp_ITP2_ps_n_pfs  = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'cl_ps_tuple':[10,210,5], 'z_var':'min_cs', 'z_list':[50,90,140]})
# ITP 2, Bebieva et al. 2019
pp_ITP2B_ps_min_cs = ahf.Plot_Parameters(x_vars=['min_cs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'cl_ps_tuple':[20,550,10], 'z_var':'maw_size', 'z_list':[50,60,140]})
pp_ITP2B_ps_l_maw  = ahf.Plot_Parameters(x_vars=['maw_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'cl_ps_tuple':[10,310,10], 'z_var':'min_cs', 'z_list':[130,210,340]})
pp_ITP2B_ps_n_pfs  = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'cl_ps_tuple':[10,210,5], 'z_var':'min_cs', 'z_list':[130,210,340]})
# ITP 3
pp_ITP3_ps_min_cs = ahf.Plot_Parameters(x_vars=['min_cs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':500, 'cl_ps_tuple':[500,50001,500]})
pp_ITP3_ps_l_maw  = ahf.Plot_Parameters(x_vars=['maw_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':500, 'cl_ps_tuple':[10,210,10]})
pp_ITP3_ps_n_pfs  = ahf.Plot_Parameters(x_vars=['n_pfs'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':500, 'cl_ps_tuple':[50,770,50]})

## Reproducing figures from Timmermans et al. 2008
## The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_T2008_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, legend=False)
## Reproducing Timmermans et al. 2008 Figure 4, with cluster coloring and 2 extra profiles
pp_T2008_fig4  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':T2008_fig4_pfs, 'plt_noise':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, ax_lims=T2008_fig4_y_lims)
## Reproducing Timmermans et al. 2008 Figure 5a, but with cluster coloring
pp_T2008_fig5a = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='cluster', extra_args={'b_a_w_plt':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, ax_lims=T2008_fig5a_ax_lims, legend=True)
## Reproducing Timmermans et al. 2008 Figure 6a, but with cluster coloring
pp_T2008_fig6a = ahf.Plot_Parameters(x_vars=['BSP'], y_vars=['aCT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, ax_lims=T2008_fig6a_ax_lims, legend=False)

## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
ITP3_min_cs = 90#300
pp_Lu2022_fig3a = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':ITP3_min_cs}, legend=False)
pp_Lu2022_fig3b = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':ITP3_min_cs}, legend=False)
pp_Lu2022_fig3c = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_SP'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':ITP3_min_cs}, legend=False)
pp_Lu2022_fig3d = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_sigma'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':ITP3_min_cs}, legend=False)

## The actual clustering done for reproducing Figure 3 from Bebieva et al. 2019
pp_B2019_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, legend=False)
## Reproducing Bebieva et al. 2019 Figure 3b
pp_B2019_fig3b = ahf.Plot_Parameters(x_vars=['cRL'], y_vars=['ca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90}, legend=False)

## Finding layer height vs. depth, reproducing Shibley et al. 2019 Figure 6
pp_S2019_fig6a = ahf.Plot_Parameters(x_vars=['CT'], y_vars=['press'], plot_type='profiles')
pp_S2019_fig6a = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':S2019_fig6a_pfs, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':230})
pp_S2019_fig6c = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':True, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':230}, legend=True)
pp_S2019_fig6b = ahf.Plot_Parameters(x_vars=['pcs_press'], y_vars=['press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':230}, legend=True, ax_lims={'x_lims':[0,10]})

## Histograms of data that's been mean centered by cluster
pp_cmc_press = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['cmc_press'], extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'plt_hist_lines':True})
pp_cmc_sigma = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['cmc_sigma'], extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'plt_hist_lines':True})
pp_cmc_temp  = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['cmc_CT'], extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'plt_hist_lines':True})
pp_cmc_salt  = ahf.Plot_Parameters(x_vars=['hist'], y_vars=['cmc_SP'], extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'min_cs':90, 'plt_hist_lines':True})

################################################################################
# Create analysis group objects
print('- Creating analysis group objects')
################################################################################

## Test Analysis Groups
# my_group0 = ahf.Analysis_Group(ds_ITP13_all, pfs_ITP13, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_test0)
# my_group1 = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_test1)
# my_group0 = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP3_all, pfs_f0, pp_test0)
# my_group0 = ahf.Analysis_Group(ds_ITP1_all, pfs_f0, pp_map)
# my_group1 = ahf.Analysis_Group(ds_ITP2_pfs, pfs_ITP2, pp_test0)

## Test figures
#
# group_cmm_SP = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_cmm_SP)
# group_cmm_press = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_cmm_press)
# Does the overlap ratio change with density? According to these plots, no
# group_press_cor  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_press_cor)
# group_press_hist = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_press_hist)
# group_sigma_cor  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_sigma_cor)
# group_sigma_hist = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_sigma_hist)
# group_temp_cor   = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_temp_cor)
# group_temp_hist  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_temp_hist)
# group_salt_cor   = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_salt_cor)
# group_salt_hist  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_salt_hist)

### Figures for paper

## Map of ITP drifts
# group_ITP_map = ahf.Analysis_Group(ds_all_ITPs, pfs_f0, pp_ITP_map)

## Parameter sweeps
## ITP2, Timmermans et al. 2008
# group_ps_min_cs = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_ps_min_cs)
# group_ps_l_maw  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_ps_l_maw)
# group_ps_n_pfs  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_ps_n_pfs)
## ITP2, Bebieva et al. 2019
group_ps_min_cs = ahf.Analysis_Group(ds_ITP2_all, pfs_B2019, pp_ITP2B_ps_min_cs)
group_ps_l_maw  = ahf.Analysis_Group(ds_ITP2_all, pfs_B2019, pp_ITP2B_ps_l_maw)
group_ps_n_pfs  = ahf.Analysis_Group(ds_ITP2_all, pfs_B2019, pp_ITP2B_ps_n_pfs)
#
# group_ps_min_cs = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_ITP3_ps_min_cs)
# group_ps_l_maw  = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_ITP3_ps_l_maw)
# group_ps_n_pfs  = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_ITP3_ps_n_pfs)

## Reproducing figures from Timmermans et al. 2008
# group_T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_clstr)
# group_T2008_fig4  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_fig4)
# group_T2008_fig5a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_fig5a)
# group_T2008_fig6a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_fig6a)

## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
# group_Lu2022_fig3a = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3a)
# group_Lu2022_fig3b = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3b)
# group_Lu2022_fig3c = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3c)
# group_Lu2022_fig3d = ahf.Analysis_Group(ds_ITP3_all, pfs_Lu2022, pp_Lu2022_fig3d)
# group_Lu2022_fig3a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3a)
# group_Lu2022_fig3b = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3b)
# group_Lu2022_fig3c = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3c)
# group_Lu2022_fig3d = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3d)

## The actual clustering done for reproducing Figure 3 from Bebieva et al. 2019
# group_B2019_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_B2019, pp_B2019_clstr)
## Reproducing Bebieva et al. 2019 Figure 3b
# group_B2019_fig3b = ahf.Analysis_Group(ds_ITP2_all, pfs_B2019, pp_B2019_fig3b)

## Finding layer height vs. depth, reproducing Shibley et al. 2019 Figure 6
# group_S2019_fig6a = ahf.Analysis_Group(ds_ITP13_all, pfs_ITP13, pp_S2019_fig6a)
# group_S2019_fig6b = ahf.Analysis_Group(ds_ITP13_all, pfs_ITP13, pp_S2019_fig6b)
# group_S2019_fig6c = ahf.Analysis_Group(ds_ITP13_all, pfs_ITP13, pp_S2019_fig6c)

## Histograms of data that's been mean centered by cluster
# group_cmc_press = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_cmc_press)
# group_cmc_sigma = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_cmc_sigma)
# group_cmc_temp  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_cmc_temp)
# group_cmc_salt  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_cmc_salt)

################################################################################
# Declare figures or summaries to output
print('- Creating outputs')
################################################################################

# ahf.make_figure([my_group0])
# ahf.make_figure([my_group0], filename='test.pickle')
# ahf.make_figure([my_group0, my_group1])
# ahf.make_figure([group_T2008_fig4])

## Test Figures
# ahf.make_figure([group_cmm_SP, group_cmm_press])
# ahf.make_figure([group_press_hist, group_press_cor, group_sigma_hist, group_sigma_cor, group_temp_hist, group_temp_cor, group_salt_hist, group_salt_cor], filename='ITP2_cor_vs_press_all_var.pickle')

### Figures for paper

## Map of ITP drifts
# ahf.make_figure([group_ITP_map])#, filename='ITP_map.pickle')

## Parameter sweeps
ahf.make_figure([group_ps_min_cs, group_ps_l_maw, group_ps_n_pfs], filename='ITP2B_sweep.pickle')

## Reproducing figures from Timmermans et al. 2008
# ahf.make_figure([group_T2008_clstr, group_T2008_fig4, group_T2008_fig5a, group_T2008_fig6a])#, filename='T2008.pickle')

## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
# ahf.make_figure([group_Lu2022_fig3a, group_Lu2022_fig3d, group_Lu2022_fig3b, group_Lu2022_fig3c], filename='Lu2022_f3_ITP2.png')#'Lu2022_f3.pickle')

## Reproducing Bebieva et al. 2019 Figure 3b
# ahf.make_figure([group_B2019_clstr, group_B2019_fig3b])

## Finding layer height vs. depth, reproducing Shibley et al. 2019 Figure 6
# ahf.make_figure([group_S2019_fig6a, group_S2019_fig6b, group_S2019_fig6c])

## Histograms of data that's been mean centered by cluster
# ahf.make_figure([group_cmc_press, group_cmc_sigma, group_cmc_temp, group_cmc_salt])
