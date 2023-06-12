"""
Author: Mikhail Schee
Created: 2022-08-18

This script is set up to make figures of Arctic Ocean profile data that has been
formatted into netcdfs by the `make_netcdf` function.

Note: You unfortunately cannot pickle plots with parasite axes. So, the pickling
functionality won't work if you try to make a plot with the box and whisker plot
in an inset.

The Jupyter notebook Create_Figures.ipynb offers detailed explanations of each 
plot created by this script.
"""

"""
NOTE: BEFORE YOU RUN THIS SCRIPT

This script expects that the following files exist:
netcdfs/ITP_2.nc
netcdfs/ITP_3.nc

and that they have been created by running the following scripts in this order:
make_netcdf.py
take_moving_average.py
cluster_data.py
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
T2008_fig5a_x_lims  = {'x_lims':[34.05,34.75]}
T2008_fig5a_ax_lims = {'x_lims':[34.05,34.75], 'y_lims':[-1.3,0.5]}
T2008_fig6a_ax_lims = {'x_lims':[0.027002,0.027042], 'y_lims':[-13e-6,3e-6]}
# The actual limits are above, but need to adjust the x lims for some reason
T2008_fig6a_ax_lims = {'x_lims':[0.026838,0.026878], 'y_lims':[-13e-6,3e-6]}

# A list of many profiles to plot from ITP2
start_pf = 1
import numpy as np
n_pfs_to_plot = 10
ITP2_some_pfs = list(np.arange(start_pf, start_pf+(n_pfs_to_plot*2), 2))

# For showing multiple layers grouped into one cluster
ITP2_some_pfs_0 = [87, 89, 95, 97, 99, 101, 103, 105, 109, 111]
ITP2_some_pfs_ax_lims_0 = {'y_lims':[245,220]}
# For showing one layer split into multiple clusters
ITP2_some_pfs_1 = [67, 69, 73, 75, 81, 83, 91, 93, 97, 99]
ITP2_some_pfs_ax_lims_1 = {'y_lims':[295,270]}

# For showing how well the clustering does near the AW core
ITP3_some_pfs_2 = [313, 315, 317, 319, 321]
ITP3_some_pfs_ax_lims_2 = {'y_lims':[370,300]}

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

# A list of many profiles to plot from ITP3
start_pf = 1331
import numpy as np
n_pfs_to_plot = 5
ITP3_some_pfs_1 = list(np.arange(start_pf, start_pf+(n_pfs_to_plot*2), 2))
ITP3_pfs1  = {'ITP_3':ITP3_some_pfs_1}
ITP3_pfs2  = {'ITP_3':ITP3_some_pfs_2}

################################################################################
# Create data filtering objects
print('- Creating data filtering objects')
################################################################################

dfs0 = ahf.Data_Filters()

################################################################################
# Create data sets by combining filters and the data to load in
print('- Creating data sets')
################################################################################

ds_all_ITPs = ahf.Data_Set(all_ITPs, dfs0)

ds_ITP2_all = ahf.Data_Set(ITP2_all, dfs0)

ds_ITP3_all = ahf.Data_Set(ITP3_all, dfs0)
ds_ITP3_some_pfs2 = ahf.Data_Set(ITP3_pfs2, dfs0)

################################################################################
# Create profile filtering objects
print('- Creating profile filtering objects')
################################################################################

pfs_0 = ahf.Profile_Filters()

pfs_ITP2  = ahf.Profile_Filters(SP_range=ITP2_S_range)
pfs_ITP3 = ahf.Profile_Filters(SP_range=Lu2022_S_range)

################################################################################
# Create plotting parameter objects

print('- Creating plotting parameter objects')
################################################################################

### Test plots
## Regrid
# Make the Profile Filters
pfs_regrid = ahf.Profile_Filters(SP_range=ITP2_S_range, regrid_TS=['CT',0.01,'SP',0.005])
# Make the Plot Parameters
pp_regrid = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='clr_all_same', legend=True)
# Make the Analysis Group
group_regrid = ahf.Analysis_Group(ds_ITP2_all, pfs_regrid, pp_regrid)
# Make the figure
ahf.make_figure([group_regrid])

### Figure 2, but with S_A instead of S_P
## Using ITP2, reproducing figures from Timmermans et al. 2008
# The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_2T2008_clstr = ahf.Plot_Parameters(x_vars=['SA'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'SA', 'cl_y_var':'la_CT', 'm_pts':170}, ax_lims={'x_lims':[34.2169,34.9171]}, legend=True, add_grid=False)
# Reproducing Timmermans et al. 2008 Figure 4, with cluster coloring and 2 extra profiles
pp_2T2008_fig4  = ahf.Plot_Parameters(x_vars=['SA'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':T2008_fig4_pfs, 'plt_noise':True, 'cl_x_var':'SA', 'cl_y_var':'la_CT', 'm_pts':170}, legend=True, ax_lims=T2008_fig4_y_lims, add_grid=True)
# Reproducing Timmermans et al. 2008 Figure 5a, but with cluster coloring
pp_2T2008_fig5a = ahf.Plot_Parameters(x_vars=['SA'], y_vars=['CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':False, 'isopycnals':0, 'place_isos':'manual', 'cl_x_var':'SA', 'cl_y_var':'la_CT', 'm_pts':170}, ax_lims={'x_lims':[34.2169,34.9171], 'y_lims':[-1.3,0.5]}, legend=True, add_grid=False)
# Reproducing Timmermans et al. 2008 Figure 6a, but with cluster coloring
pp_2T2008_fig6a = ahf.Plot_Parameters(x_vars=['BSA'], y_vars=['aCT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'isopycnals':True, 'cl_x_var':'SA', 'cl_y_var':'la_CT', 'm_pts':170}, ax_lims={'x_lims':[0.026967,0.027007], 'y_lims':[-13e-6,3e-6]}, legend=True, add_grid=False)

### Figure 4, but with S_A
## Evaluating clusterings with the overlap ratio and lateral density ratio
pp_2salt_nir = ahf.Plot_Parameters(x_vars=['nir_SA'], y_vars=['ca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plt_noise':False, 'cl_x_var':'SA', 'cl_y_var':'la_CT', 'm_pts':170}, legend=False)
pp_2salt_R_L = ahf.Plot_Parameters(x_vars=['cRL'], y_vars=['ca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plt_noise':False, 'plot_slopes':True, 'cl_x_var':'SA', 'cl_y_var':'la_CT', 'm_pts':170}, legend=False)

## Reproducing figures from Timmermans et al. 2008
## The actual clustering done for reproducing figures from Timmermans et al. 2008
# pp_T2008_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False},  legend=False)
# pp_T2008_clstr2 = ahf.Plot_Parameters(x_vars=['v1_SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'cl_x_var':'v1_SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts}, legend=True)

################################################################################
### Figures for paper

### Figure 1
## Map of ITP drifts for ITP2 and ITP3
pp_ITP_map = ahf.Plot_Parameters(plot_type='map', clr_map='clr_by_instrmt')
pp_ITP_map_full_Arctic = ahf.Plot_Parameters(plot_type='map', clr_map='clr_by_instrmt', extra_args={'map_extent':'Full_Arctic'}, legend=False, add_grid=False)

### Figure 2
## Using ITP2, reproducing figures from Timmermans et al. 2008
# The actual clustering done for reproducing figures from Timmermans et al. 2008
pp_T2008_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False}, ax_lims=T2008_fig5a_x_lims, legend=False, add_grid=False)
# Reproducing Timmermans et al. 2008 Figure 4, with cluster coloring and 2 extra profiles
pp_T2008_fig4  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':T2008_fig4_pfs, 'plt_noise':True}, legend=True, ax_lims=T2008_fig4_y_lims, add_grid=True)
# Reproducing Timmermans et al. 2008 Figure 5a, but with cluster coloring
pp_T2008_fig5a = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':False, 'isopycnals':0, 'place_isos':'auto'}, ax_lims=T2008_fig5a_ax_lims, legend=False, add_grid=False)
# Reproducing Timmermans et al. 2008 Figure 6a, but with cluster coloring
pp_T2008_fig6a = ahf.Plot_Parameters(x_vars=['BSP'], y_vars=['aCT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':True, 'isopycnals':True}, ax_lims=T2008_fig6a_ax_lims, legend=False, add_grid=False)
## For ITP3
# The actual clustering done for reproducing figures from Lu et al. 2022
pp_Lu2022_clstr = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['la_CT'], clr_map='cluster', extra_args={'b_a_w_plt':True}, legend=True, add_grid=False)
# Plotting clusters back in regular CT vs SP space
pp_Lu2022_CT_SP = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['CT'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plot_slopes':False, 'isopycnals':0, 'place_isos':'auto'}, legend=False, add_grid=False)

### Figure 3
## Parameter sweep across \ell and m_pts for ITP2, Timmermans et al. 2008
pp_ITP2_ps_m_pts = ahf.Plot_Parameters(x_vars=['m_pts'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts, 'cl_ps_tuple':[20,455,10]}, legend=False, add_grid=False)
pp_ITP2_ps_ell  = ahf.Plot_Parameters(x_vars=['ell_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':T2008_m_pts, 'cl_ps_tuple':[10,271,10]}, legend=False, add_grid=False)
## Parameter sweep across \ell and m_pts for ITP3, Lu et al. 2022
pp_ITP3_ps_m_pts = ahf.Plot_Parameters(x_vars=['m_pts'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts, 'cl_ps_tuple':[800,1001,10]}, legend=False, add_grid=False)
pp_ITP3_ps_ell  = ahf.Plot_Parameters(x_vars=['ell_size'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'cl_x_var':'SP', 'cl_y_var':'la_CT', 'm_pts':Lu2022_m_pts, 'cl_ps_tuple':[10,210,10]}, legend=False, add_grid=False)

### Figure 4
## Evaluating clusterings with the lateral density ratio and the normalized inter-cluster range
pp_salt_R_L = ahf.Plot_Parameters(x_vars=['cRL'], y_vars=['ca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plt_noise':False, 'plot_slopes':True}, legend=False)
pp_salt_nir = ahf.Plot_Parameters(x_vars=['nir_SP'], y_vars=['ca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False, 'plt_noise':False}, legend=False)

### Figure 5
## Tracking clusters across a subset of profiles
# For ITP2
pp_ITP2_some_pfs_0  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':ITP2_some_pfs_0, 'plt_noise':True}, legend=False, ax_lims=ITP2_some_pfs_ax_lims_0)
pp_ITP2_some_pfs_1  = ahf.Plot_Parameters(x_vars=['SP'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'pfs_to_plot':ITP2_some_pfs_1, 'plt_noise':True}, legend=False, ax_lims=ITP2_some_pfs_ax_lims_1)

### Figure 6
## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
pp_Lu2022_fig3a = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_press'], clr_map='cluster', extra_args={'b_a_w_plt':False}, legend=False)
pp_Lu2022_fig3b = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_CT'], clr_map='cluster', extra_args={'b_a_w_plt':False}, legend=False)
pp_Lu2022_fig3c = ahf.Plot_Parameters(x_vars=['dt_start'], y_vars=['pca_SP'], clr_map='cluster', extra_args={'b_a_w_plt':False}, legend=False)

### Figure 8
## Tracking clusters across a subset of profiles
# For ITP3
pp_ITP3_some_pfs_2 = ahf.Plot_Parameters(x_vars=['SP','CT'], y_vars=['press'], plot_type='profiles', clr_map='cluster', extra_args={'plt_noise':True}, legend=False, ax_lims=ITP3_some_pfs_ax_lims_2)

################################################################################
# Create analysis group objects
print('- Creating analysis group objects')
################################################################################

## Test figures

### Figure 2, but using S_A
## Using ITP2, reproducing figures from Timmermans et al. 2008
if False:
    print('')
    print('- Creating Figure 2')
    # Make the subplot groups
    group_2T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_2T2008_clstr, plot_title='')
    group_2T2008_fig4  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_2T2008_fig4, plot_title='')
    group_2T2008_fig5a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_2T2008_fig5a, plot_title='')
    group_2T2008_fig6a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_2T2008_fig6a, plot_title='')
    # Make the figure
    #   Remember, adding isopycnals means it will prompt you to place the in-line labels manually
    ahf.make_figure([group_2T2008_fig5a, group_2T2008_fig4, group_2T2008_clstr, group_2T2008_fig6a])#, filename='Figure_2_SA.png')
    # ahf.make_figure([group_2T2008_clstr])
    # ahf.make_figure([group_2T2008_fig4])
    # ahf.make_figure([group_2T2008_fig5a])
    # ahf.make_figure([group_2T2008_fig6a])
### Figure 4, but using S_A
## Evaluating clusterings with the overlap ratio and lateral density ratio
# For the reproduction of Timmermans et al. 2008
if False:
    print('')
    print('- Creating Figure 4')
    # Make the subplot groups
    group_2salt_nir = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_2salt_nir, plot_title='')
    group_2salt_R_L = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_2salt_R_L, plot_title='')
    # Make the figure
    ahf.make_figure([group_2salt_nir, group_2salt_R_L], filename='Figure_4_SA.png')

################################################################################
### Figures for paper

### Figure 1
## Map of ITP drifts for ITP2 and ITP3
if False:
    print('')
    print('- Creating Figure 1')
    # Make the subplot groups
    group_ITP_map = ahf.Analysis_Group(ds_all_ITPs, pfs_0, pp_ITP_map, plot_title='')
    group_ITP_map_full_Arctic = ahf.Analysis_Group(ds_all_ITPs, pfs_0, pp_ITP_map_full_Arctic, plot_title='')
    # Make the figure
    ahf.make_figure([group_ITP_map_full_Arctic, group_ITP_map], use_same_x_axis=False, use_same_y_axis=False, filename='Figure_1.pickle')
    # Find the maximum distance between any two profiles for each data set in the group
    # ahf.find_max_distance([group_ITP_map])

### Figure 2
## Using ITP2, reproducing figures from Timmermans et al. 2008
if False:
    print('')
    print('- Creating Figure 2')
    # Make the subplot groups
    group_T2008_clstr = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_clstr, plot_title='')
    group_T2008_fig4  = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_fig4, plot_title='')
    group_T2008_fig5a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_fig5a, plot_title='')
    group_T2008_fig6a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_T2008_fig6a, plot_title='')
    # Make the figure
    #   Remember, adding isopycnals means it will prompt you to place the in-line labels manually
    ahf.make_figure([group_T2008_fig5a, group_T2008_fig4, group_T2008_clstr, group_T2008_fig6a])#, filename='Figure_2.png')
    # ahf.make_figure([group_T2008_clstr])
    # ahf.make_figure([group_T2008_fig4])
    # ahf.make_figure([group_T2008_fig5a])
    # ahf.make_figure([group_T2008_fig6a])
## The actual clustering done for reproducing figures from Lu et al. 2022
if False:
    print('')
    print('- Creating Figure 2, for ITP3')
    # Make the subplot groups
    group_Lu2022_clstr = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_Lu2022_clstr)
    group_Lu2022_CT_SP = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_Lu2022_CT_SP)
    # Make the figure
    ahf.make_figure([group_Lu2022_clstr,group_Lu2022_CT_SP], use_same_y_axis=False)
    # ahf.make_figure([group_Lu2022_CT_SP])

### Figure 3
## Parameter sweep across \ell and m_pts for ITP2
if False:
    print('')
    print('- Creating Figure 3')
    # Make the subplot groups
    group_ps_ell = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_ps_ell, plot_title='')
    group_ps_m_pts = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_ps_m_pts, plot_title='')
    # Make the figure
    ahf.make_figure([group_ps_ell, group_ps_m_pts], use_same_y_axis=False, filename='Figure_3.pickle')
## Parameter sweep across \ell and m_pts for ITP3
if False:
    print('')
    # Make the subplot groups
    group_ps_ell = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_ITP3_ps_ell)
    group_ps_m_pts = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_ITP3_ps_m_pts)
    # Make the figure
    ahf.make_figure([group_ps_ell, group_ps_m_pts], use_same_y_axis=False, filename='ITP3_param_sweep.pickle')

### Figure 4
## Evaluating clusterings with the overlap ratio and lateral density ratio
# For the reproduction of Timmermans et al. 2008
if False:
    print('')
    print('- Creating Figure 4')
    # Make the subplot groups
    group_salt_R_L = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_salt_R_L, plot_title='')
    group_salt_nir = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_salt_nir, plot_title='')
    # Make the figure
    ahf.make_figure([group_salt_R_L, group_salt_nir], filename='Figure_4.pickle')
# For the reproduction of Lu et al. 2022
if False:
    print('')
    print('- Creating Figure 4, for ITP3')
    group_salt_R_L = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_salt_R_L)
    group_salt_nir = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_salt_nir)
    # Make the figure
    #   Confirmed
    ahf.make_figure([group_salt_R_L, group_salt_nir])

### Figure 5
## Tracking clusters across a subset of profiles
# For ITP2
if False:
    print('')
    print('- Creating Figure 5')
    # Make the subplot groups
    group_ITP2_some_pfs_0 = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_some_pfs_0, plot_title='')
    group_ITP2_some_pfs_1 = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_ITP2_some_pfs_1, plot_title='')
    # Make the figure
    ahf.make_figure([group_ITP2_some_pfs_0, group_ITP2_some_pfs_1], use_same_x_axis=False, use_same_y_axis=False, row_col_list=[2,1, 0.3, 1.70], filename='Figure_5.pickle')

### Figure 6
## Tracking clusters across profiles, reproducing Lu et al. 2022 Figure 3
if False:
    print('')
    print('- Creating Figure 6')
    # Make the subplot groups
    group_Lu2022_fig3a = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_Lu2022_fig3a, plot_title='')
    group_Lu2022_fig3b = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_Lu2022_fig3b, plot_title='')
    group_Lu2022_fig3c = ahf.Analysis_Group(ds_ITP3_all, pfs_ITP3, pp_Lu2022_fig3c, plot_title='')
    # Make the figure
    #   Confirmed that it will pull from netcdf. Takes a long time to run still
    ahf.make_figure([group_Lu2022_fig3a, group_Lu2022_fig3b, group_Lu2022_fig3c], filename='Figure_6.pickle')
# For ITP2
if False:
    print('')
    print('- Creating Figure 6, for ITP2')
    # Make the subplot groups
    group_Lu2022_fig3a = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3a)
    group_Lu2022_fig3b = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3b, plot_title='')
    group_Lu2022_fig3c = ahf.Analysis_Group(ds_ITP2_all, pfs_ITP2, pp_Lu2022_fig3c, plot_title='')
    # Make the figure
    ahf.make_figure([group_Lu2022_fig3a, group_Lu2022_fig3b, group_Lu2022_fig3c])

### Figure 8
## Tracking clusters across a subset of profiles
# For ITP3
if False:
    print('')
    print('- Creating Figure 8')
    # Make the subplot group
    group_ITP3_some_pfs_2 = ahf.Analysis_Group(ds_ITP3_some_pfs2, pfs_ITP3, pp_ITP3_some_pfs_2)
    # Make the figure
    ahf.make_figure([group_ITP3_some_pfs_2], row_col_list=[1,1, 0.27, 0.90], filename='Figure_8.pickle')

################################################################################