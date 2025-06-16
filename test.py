import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# For custom analysis functions
import plotHelper as ahf

# current working codes:

dfs0 = ahf.Data_Filters()

# Some specific profiles
ITP2_pfs  = {'ITP_2':[183, 185, 187]}
ds_ITP2_pfs = ahf.Data_Set(ITP2_pfs, dfs0)
# tests for full analysis:
# change SP_range to CT range for better scope 
test_ITP2_CT_range = [0,1]
test_pfs_ITP2 = ahf.Profile_Filters(CT_range=test_ITP2_CT_range)
test_T2008_fig4_y_lims = {'y_lims':[400,220]}
test_pp_T2008_fig4  = ahf.Plot_Parameters(x_vars=['CT'], y_vars=['press'], plot_type='profiles', clr_map='cluster', 
                                          extra_args={'b_a_w_plt':False, 'cl_x_var':'CT', 'cl_y_var':'depth', 'm_pts':170},ax_lims=test_T2008_fig4_y_lims, legend=False)
#   
test_group_T2008_fig4  = ahf.Analysis_Group(ds_ITP2_pfs, test_pfs_ITP2, test_pp_T2008_fig4,  plot_title='')

ahf.make_figure([test_group_T2008_fig4], filename='test.png')

# functions called:
# 1,2,3,8,4,7,6,15,16,24,28,19,29,30,23,17
# 1,2,3,4,6,7,8,15,16,17,19,23,24,28,29,30