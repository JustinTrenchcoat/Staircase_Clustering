"""
Author: Mikhail Schee
Created: 2023-02-27

This script will take in the name of a csv file and plot the data within it
I specifically made this to do parameter sweeps that take too long to plot
with the functions I made in `analysis_helper_functions.py`

Usage:
    plot_from_csv.py CSV

Options:
    CSV             # filepath of the csv to plot
"""
import matplotlib.pyplot as plt
# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
my_csv = args['CSV']       # filename of the pickle to unpickle

import pandas as pd
# Try to load the specified csv
print('- Loading '+my_csv)
try:
    df = pd.read_csv(my_csv)
except:
    print('Could not load '+my_csv)
    exit(0)

# Display the figure in the interactive matplotlib GUI
print('- Displaying figure')
print(df)

exit(0)

pp_ITP3_ps_m_pts = ahf.Plot_Parameters(x_vars=['m_pts'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'z_var':'maw_size', 'z_list':[50]})
ahf.make_figure(df, pp_ITP3_ps_m_pts, filename='ITP3_test_sweep.pickle')

def make_figure(df, pp, filename=None, use_same_y_axis=None):
    """
    Takes in a list of Analysis_Group objects, one for each subplot. Determines
    the needed arrangement of subplots, then passes one Analysis_Group object to
    each axis for plotting

    df              A pandas DataFrame
    pp              A custom Plot Parameters object
    filename        A string of the file to which to output the plot 
    """
    fig, ax = set_fig_axes([1], [1], fig_ratio=0.8, fig_size=1.25)
    xlabel, ylabel, plt_title, ax = make_subplot(ax, df, pp, fig, 111)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # If axes limits given, limit axes
    if not isinstance(pp.ax_lims, type(None)):
        try:
            ax.set_xlim(pp.ax_lims['x_lims'])
            print('\tSet x_lims to',pp.ax_lims['x_lims'])
        except:
            foo = 2
        try:
            ax.set_ylim(pp.ax_lims['y_lims'])
            print('\tSet y_lims to',pp.ax_lims['y_lims'])
        except:
            foo = 2
    ax.set_title(plt_title)
    #
    plt.tight_layout()
    #
    if filename != None:
        print('- Saving figure to outputs/'+filename)
        if '.png' in filename:
            plt.savefig('outputs/'+filename, dpi=400)
        elif '.pickle' in filename:
            pl.dump(fig, open('outputs/'+filename, 'wb'))
        else:
            print('File extension not recognized in',filename)
    else:
        print('- Displaying figure')
        plt.show()

def make_subplot(ax, df, pp, fig, ax_pos):
    """
    Takes in an Analysis_Group object which has the data and plotting parameters
    to produce a subplot. Returns the x and y labels and the subplot title

    ax              The axis on which to make the plot
    df              A pandas DataFrame
    pp              A custom Plot Parameters object
    fig             The figure in which ax is contained
    ax_pos          A tuple of the ax (rows, cols, linear number of this subplot)
    """
    ## Make a standard x vs. y scatter plot
    # Set the main x and y data keys
    x_key = pp.x_vars[0]
    y_key = pp.y_vars[0]
    # Check for twin data key
    try:
        tw_y_key = pp.y_vars[1]
        tw_ax_x  = ax.twinx()
    except:
        tw_y_key = None
        tw_ax_x  = None
    # Add a standard title
    plt_title = 'ITP3'
    # Plot the parameter sweep
    xlabel, ylabel = plot_clstr_param_sweep(ax, tw_ax_x, df, pp, plt_title)

def plot_clstr_param_sweep(ax, tw_ax_x, df, pp, plt_title=None):
    """
    Plots the number of clusters found by HDBSCAN vs. the number of profiles
    included in the data set

    ax              The axis on which to make the plot
    df              A pandas DataFrame
    pp              A custom Plot Parameters object
    """
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
    # Check for a z variable
    try:
        z_key = cluster_plt_dict['z_var']
        z_list = cluster_plt_dict['z_list']
    except:
        z_key = None
        z_list = [0]
    # Get a list of the column headings from the dataframe
    # Check to make sure the x, y, and z keys are in the dataframe
    # Set x and y labels
    if x_key == 'maw_size':
        xlabel = r'Moving average window $\ell_{maw}$ (dbar)'
    elif x_key == 'm_pts':
        xlabel = r'Minimum density threshold $m_{pts}$'
    if y_key == 'DBCV':
        ylabel = 'DBCV'
    elif y_key == 'n_clusters':
        ylabel = 'Number of clusters'
    if tw_y_key:
        if tw_y_key == 'DBCV':
            tw_ylabel = 'DBCV'
        elif tw_y_key == 'n_clusters':
            tw_ylabel = 'Number of clusters'
        #
    #
    for i in range(len(z_list)):
        if z_key == 'maw_size':
            # Need to apply moving average window to original data, before
            #   the data filters were applied, so make a new Analysis_Group
            a_group.profile_filters.m_avg_win = z_list[i]
            new_a_group = Analysis_Group(a_group.data_set, a_group.profile_filters, a_group.plt_params)
            this_df = pd.concat(new_a_group.data_frames)
            zlabel = r'$\ell_{maw}=$'+str(z_list[i])+' dbar'
        
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
        new_df, rel_val = HDBSCAN_(this_df, cl_x_var, cl_y_var, m_pts, min_samp=min_s)
        # Record outputs to plot
        ax.plot(x_var_array, y_var_array, color=std_clr, linestyle=l_styles[i], label=zlabel)
        if tw_y_key:
            tw_ax_x.plot(x_var_array, tw_y_var_array, color=alt_std_clr, linestyle=l_styles[i])
            tw_ax_x.set_ylabel(tw_ylabel)
            # Change color of the axis label on the twin axis
            tw_ax_x.yaxis.label.set_color(alt_std_clr)
            # Change color of the ticks on the twin axis
            tw_ax_x.tick_params(axis='y', colors=alt_std_clr)
    if z_key:
        ax.legend()
    return xlabel, ylabel
