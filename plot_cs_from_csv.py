"""
Author: Mikhail Schee
Created: 2023-02-27

This script will take in the name of a csv file and plot the data within it
I specifically made this to do parameter sweeps that take too long to plot
with the functions I made in `analysis_helper_functions.py`

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions in source code must retain the accompanying copyright notice, this list of conditions, and the following disclaimer.
    2. Redistributions in binary form must reproduce the accompanying copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3. Names of the copyright holders must not be used to endorse or promote products derived from this software without prior written permission from the copyright holders.
    4. If any files are modified, you must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Usage:
    plot_cs_from_csv.py CSV

Options:
    CSV             # filepath of the csv to plot
"""
import matplotlib.pyplot as plt
import pandas as pd
import analysis_helper_functions as ahf
import dill as pl
# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
my_csv = args['CSV']       # filename of the pickle to unpickle

# Try to load the specified csv
print('- Loading '+my_csv)
try:
    df = pd.read_csv(my_csv)
except:
    print('Could not load '+my_csv)
    exit(0)

################################################################################
# Declare variables for plotting
################################################################################
dark_mode = False

# Enable dark mode plotting
if dark_mode:
    plt.style.use('dark_background')
    std_clr = 'w'
    alt_std_clr = 'yellow'
    clr_ocean = 'k'
    clr_land  = 'grey'
    clr_lines = 'w'
else:
    std_clr = 'k'
    alt_std_clr = 'olive'
    clr_ocean = 'w'
    clr_land  = 'grey'
    clr_lines = 'k'

# Define array of linestyles to cycle through
l_styles = ['-', '--', '-.', ':']

################################################################################
# Plotting functions ###########################################################
################################################################################

def make_figure(df, pp, filename=None, use_same_y_axis=None):
    """
    Takes in a list of Analysis_Group objects, one for each subplot. Determines
    the needed arrangement of subplots, then passes one Analysis_Group object to
    each axis for plotting

    df              A pandas DataFrame
    pp              A custom Plot Parameters object
    filename        A string of the file to which to output the plot 
    """
    fig, ax = ahf.set_fig_axes([1], [1], fig_ratio=0.8, fig_size=1.25)
    xlabel, ylabel, plt_title, ax = make_subplot(ax, df, pp, fig)
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

################################################################################

def make_subplot(ax, df, pp, fig):
    """
    Takes in an Analysis_Group object which has the data and plotting parameters
    to produce a subplot. Returns the x and y labels and the subplot title

    ax              The axis on which to make the plot
    df              A pandas DataFrame
    pp              A custom Plot Parameters object
    fig             The figure in which ax is contained
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
    return xlabel, ylabel, plt_title, ax

################################################################################

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
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # Loop through the different z values in the list
    for i in range(len(z_list)):
        # Set z label
        if z_key == 'maw_size':
            zlabel = r'$\ell_{maw}=$'+str(z_list[i])+' dbar'
        elif z_key == 'm_pts':
            zlabel = r'$m_{pts}=$: '+str(m_pts)
        # Plot main axis
        ax.plot(df[x_key], df[y_key], color=std_clr, linestyle=l_styles[i], label=zlabel)
        if tw_y_key:
            # Plot twin axis
            tw_ax_x.plot(df[x_key], df[tw_y_key], color=alt_std_clr, linestyle='--')#l_styles[i])
            tw_ax_x.set_ylabel(tw_ylabel)
            # Change color of the axis label on the twin axis
            tw_ax_x.yaxis.label.set_color(alt_std_clr)
            # Change color of the ticks on the twin axis
            tw_ax_x.tick_params(axis='y', colors=alt_std_clr)
    if z_key:
        if False:
            ax.legend()
    return xlabel, ylabel

################################################################################

pp_ITP3_ps_m_pts = ahf.Plot_Parameters(x_vars=['m_pts'], y_vars=['n_clusters','DBCV'], clr_map='clr_all_same', extra_args={'z_var':'maw_size', 'z_list':[100]})
make_figure(df, pp_ITP3_ps_m_pts, filename='ITP3_test_sweep.pickle')