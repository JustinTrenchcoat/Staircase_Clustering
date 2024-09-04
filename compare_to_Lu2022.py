"""
Author: Mikhail Schee
Created: 2023-03-13

This script will plot average temperature and salinity data from two csv files.
I specifically made this to compare my clustering results to those of Lu et al.
2022 with the functions I made in `analysis_helper_functions.py`

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions in source code must retain the accompanying copyright notice, this list of conditions, and the following disclaimer.
    2. Redistributions in binary form must reproduce the accompanying copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3. Names of the copyright holders must not be used to endorse or promote products derived from this software without prior written permission from the copyright holders.
    4. If any files are modified, you must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import matplotlib.pyplot as plt
import pandas as pd
import analysis_helper_functions as ahf
# Import the Thermodynamic Equation of Seawater 2010 (TEOS-10) from GSW
import gsw
import dill as pl

Lu2022_csv = 'outputs/Lu2022_Table_A1.csv'
my_csv = 'outputs/ITP3_cluster_table.csv'
# Try to load the specified csvs
print('- Loading '+Lu2022_csv)
try:
    Lu2022_df = pd.read_csv(Lu2022_csv)
except:
    print('Could not load '+Lu2022_csv)
    exit(0)
print('- Loading '+my_csv)
try:
    my_df = pd.read_csv(my_csv)
except:
    print('Could not load '+my_csv)
    exit(0)


################################################################################
# Declare variables for plotting
################################################################################
# Be sure to set dark_mode similiarly in the other scripts
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

# Set some plotting styles
Lu_mrk_size = 20
my_mrk_size = 30
std_marker = '.'
mpl_clrs = ahf.distinct_clrs
mpl_mrks = ahf.mpl_mrks

################################################################################
# Plotting functions ###########################################################
################################################################################

def make_figure(Lu2022_df, my_df, pp, filename=None, use_same_y_axis=None):
    """
    Takes in a list of Analysis_Group objects, one for each subplot. Determines
    the needed arrangement of subplots, then passes one Analysis_Group object to
    each axis for plotting

    Lu2022_df       A pandas DataFrame
    my_df           A pandas DataFrame
    pp              A custom Plot Parameters object
    filename        A string of the file to which to output the plot 
    """
    fig, ax = ahf.set_fig_axes([1], [1], fig_ratio=0.8, fig_size=1.25)
    xlabel, ylabel, invert_y_axis = plot_comparison(ax, Lu2022_df, my_df, pp)
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
    if invert_y_axis:
        ax.invert_yaxis()
    ax.set_title('Comparing to Lu et al. 2022')
    #
    plt.tight_layout()
    #
    if filename != None:
        print('- Saving figure to outputs/'+filename)
        if '.png' in filename:
            plt.savefig('outputs/'+filename, dpi=400, transparent=True)
        elif '.pickle' in filename:
            pl.dump(fig, open('outputs/'+filename, 'wb'))
        else:
            print('File extension not recognized in',filename)
    else:
        print('- Displaying figure')
        plt.show()

################################################################################

def get_axis_label(var_key):
    """
    Takes in a variable name and returns a nicely formatted string for an axis label

    var_key             A string of the variable name
    """
    if 'SP' in var_key:
        ax_label = r'$S_P$ (g/kg)'
    elif 'SA' in var_key:
        ax_label = r'$S_A$ (g/kg)'
    elif 'CT' in var_key:
        ax_label = r'$\Theta$ ($^\circ$C)'
    elif 'PT' in var_key:
        ax_label = r'$\theta$ ($^\circ$C)'
    elif 'depth' in var_key:
        ax_label = r'depth (m)'
    elif 'press' in var_key:
        ax_label = r'pressure (dbar)'
    # Check for certain modifications to variables,
    #   check longer strings first to avoid mismatching
    # Check for cluster average variables
    if 'ca_' in var_key:
        return 'Cluster average of '+ ax_label
    # Check for cluster span variables
    elif 'cs_' in var_key:
        return 'Cluster span of '+ ax_label
    else:
        return 'None'

################################################################################

def plot_comparison(ax, Lu2022_df, my_df, pp):
    """
    Plots the average values of clusters from Lu2022 and also from my results

    ax              The axis on which to make the plot
    df              A pandas DataFrame
    pp              A custom Plot Parameters object
    """
    # Set the main x and y data keys
    x_key = pp.x_vars[0]
    y_key = pp.y_vars[0]
    # Get the axis labels
    xlabel = get_axis_label(x_key)
    ylabel = get_axis_label(y_key)
    # print('xlabel:',xlabel,'ylabel:',ylabel)
    # Check to make sure those values are available
    #   Assume a general position of 76 N, 135 W
    temp_lat = [75]*len(Lu2022_df)
    temp_lon = [-135]*len(Lu2022_df)
    if 'ca_press' in [x_key, y_key]:
        # Depth values must be negative
        Lu2022_df['ca_press'] = gsw.p_from_z(-Lu2022_df['ca_depth'], temp_lat)
    if 'ca_CT' in [x_key, y_key]:
        Lu2022_df['ca_press'] = gsw.p_from_z(-Lu2022_df['ca_depth'], temp_lat)
        Lu2022_df['ca_SA'] = gsw.SA_from_SP(Lu2022_df['ca_SP'].values, Lu2022_df['ca_press'].values, temp_lon, temp_lat)
        Lu2022_df['ca_CT'] = gsw.CT_from_pt(Lu2022_df['ca_SA'], Lu2022_df['ca_PT'])
    # print(Lu2022_df)
    # Set variable as to whether to invert the y axis
    if 'press' in y_key or 'depth' in y_key:
        invert_y_axis = True
    else:
        invert_y_axis = False
    # Plot every point the same color, size, and marker for Lu et al. 2022
    ax.scatter(Lu2022_df[x_key], Lu2022_df[y_key], color=std_clr, s=Lu_mrk_size, marker=std_marker, zorder=5, label='Lu et al. 2022')
    # Plot my clusters, selecting color and marker shape by cluster id
    for i in my_df['cluster'].values:
        # Decide on the color and symbol, don't go off the end of the arrays
        my_clr = mpl_clrs[i%len(mpl_clrs)]
        my_mkr = mpl_mrks[i%len(mpl_mrks)]
        # Find the data from this cluster
        df_this_cluster = my_df[my_df['cluster']==i]
        ax.scatter(df_this_cluster[x_key], df_this_cluster[y_key], color=my_clr, s=my_mrk_size, marker=my_mkr, zorder=3)
        # Mark outliers
        # if df_this_cluster['out_1st'].values[0]:
        #     ax.scatter(df_this_cluster[x_key], df_this_cluster[y_key], edgecolors='r', s=my_mrk_size*5, marker='o', facecolors='none', zorder=2)
        # if df_this_cluster['out_2nd'].values[0]:
        #     ax.scatter(df_this_cluster[x_key], df_this_cluster[y_key], edgecolors='b', s=my_mrk_size*5, marker='o', facecolors='none', zorder=2)
        if df_this_cluster['out_nir_SP'].values[0] or df_this_cluster['out_cRL'].values[0]:
            ax.scatter(df_this_cluster[x_key], df_this_cluster[y_key], edgecolors='r', s=my_mrk_size*5, marker='o', facecolors='none', zorder=2)
    # Add a standard legend
    ax.legend()
    # Add grid
    ax.grid(color=std_clr, linestyle='--', alpha=0.5)
    return xlabel, ylabel, invert_y_axis

################################################################################

# pp_Lu2022 = ahf.Plot_Parameters(x_vars=['ca_SP'], y_vars=['ca_CT'], clr_map='clr_all_same')
# make_figure(Lu2022_df, my_df, pp_Lu2022)#, filename='Figure_7.pickle')