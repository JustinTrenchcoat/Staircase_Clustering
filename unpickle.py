"""
Author: Mikhail Schee
Created: 2023-01-18

This script will take in the name of a pickle file in the `figures` folder and
unpickle it to display the contained figure in the interactive matplotlib GUI

Usage:
    unpickle.py PICKLE

Options:
    PICKLE          # filename of the pickle to unpickle in the `figures` folder
"""
import matplotlib.pyplot as plt
# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
my_pickle   = args['PICKLE']       # filename of the pickle to unpickle

import pickle as pl
# Try to unpickle the specified figure
try:
    print('Loading figures/'+my_pickle+'.pickle')
    fig = pl.load(open('figures/'+my_pickle+'.pickle', 'rb'))
except:
    print('Could not find figures/'+my_pickle+'.pickle')
    exit(0)

# Display the figure in the interactive matplotlib GUI
plt.show()

exit(0)

# Get the data from the unpickled figure
data = fig.axes[0].lines[0].get_data()
# data = fig.axes[0].images[0].get_data()
print(data)
