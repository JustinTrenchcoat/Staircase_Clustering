"""
Author: Mikhail Schee
Created: 2023-01-18

This script will take in the name of a pickle file in the `figures` folder and
unpickle it to display the contained figure in the interactive matplotlib GUI

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions in source code must retain the accompanying copyright notice, this list of conditions, and the following disclaimer.
    2. Redistributions in binary form must reproduce the accompanying copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
    3. Names of the copyright holders must not be used to endorse or promote products derived from this software without prior written permission from the copyright holders.
    4. If any files are modified, you must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Usage:
    unpickle.py PICKLE

Options:
    PICKLE          # filepath of the pickle to unpickle
"""
import matplotlib.pyplot as plt
# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
my_pickle   = args['PICKLE']       # filename of the pickle to unpickle

import dill as pl
# Try to unpickle the specified figure
print('- Loading '+my_pickle)
try:
    fig = pl.load(open(my_pickle, 'rb'))
except:
    print('Could not load '+my_pickle)
    exit(0)

# Display the figure in the interactive matplotlib GUI
print('- Displaying figure')
plt.show()

exit(0)

# Get the data from the unpickled figure
data = fig.axes[0].lines[0].get_data()
# data = fig.axes[0].images[0].get_data()
print(data)
