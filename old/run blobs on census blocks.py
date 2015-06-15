### census blocks code
### contributors: jgiuffrida@uchicago.edu
### 5/11/15

# runs blobs at the census block level; still very exploratory
# by default uses np.ones as the floor variable (so min_pop=5 means
# all blobs must have at least five census blocks in them)

import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
import matplotlib.pyplot as plt
import time
import cmd
import datetime
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans


%cd 'Documents/tech/git/blobs'

shp_link = 'blocks/CensusBlockTIGER2010.shp'
w=ps.open('blocks/CensusBlockTIGER2010.gal').read()
calls = pd.read_csv('final block data.csv')


# histogram helper function
def hist(data, title='Histogram of Values', bins=20, range=None):
    """Create a nice-looking histogram.

    Parameters
    ----------

    data            : array
                      n*1 vector of observations on variable of interest

    title           : string
                      title for the chart

    bins            : int
                      number of bins in which to group observations

    Attributes
    ----------

    none

    Examples
    --------

    >>> import numpy as np
    >>> x = np.random.rand(100,1)
    >>> hist(x, title='Random Uniform Distribution', bins=8)

    """

    if not range:
        range = (data.min(), data.max())

    hist, bins = np.histogram(data, bins=bins, range=range)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.title(title)
    plt.show()


# helper function to determine how to color choropleth
def sort_regions(self, method='objective'):
    sr = np.zeros([self.k,2])
    for region in range(0,self.k):
        sr[region][0] = region
        selectionIDs = [self.w.id_order.index(i) for i in self.regions[region]]
        m = self.z[selectionIDs, :]
        if method == 'objective':
            var = m.var(axis=0)
            sr[region][1] = sum(np.transpose(var)) * len(self.regions[region])
        elif method == 'mean':
            sr[region][1] = m.mean()  # simple mean of all variables
    srdf = pd.DataFrame(sr)
    srdf = srdf.sort(columns=1)
    self.sorted_regions = dict()
    for i in range(0,self.k):
        self.sorted_regions[int(srdf[i:i+1][0])] = i

ps.Maxp.sort_regions = sort_regions

# helper function to assign weights to variables
def format_blobs(data, option='equal votes', weights=None):
    if option == 'default':
        # use max p as originally designed
        # variables will be implicitly weighted in proportion to their means
        return data
    elif option == 'equal votes':
        # give equal weight to all variables by standardizing them
        x = np.zeros(data.shape)
        for v in range(data.shape[1]):
            x[:,v] = (data[:,v] - np.mean(data[:,v])) / np.std(data[:,v])
        return x
    elif option == 'weighted':
        # assign explicit weights to standardized variables
        x = np.zeros(data.shape)
        for v in range(data.shape[1]):
            x[:,v] = ((data[:,v] - np.mean(data[:,v])) / \
                np.std(data[:,v])) * np.sqrt(weights[v])
        return x

# helper function to retrieve original, non-standardized data about a blob
def retrieve_raw_data():
    pass  # todo


# main blobs method
def blobs(v, min_pop, floor_var='pop', iterations=10, method='equal votes', weights=[], 
    initial=10, plot=True, savedata=False, plot_values=False, verbose=True):
    """Create a max-p regions solution for a given shapefile and associated 
    dataset. Builds on pysal.Maxp with improvements to the user interface, 
    verbosity, and mapping. 

    Original problem from "The Max-p-Regions Problem," Duque, Anselin, and Rey, 
    JRS, October 2010, available at http://geography.sdsu.edu/Research/
    Projects/IPC/publication/MaxP_authored.pdf.

    Parameters
    ----------
    v           : array
                  array of variables on which to create blobs (for all 
                    variables, use ['all'])

    min_pop     : int
                  minimum population in each blob

    iterations  : int
                  number of blobs solutions to create (will return best): 10 by 
                    default

    method      : {'equal votes', 'default', 'weighted'}
                  equal votes' by default, can change to 'weighted'

    weights     : array
                  if method='weighted', add weights for variables as an array

    initial     : int
                  number of times to revise each solution (10 by default)

    plot        : boolean
                  will plot the best solution (True by default)

    savedata    : boolean
                  will save a CSV of the blobs data to the root folder (False 
                    by default)

    plot_values : boolean
                  will color-code the plot by the mean of the underlying 
                    variables. only makes sense with one variable. default 
                    False (plots by ID of the blob)
    
    Sample usage
    ------------

    >>> blobs(['all_calls_per1000'], min_pop=10000, plot_values=True)

    """
    
    solutions = []
    top_scores = []
    times = []
    num_blobs = []
    current_time = []
    iteration = []
    best_score = -1
    best_solution = None
    floor_var_array = np.ones((calls.shape[0],1))   #################### changed this
    blob_vars = np.array(calls.loc[:,v], np.float64)
    if len(v) == 1:
        # add shape to the array
        blob_vars.shape = (blob_vars.shape[0], 1)
    print('\n### CREATING BLOBS FROM ' + str(len(v)) + 
        ' VARIABLES ###\n    PARAMETERS:\n     # Minimum ' + floor_var + ' in each blob: ' + 
        str(int(min_pop)) + '\n     # Iterations: ' + str(iterations) +
        '\n     # Method: ' + method + '\n     # Plot blobs: ' + str(plot) + 
        '\n     # Save blobs data: ' + str(savedata) + '\n')
    for i in range(0,iterations):
        start = time.time()
        r=ps.Maxp(w, format_blobs(blob_vars, method, weights=weights), 
            floor=min_pop, floor_variable=floor_var_array, initial=initial, verbose=verbose)
        end = time.time()
        times.append(end - start)
        current_time.append(end)
        current_score = r.objective_function()
        solutions.append(current_score)
        num_blobs.append(r.k)
        if (best_score == -1 or current_score < best_score):
            best_score = current_score
            best_solution = r
        top_scores.append(best_score)
        iteration.append(i)
        msg = '\r# ITERATION '+str(i+1)+'                 \n  Score: ' + \
            str(round(current_score,2)) + '\n  Created '+str(r.k)+' blobs (' + \
            str(int(calls.shape[0]/r.k)) + ' tracts per blob)\n  Best solution so far: ' + \
            str(round(best_score,2))
        msg += '\n  Time taken: '+str(round(end-start,1))+' seconds ('+ \
            str(int(np.mean(times)*(iterations-i-1)))+' seconds remaining)\n'
        print msg
    
    r = best_solution
    print('\r# BEST SOLUTION:                      \n  Score: '+
        str(round(r.objective_function(),2)) + 
        '\n  '+str(r.k)+' blobs ('+str(int(calls.shape[0]/r.k))+
        ' blocks per blob)')
    if plot:
        print('  Plotting...'),
        # prep for plotting
        ids=np.array(calls['block_id']).astype(str)
        if plot_values:
            r.sort_regions(method='mean')  # sort regions by intensity of the variable
        regions=np.empty(calls.shape[0])
        for j in range(0,calls.shape[0]):
            reg=r.area2region[ids[j]]
            regions[j]=reg
        # show blobs we created
        maps.plot_choropleth(shp_link, regions, type='quantiles',
            title='Chicago blobs from census tracts\n(min ' + 
                str(int(r.floor)) +' population per blob, ' + 
                str(r.p)+' blobs)', k=r.p, figsize=(6,9))
        print('\r             \n')
    
    #build data structure
    sr = np.zeros([r.k, len(v)*2+4])
    for region in range(0,r.k):
        # blob ID
        sr[region][0] = region
        selectionIDs = [r.w.id_order.index(i) for i in r.regions[region]]
        m = r.z[selectionIDs, :]
        # objective function
        var = m.var(axis=0)
        sr[region][1] = sum(np.transpose(var)) * len(r.regions[region])
        # blob size (number of places in blob)
        sr[region][2] = len(r.regions[region])
        # blob population
        sr[region][3] = calls.loc[selectionIDs, floor_var].sum()
        # variable means and standard deviations
        for j in range(0,len(v)):
            sr[region][4+j*2] = m[:,j].mean()
            sr[region][5+j*2] = m[:,j].std()
    srdf = pd.DataFrame(sr)
    cols = ['Blob', 'Score', 'Size', floor_var]
    for j in range(0, len(v)):
        cols.append(v[j]+'_mean')
        cols.append(v[j]+'_stdev')
    srdf.columns = cols
    if savedata:
        srdf.to_csv('Blobs data ' + datetime.datetime.now().strftime('%Y%m%d %H%M') + \
            '.csv', index=False)
    return dict(best=r, data=srdf, regions=r.area2region)


r = blobs(['graffiti_removal', 'rodent_baiting', 'tree_trims'], 20, iterations=1, 
  initial=1, savedata=True)

