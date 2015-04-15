### blobs code
### jgiuffrida@uchicago.edu
### April 2015

###########
### PREP (run to end)
###########

import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
import matplotlib.pyplot as plt
import time
shp_link = 'tracts/CensusTractsTIGER2010.shp'
dbf = ps.open('tracts/CensusTractsTIGER2010.dbf')
cols = np.array([dbf.by_col(col) for col in dbf.header]).T
df = pd.DataFrame(cols)
df.columns = dbf.header
df.columns = df.columns.map(lambda x: x.lower())
df.commarea = df.commarea.astype('int')
df['order'] = df.index
w=ps.open('tracts/CensusTractsTIGER2010_fixed.gal').read()
calls = pd.read_csv('master311.csv', dtype=object)
for c in calls.columns[1:]:
    calls[c] = calls[c].astype('float')

ordered_tracts = pd.DataFrame(df.loc[:,['tractce10', 'commarea', 'order']])
calls = pd.merge(calls, ordered_tracts, how='right', left_on='tract', 
    right_on='tractce10', sort=False).fillna(0)
calls = calls.sort(['order'])

def hist(data, title='', bins=20):
    hist, bins = np.histogram(data, bins=20)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.title(title)
    plt.show()

def sort_regions(self, method='objective'):
    sr = np.zeros([self.k,2])
    for region in range(0,r.k):
        sr[region][0] = region
        selectionIDs = [r.w.id_order.index(i) for i in r.regions[region]]
        m = r.z[selectionIDs, :]
        if method == 'objective':
            var = m.var(axis=0)
            sr[region][1] = sum(np.transpose(var)) * len(r.regions[region])
        elif method == 'mean':
            sr[region][1] = m.mean()  # simple mean of all variables
    srdf = pd.DataFrame(sr)
    srdf = srdf.sort(columns=1)
    self.sorted_regions = dict()
    for i in range(0,self.k):
        self.sorted_regions[int(srdf[i:i+1][0])] = i

ps.Maxp.sort_regions = sort_regions
def format_blobs(data, option='default', weights=None):
    if option == 'default':
        # use max p as originally designed
        # variables will be implicitly weighted in proportion to their means
        return data
    elif option == 'equal votes':
        # give equal weight to all variables by standardizing them
        x = np.zeros(blob_vars.shape)
        for v in range(0, blob_vars.shape[1]):
            x[:,v] = (blob_vars[:,v] - np.mean(blob_vars[:,v])) / np.std(blob_vars[:,v])
        return x
    elif option == 'weighted':
        # assign explicit weights to standardized variables
        x = np.zeros(blob_vars.shape)
        for v in range(0, blob_vars.shape[1]):
            x[:,v] = ((blob_vars[:,v] - np.mean(blob_vars[:,v])) / \
                np.std(blob_vars[:,v])) * np.sqrt(weights[v])
        return x

v = ['vehicles_per1000', 'alley_lights_per1000', 'garbage_per1000',
    'graffiti_per1000', 'potholes_per1000', 'rodents_per1000', 
    'sanitation_per1000', 'street_lights_one_per1000', 'street_lights_all_per1000', 
    'tree_debris_per1000', 'tree_trims_per1000', 'buildings_per1000']
blob_vars = np.array(calls.loc[:,v], np.float64)

def blobs(vars, min_pop, iterations, method='equal votes', weights=[], 
    initial=10, plot=True):
    solutions = []
    top_scores = []
    times = []
    num_blobs = []
    current_time = []
    iteration = []
    best_score = 10**12
    best_solution = None
    for i in range(0,iterations):
        start = time.time()
        r=ps.Maxp(w, format_blobs(blob_vars, method, weights=weights), 
            floor=min_pop, floor_variable=calls['pop'], initial=initial)
        end = time.time()
        times.append(end - start)
        current_time.append(end)
        solutions.append(r.objective_function())
        num_blobs.append(r.k)
        if (r.objective_function() < best_score):
            best_score = r.objective_function()
            best_solution = r
        top_scores.append(best_score)
        iteration.append(i)
        print('iteration '+str(i+1)+' - score: '+str(round(r.objective_function(),2))+
        ' ('+str(r.k)+' blobs), best: '+str(round(best_score,2))+', '+
        str(round(end-start,1))+'s')
    if plot:
        # prep for plotting
        r = best_solution
        ids=np.array(calls['tractce10']).astype(str)
        # r.sort_regions(method='mean')  # uncomment to sort regions by intensity of the variable
        regions=np.empty(calls.shape[0])
        for j in range(0,calls.shape[0]):
            reg=r.area2region[ids[j]]
            regions[j]=reg
        # show blobs we created
        maps.plot_choropleth(shp_link, regions, type='quantiles',
            title='Chicago blobs from census tracts\n(min ' + 
                str(r.floor) +' population per blob, ' + 
                str(r.p)+' blobs)', k=r.p, figsize=(6,8))
    return dict(times=times, solutions=solutions, top_scores=top_scores, 
        current_time=current_time, iteration=iteration, best=r)


#############
#### END PREP
#############

# blobs() is our function that creates blobs
# usage:
#   vars: variables on which to create blobs; all variables by default. to 
#         use only certain variables, set this to 
#         np.array(calls.loc[:,[x, y, z]], np.float64),
#         where x, y, z are the string names of variables desired, e.g.
#         np.array(calls.loc[:,['tree_trims_per1000', 'buildings_per1000']], np.float64)
#   min_pop: minimum population in each blob
#   iterations: number of blobs solutions to create (will return best)
#   method: 'equal votes' by default, can change to 'weighted'
#   weights: if method='weighted', add weights for variables as an array
#   initial: number of times to revise each solution (10 by default)
#   plot: will plot the best solution (True by default)

# the script will print the scores for each solution and the time it took.
# time increases quickly with initial and min_pop.
# interesting things to try: change min_pop to numbers between 10000 and 500000,
# use only 1 or 2 vars, or change method to 'weighted' and add weights such as
# [100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] (remember, as many weights as variables)


results_10k = blobs(blob_vars, 10000, 1)
results_100k = blobs(blob_vars, 100000, 1)
results_trees = blobs(np.array(calls.loc[:,['tree_trims_per1000', 
    'tree_debris_per1000']], np.float64), 100000, 1)
results_dirtiness = blobs(np.array(calls.loc[:,['garbage_per1000', 
    'rodents_per1000', 'sanitation_per1000']], np.float64), 100000, 1)



















