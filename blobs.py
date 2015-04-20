### blobs code, v0.2
### contributors: jgiuffrida@uchicago.edu
### 4/20/15

###########
### PREP (run to end)
###########

# import modules and data
import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
import matplotlib.pyplot as plt
import time
import cmd
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

# format data and merge on shapefile IDs
ordered_tracts = pd.DataFrame(df.loc[:,['tractce10', 'commarea', 'order']])
calls = pd.merge(calls, ordered_tracts, how='right', left_on='tract', 
    right_on='tractce10', sort=False).fillna(0)
calls = calls.sort(['order'])

all_vars = ['vehicles_per1000', 'alley_lights_per1000', 'garbage_per1000',
    'graffiti_per1000', 'potholes_per1000', 'rodents_per1000', 
    'sanitation_per1000', 'street_lights_one_per1000', 'street_lights_all_per1000', 
    'tree_debris_per1000', 'tree_trims_per1000', 'buildings_per1000']


# histogram helper function
def hist(data, title='', bins=20):
    hist, bins = np.histogram(data, bins=20)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.title(title)
    plt.show()

# helper function to determine how to color choropleth
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

# helper function to assign weights to variables
def format_blobs(data, option='default', weights=None):
    if option == 'default':
        # use max p as originally designed
        # variables will be implicitly weighted in proportion to their means
        return data
    elif option == 'equal votes':
        # give equal weight to all variables by standardizing them
        x = np.zeros(data.shape)
        for v in range(0, data.shape[1]):
            x[:,v] = (data[:,v] - np.mean(data[:,v])) / np.std(data[:,v])
        return x
    elif option == 'weighted':
        # assign explicit weights to standardized variables
        x = np.zeros(data.shape)
        for v in range(0, data.shape[1]):
            x[:,v] = ((data[:,v] - np.mean(data[:,v])) / \
                np.std(data[:,v])) * np.sqrt(weights[v])
        return x

# main blobs method

def blobs(v, min_pop, floor_var='pop', iterations=10, method='equal votes', weights=[], 
    initial=10, plot=True):
    # usage:
    #   v: array of variables on which to create blobs (for all variables, use ['all'])
    #   min_pop: minimum population in each blob
    #   iterations: number of blobs solutions to create (will return best): 10 by default
    #   method: 'equal votes' by default, can change to 'weighted'
    #   weights: if method='weighted', add weights for variables as an array
    #   initial: number of times to revise each solution (10 by default)
    #   plot: will plot the best solution (True by default)
    solutions = []
    top_scores = []
    times = []
    num_blobs = []
    current_time = []
    iteration = []
    best_score = -1
    best_solution = None
    floor_var_array = calls[floor_var]
    if v == ['all']:
        v = all_vars
    blob_vars = np.array(calls.loc[:,v], np.float64)
    print('\n### CREATING BLOBS FROM ' + str(len(v)) + 
        ' VARIABLES ###\n    PARAMETERS:\n     # Minimum ' + floor_var + ' in each blob: ' + 
        str(int(min_pop)) + '\n     # Iterations: ' + str(iterations) +
        '\n     # Method: ' + method + '\n')
    for i in range(0,iterations):
        start = time.time()
        r=ps.Maxp(w, format_blobs(blob_vars, method, weights=weights), 
            floor=min_pop, floor_variable=floor_var_array, initial=initial)
        end = time.time()
        times.append(end - start)
        current_time.append(end)
        solutions.append(r.objective_function())
        num_blobs.append(r.k)
        if (best_score == -1 or r.objective_function() < best_score):
            best_score = r.objective_function()
            best_solution = r
        top_scores.append(best_score)
        iteration.append(i)
        print('\r# ITERATION '+str(i+1)+'                 \n  Score: '+
            str(round(r.objective_function(),2))+
            '\n  Created '+str(r.k)+' blobs ('+str(int(calls.shape[0]/r.k))+
            ' tracts per blob)\n  Best solution so far: '+str(round(best_score,2))+
            '\n  Time taken: '+str(round(end-start,1))+' seconds ('+
            str(int(np.mean(times)*(iterations-i-1)))+' seconds remaining)\n')
    r = best_solution
    print('\r# BEST SOLUTION:                      \n  Score: '+
        str(round(r.objective_function(),2))+
        '\n  '+str(r.k)+' blobs ('+str(int(calls.shape[0]/r.k))+
        ' tracts per blob)')
    if plot:
        print('  Plotting...'),
        # prep for plotting
        ids=np.array(calls['tractce10']).astype(str)
        # r.sort_regions(method='mean')  # uncomment to sort regions by intensity of the variable
        regions=np.empty(calls.shape[0])
        for j in range(0,calls.shape[0]):
            reg=r.area2region[ids[j]]
            regions[j]=reg
        # show blobs we created
        maps.plot_choropleth(shp_link, regions, type='quantiles',
            title='Chicago blobs from census tracts\n(min ' + 
                str(int(r.floor)) +' population per blob, ' + 
                str(r.p)+' blobs)', k=r.p, figsize=(6,9))
        print('\r           \n')
    return dict(times=times, solutions=solutions, top_scores=top_scores, 
        current_time=current_time, iteration=iteration, best=r)

# command line processor
class CmdBlobs(cmd.Cmd):
    """Command line interface for blobs"""
    
    variables = list(calls.columns)
    picked = []
    floor_var = ''
    floor_size = 0

    step = 1
    
    def do_select(self, name):
        "Step 1: select a variable"
        if name and name in self.variables:
            response = '\n  Added %s.\n  Add more variables or enter command \'next\'.\n' % name
            self.picked.append(name)
        elif name:
            response = '\nError: could not find variable %s\n' % name
        else:
            response = '\nError: please give a variable\n'
        print response
    
    def complete_select(self, text, line, begidx, endidx):
        "Autocomplete variable selection"
        if not text:
            completions = [ f
                            for f in self.variables
                            if f not in self.picked
                            ]
        else:
            completions = [ f
                            for f in self.variables
                            if (f.startswith(text)
                            and f not in self.picked)
                            ]
        return completions
    
    def do_floor(self, name):
        "Step 2 part 1: Select floor variable"
        if name and name in self.variables:
            response = '\n  Set %s as floor variable.\n' % name
            self.floor_var = name
            response += '\n  Now, please enter command '+\
                '\'size\' followed by the minimum floor size\n  you would like to set. '+\
                'To help, we\'ve provided you with a histogram \n  of %s values.\n' % name
            print response
            hist(calls[name], title='Histogram of '+name)
        elif name:
            response = '\nError: could not find variable %s\n' % name
            print response
        else:
            response = '\nError: please give a variable\n'
            print response

    def complete_floor(self, text, line, begidx, endidx):
        "Autocomplete floor variable selection"
        if not text:
            completions = self.variables[:]
        else:
            completions = [ f
                            for f in self.variables
                            if f.startswith(text)
                            ]
        return completions

    def do_size(self, size):
        "Step 2 part 2: Set size of floor"
        if float(size) > 0:
            response = '\n  Set %s as minimum floor\n' % str(size)
            self.floor_size = float(size)
            print response
            self.do_next('')
        elif size:
            response = '\nError: size should not be 0\n'
            print response
        else:
            response = '\nPlease set the size\n'
            print response

    def do_next(self, line):
        self.step += 1
        if self.step == 2:
                self.variables = list(calls.columns)
                print '\n## Step 2: Set Blob Size'
                print('  Enter command \'floor\' followed by the variable you want '+
                    'to use\n  as the \'floor\' variable. Use tab key to autocomplete.\n')
        if self.step == 3:
            response = '\n## Step 3: Run Blobs\n  Ready to run blobs. Parameters:'+\
                '\n  Variables: '
            for v in self.picked:
                response += '\n    %s' % v
            response += '\n  Floor Variable: %s' % self.floor_var
            response += '\n  Floor Size: %s' % str(self.floor_size)
            response += '\n\n  To run blobs using these parameters, enter command '+\
                '\'run\'.\n  Otherwise, enter command \'exit\'.\n'
            print response
    
    def do_exit(self, line):
        return True

    def do_run(self, line):
        blobs(self.picked, self.floor_size, self.floor_var)
        self.do_exit('')
   


def interface():
    print '\nThis is a command line interface for Blobs.'
    print '\n## Step 1: Select Variables'
    print('  Enter command \'select\' followed by a variable you want. Use tab key '+
        '\n  to see options or autocomplete. When finished, enter command \'next\'.\n'+
        '  At any time, you can enter command \'exit\' to exit.\n')
    CmdBlobs().cmdloop()


#############
#### END PREP
#############

# blobs() is our function that creates blobs.
# can also run interface() for a command-line user interface to guide 
# you through blob creation.

# the script will print the scores for each solution and the time it took.
# time increases quickly with initial and min_pop.
# interesting things to try: change min_pop to numbers between 10000 and 500000,
# use only 1 or 2 vars, or change method to 'weighted' and add weights such as
# [100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] (remember to use as many weights as variables)


# example usage:
# interface()
# results_10k = blobs(['all'], 10000)
# results_50k = blobs(['all'], 50000)
# results_trees = blobs(['tree_trims_per1000', 'tree_debris_per1000'], 10000)















