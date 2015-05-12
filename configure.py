### blobs configuration code, v0.1
### contributors: jgiuffrida@uchicago.edu
### 5/10/15

# code for configuring a new blobs solution
# including matching up the shapefile with the raw data

import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
import matplotlib.pyplot as plt

%cd "/Users/jcgiuffrida/Documents/tech/git/blobs"

# get shapefile
shp_link = 'blocks/CensusBlockTIGER2010.shp'
maps.plot_poly_lines(shp_link)  # test shapefile

# get associated data
dbf = ps.open('blocks/CensusBlockTIGER2010.dbf')
cols = np.array([dbf.by_col(col) for col in dbf.header]).T
df = pd.DataFrame(cols)
df.columns = dbf.header
df.columns = df.columns.map(lambda x: x.lower())

# if duplicates, need to remove
len(df.ix[df.duplicated('geoid10'),:])  # number of duplicate pairs

## create weights (only need to run once)
w = ps.rook_from_shapefile(shp_link)
w.n == df.shape[0] # should be true
gal = ps.open('blocks/CensusBlockTIGER2010.gal','w')
gal.write(w)
gal.close()

df.tractce10 = df.tractce10.astype('int')
df['order'] = df.index

# plot community areas
maps.plot_choropleth(shp_link, np.array(df.tractce10), type='equal_interval',
     title='Initial Map', k=80)

# get spatial weights
w=ps.open('blocks/CensusBlockTIGER2010.gal').read()
# need to fix the ohare island (tracts 980000 and 770602)
# the following was saved to X_fixed.gal
# w.neighbors['770602'] = ['980000', '090100']
# w.weights['770602'] = [1.0, 1.0]
# w.neighbors['980000'] = ['770602', '760802']
# w.weights['980000'] = [1.0, 1.0]
# w.neighbors['090100'] = ['770602', '090200']
# w.weights['090100'] = [1.0, 1.0]
# w.neighbors['760802'] = ['980000', '760801', '170500', '170600', '760803', '770902']
# w.weights['760802'] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


# get data
calls = pd.read_csv('plenario data by block.csv', dtype=object)
# convert all columns to float
# first, replace invalid pop values with 0
for i in range(calls.shape[0]):
  try:
    calls.Total[i] = float(calls.Total[i])
  except:
    calls.Total[i] = 0
for c in calls.columns[3:]:
    calls[c] = calls[c].astype('float')

# reorder calls to have same ordering as shapefile
ordered_blocks = pd.DataFrame(df.loc[:,['geoid10', 'order']], dtype=float)
calls = pd.merge(calls, ordered_blocks, how='right', left_on='block_id', 
    right_on='geoid10', sort=False).fillna(0)
calls = calls.sort(['order'])

# double-check
tracts = calls['block_id'].map(lambda x: str(x)[6:12])
maps.plot_choropleth(shp_link, np.array(df.tractce10), type='equal_interval',
     title='Initial Map', k=80)
maps.plot_choropleth(shp_link, np.array(tracts).astype(int), type='quantiles',
     title='Initial Map', k=80)

# population-standardized (per 1000)
for c in calls.columns[5:17]:
  calls[c[21:]] = np.where(calls['Total']>0, 
    calls[c] / calls['Total']*1000, 0)


calls.to_csv('final block data.csv', index=False)

calls = pd.read_csv('final block data.csv')

# plot one example
maps.plot_choropleth(shp_link, np.array(calls.graffiti), type='quantiles',
     title='Graffiti, Per 1000', k=80)


# plot community areas from data portal
maps.plot_choropleth(shp_link, np.array(calls['311_service_requests_graffiti_removal']), 
  type='fisher_jenks', title='Blocks', k=80)

# # plot community areas from 311 data
# maps.plot_choropleth(shp_link, np.array(calls.ca), type='equal_interval',
#      title='Community Areas', k=80)

# plot population
maps.plot_choropleth(shp_link, np.array(calls['Total']), type='fisher_jenks',
     title='Census Block Population', k=80)

# # plot 311 potholes calls: absolute numbers
# maps.plot_choropleth(shp_link, np.array(calls.potholes), type='fisher_jenks',
#      title='Community Areas', k=5)

# # plot 311 potholes calls: absolute numbers
maps.plot_choropleth(shp_link, np.array(calls.rodents_per1000), type='classless',
     title='Rodent Calls by Census Area, 2011-2015, Population-Adjusted', k=6, figsize=(6,9))

# simpler histogram function
def hist(data, title='', bins=20):
    hist, bins = np.histogram(data, bins=20)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.title(title)
    plt.show()

# extend max-P with this method which reorders regions in ascending order of 
# the objective function (to show results better)
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



# this lets us easily change how blobs calculates the objective function
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



#######################
######## STUFF TO CHANGE
#######################
# vars on which to create blobs (minimum of 1
# type sdh.describe().T to see all columns
v = ['vehicles_per1000', 'alley_lights_per1000', 'garbage_per1000',
    'graffiti_per1000', 'potholes_per1000', 'rodents_per1000', 
    'sanitation_per1000', 'street_lights_one_per1000', 'street_lights_all_per1000', 
    'tree_debris_per1000', 'tree_trims_per1000', 'buildings_per1000']

method = 'equal votes'    # default, equal votes, weighted
weights = [2,1]   # only for weighted method; 0 nulls the variable
min_pop = 200000   # minimum pop in each blob (50k-500k is best)
iterations = 3   # num iterations to find best solution
# num_seeds = 6  # num initial seed census tracts


##############
#### Run the following code to the end
##############


# seeds = np.random.choice(calls['tractce10'], replace=False, size=num_seeds)
# iterate to find the best solution
blob_vars = np.array(calls.loc[:,v], np.float64)
best_score = 10**9
best_solution = None
for i in range(0,iterations):
    r=ps.Maxp(w, format_blobs(blob_vars, method, weights=weights), 
        floor=min_pop, floor_variable=calls['pop'], initial=10, myverbose=True)
    if (r.objective_function() < best_score):
        best_score = r.objective_function()
        best_solution = r
    print('iteration '+str(i+1)+' - score: '+str(round(r.objective_function(),2))+
    ' ('+str(r.k)+' blobs), best: '+str(round(best_score,2)))

# prep for plotting
r = best_solution
ids=np.array(calls['tractce10']).astype(str)
r.sort_regions(method='mean')  # uncomment to sort regions by intensity of the variable
regions=np.empty(calls.shape[0])
for j in range(0,calls.shape[0]):
    reg=r.area2region[ids[j]]
    regions[j]=reg
# show blobs we created
maps.plot_choropleth(shp_link, regions, type='quantiles',
    k=r.p, figsize=(6,8), title='Chicago Blobs\n(Based on Unemployment)')
    #title='Chicago blobs from census tracts\n(min ' + 
    #    str(r.floor) +' population per blob, ' + 
    #    str(r.p)+' blobs)', 


# check against actual values
maps.plot_choropleth(shp_link, np.array(calls.ca), type='equal_interval',
     title='Community Areas', k=80)
    
# print some stuff
print('\n\nsolution: ' + str(r.p) + ' blobs, ' + method + ' method. total score: ' + \
    str(round(r.objective_function(),1))),
# for i in range(0, len(r.regions)):
#     print('\n\nBLOB ' + str(i) + '\ntracts:'),
#     for j in r.regions[i]: 
#         print(j),
#     selectionIDs = [r.w.id_order.index(k) for k in r.regions[i]]
#     m = r.z[selectionIDs, :]
#     total_pop = sum(calls.loc[selectionIDs, 'pop'])
#     var = m.var(axis=0)
#     region_score = sum(np.transpose(var)) * len(r.regions[i])
#     print('\nblob score: ' + str(round(region_score,1)) + \
#         ' (' + str(round(region_score / r.objective_function() * 100,1)) + '% of total)')
#     print('blob population: ' + str(sum(calls.loc[selectionIDs, 'pop'])))
#     for m in v:
#         print('  ' + m + ': avg ' + str(round(np.mean(calls.loc[selectionIDs, m]),1)) + \
#             ', sd ' + str(round(np.std(calls.loc[selectionIDs, m]),1)) + \
#             ' (' + str(round(np.transpose(var)[v.index(m)]*len(r.regions[i]),1)) + \
#             ' score)')



#############
#### END
#############


###########
##### AUTOMATION
###########



# run the following
import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
import matplotlib.pyplot as plt
import time
%cd "/Users/jcgiuffrida/Documents/chicago docs/Brett/Blobs/311"
shp_link = './tracts/CensusTractsTIGER2010.shp'
dbf = ps.open('./tracts/CensusTractsTIGER2010.dbf')
cols = np.array([dbf.by_col(col) for col in dbf.header]).T
df = pd.DataFrame(cols)
df.columns = dbf.header
df.columns = df.columns.map(lambda x: x.lower())
df.commarea = df.commarea.astype('int')
df['order'] = df.index
w=ps.open('./tracts/CensusTractsTIGER2010_fixed.gal').read()
calls = pd.read_csv('./master311.csv', dtype=object)
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
    initial=10, plot=False):
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

results_100k = blobs(blob_vars, 100000, 100)

plt.plot(results_100k['times'], results_100k['solutions'], 'bs')
plt.title('Solution score vs. time taken (s)')
plt.show()












