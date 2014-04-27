# PySAL max-p-regions demo (blobs)
import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps  # might require basemap

# extend max-P with this method which reorders regions in ascending order of 
# the objective function (to show results better)
def sort_regions(self):
    sr = np.zeros([self.k,2])
    for region in range(0,self.k):
        sr[region][0] = region
        selectionIDs = [self.w.id_order.index(i) for i in self.regions[region]]
        sr[region][1] = self.z[selectionIDs,:].mean() # this is really really simplistic
    srdf = pd.DataFrame(sr)
    srdf = srdf.sort(columns=1)
    self.sorted_regions = dict()
    for i in range(0,self.k):
        self.sorted_regions[int(srdf[i:i+1][0])] = i
        
ps.Maxp.sort_regions = sort_regions


############################
## blobs proof of concept ##
############################

# our example (which comes with PySAL) is the rate of sudden infant death
# syndrome by counties in north carolina, 1974-1979
# the rate seems to vary spatially, especially by urbanicity

shp_link = ps.examples.get_path('sids2.shp')  # shapefile
w=ps.open(ps.examples.get_path('sids2.gal')).read()  # spatial weights
dbf = ps.open(ps.examples.get_path('sids2.dbf'))  # attributes
cols = np.array([dbf.by_col(col) for col in dbf.header]).T
numeric_cols = np.array([dbf.by_col(dbf.header[i]) for i in range(len(dbf.header)) 
    if dbf.field_spec[i][0] == 'N']).T
sidr = numeric_cols[:,12:14]  # SIDS rate in 1974 and 1979
births74 = numeric_cols[:,6]
births7479 = numeric_cols[:,[6,9]]

# generate map of birth rates in NC (proxy for population density)
maps.plot_choropleth(shp_link, births74, type='quantiles',
    title='Births', k=20)

# run blobs to create a maximal number of homogeneous regions that each have 
# at least 12,000 births
r=ps.Maxp(w, births7479, floor=12000, floor_variable=births74, initial=200)

# prep for plotting
fips=np.array(dbf.by_col('FIPS'))
r.sort_regions()
regions=np.empty(100)
for j in range(0,100):
    reg=r.area2region[fips[j]]
    regions[j]=r.sorted_regions[reg]
# show blobs we created
maps.plot_choropleth(shp_link, regions, type='quantiles',
    title='Birth blobs (min '+str(r.floor)+' births per region, '+
    str(r.p)+' regions)', k=r.p)
# best if you full-screen both maps and flip between them
# you can see how blobs create regions that are maximally homogeneous internally
# and maximally heterogeneous to other regions
# blobs also creates as many regions as possible given the 12k births constraint
# some regions are a single county, some are 12 or 13


#####################
## blobs live demo ##
#####################

# what we just did is basically intelligent clustering
# more interesting: can we create homogeneous regions based on some objective
# function (trend in SIDS rate) such that all the regions still meet that 12k 
# births constraint, to avoid small sample bias?

# PySAL's maxP implementation can be subject to randomness - it chooses one initial
# solution and builds incrementally from there. so to ensure we have a good
# solution we run it 50 times, and pick the solution with the smallest weighted
# sum of squares (sum of heterogeneity in SIDS rates within all regions)
best_wss = 10**8
best_solution = r
for i in range(0,50):
    r=ps.Maxp(w, sidr, floor=12000, floor_variable=births74, initial=200)
    if (r.objective_function() < best_wss):
        best_wss = r.objective_function()
        best_solution = r
    print('iteration '+str(i+1)+' - wss: '+str(round(r.objective_function(),2))+
    ' ('+str(r.k)+' regions), best: '+str(round(best_wss,2)))

# plot best solution
r=best_solution
fips=np.array(dbf.by_col('FIPS'))
r.sort_regions()
regions=empty(100)
for j in range(0,100):
    reg=r.area2region[fips[j]]
    regions[j]=r.sorted_regions[reg]
# show SIDS rate in 1974
maps.plot_choropleth(shp_link, numeric_cols[:,12], type='quantiles',
    title='1974 SIDS rate (per 1,000 live births)', k=r.p)
# show blobs based on SIDS rate trend from 1974-1979
maps.plot_choropleth(shp_link, regions, type='quantiles',
    title='SIDS blobs (min '+str(r.floor)+' births per region, '+
    str(r.p)+' regions)', k=r.p)
    
# it's very easy to change the minimum region size, or to change the mechanics
# of the objective function (which defines homogeneity)
# our challenges are primarily computational - ensuring the best feasible solution

