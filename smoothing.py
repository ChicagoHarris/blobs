# learning spatial autocorrelation

import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
%cd 'Documents/tech/git/blobs'


shp_link = 'tracts/CensusTractsTIGER2010.shp'
dbf = ps.open('tracts/CensusTractsTIGER2010.dbf')
w=ps.open('tracts/CensusTractsTIGER2010_fixed.gal').read()

cols = np.array([dbf.by_col(col) for col in dbf.header]).T
df = pd.DataFrame(cols)
df.columns = dbf.header
df.columns = df.columns.map(lambda x: x.lower())
df.commarea = df.commarea.astype('int')
df['order'] = df.index
calls = pd.read_csv('master311.csv', dtype=object)
for c in calls.columns[1:]:
    calls[c] = calls[c].astype('float')
ordered_tracts = pd.DataFrame(df.loc[:,['tractce10', 'commarea', 'order']])
calls = pd.merge(calls, ordered_tracts, how='right', left_on='tract', 
    right_on='tractce10', sort=False).fillna(0)
calls = calls.sort(['order'])

# all calls by census tract
y = np.array(calls['all_calls_per1000'])
# map values
maps.plot_choropleth(shp_link, np.array(calls.all_calls_per1000), type='fisher_jenks',
     title='All 311 Calls by Census Area, 2011-2015\nUnsmoothed', 
     k=20, figsize=(6,9))

# Global Moran's I
mi = ps.Moran(y, w)
mi.I
mi.EI
mi.p_norm


# Geary's C
gc = ps.Geary(y, w)
gc.C
gc.EC
gc.z_norm


# Local Moran's I
lm = ps.Moran_Local(y, w)
sigs = lm.p_sim
for i in range(lm.n):
  if sigs[i] > 0.1:
    sigs[i] = 0
  elif sigs[i] > 0.05:
    sigs[i] = 1
  elif sigs[i] > 0.01:
    sigs[i] = 2
  elif sigs[i] > 0.001:
    sigs[i] = 3
  elif sigs[i] > 0.0:
    sigs[i] = 4

# plot significant autocorrelation
maps.plot_choropleth(shp_link, np.array(sigs), type='equal_interval',
     title='Significant Spatial Autocorrelation in 311 Calls', k=10, figsize=(6,9))


# spatial smoothing
from pysal.esda import smoothing as sm

pop = np.array(calls['pop'])
all311 = np.array(calls['all_calls'])

# Locally weighted median smoothing, weighted by pop
rate = sm.Spatial_Median_Rate(all311, pop, w, aw=pop)
# weights are populations
y = rate.r
# check for nans
len(y[np.isnan(y)])
y[np.isnan(y)] = 0

maps.plot_choropleth(shp_link, y, type='fisher_jenks',
     title='All Calls by Census Area, 2011-2015\nSpatially Smoothed', 
     k=20, figsize=(6,9))


# Locally weighted median smoothing, weighted by pop, iterated three times
rate = sm.Spatial_Median_Rate(all311, pop, w, aw=pop, iteration=3)
# weights are populations
y = rate.r
# check for nans
len(y[np.isnan(y)])
y[np.isnan(y)] = 0

maps.plot_choropleth(shp_link, y, type='fisher_jenks',
     title='All Calls by Census Area, 2011-2015\nSpatially Smoothed 3 Times', 
     k=20, figsize=(6,9))


lm = ps.Moran_Local(y, w)
sigs = lm.p_sim
for i in range(lm.n):
  if sigs[i] > 0.1:
    sigs[i] = 0
  elif sigs[i] > 0.05:
    sigs[i] = 1
  elif sigs[i] > 0.01:
    sigs[i] = 2
  elif sigs[i] > 0.001:
    sigs[i] = 3
  elif sigs[i] > 0.0:
    sigs[i] = 4

# plot significant autocorrelation
maps.plot_choropleth(shp_link, np.array(sigs), type='equal_interval',
     title='Significant Spatial Autocorrelation in 311 Calls\nAfter Smoothing', 
     k=20, figsize=(6,9))


# Locally weighted median smoothing, weighted by pop, iterated ten times
rate = sm.Spatial_Median_Rate(all311, pop, w, aw=pop, iteration=10)
# weights are populations
y = rate.r
# check for nans
len(y[np.isnan(y)])
y[np.isnan(y)] = 0

maps.plot_choropleth(shp_link, y, type='fisher_jenks',
     title='All Calls by Census Area, 2011-2015\nSpatially Smoothed 10 Times', 
     k=20, figsize=(6,9))

































