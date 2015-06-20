Blobs: Smart Clustering at the Urban Level
=====

*v0.9: June 2015*

Spatial analytics is often hampered by the arbitrary choice of units, allowing local heterogeneity to obscure true patterns. Blobs is a new “smart clustering” technique that lets us use large quantities of open municipal data from [Plenario](http://plenar.io) to redraw city maps to reflect facts on the ground, not administrative boundaries. 

The algorithm, built on the *max-p regions* implementation by researchers at Arizona State University, creates spatial clusters using only one input parameter, the minimum size of each cluster (defined in any of several ways). This nonparametric approach creates clusters that fit the data as closely as possible, fully isolating regions based only on the variables of interest.  

This package allows a user to create "blobs" from start to finish using any dataset in Plenario. This involves the following:

* Downloading the relevant data from Plenario, at the desired unit of analysis (currently supports Census blocks and tracts)
* Building blobs from those units of analysis using counts of observations from any combination of the requested datasets
* Map the blobs solution along any variable
* Cluster blobs using k-means
* Save data on blobs to use in any research or administrative solution

More information on the basic approach can be found in the [documentation for the PySAL Maxp class](http://www.pysal.org/library/region/maxp.html), upon which blobs is built. 


Dependencies
------------

* pysal
* numpy
* pandas
* matplotlib.pyplot with mpl_toolkits.mplot3d
* sklearn.cluster

Repository Structure
------------

##### Primary

* `blobs.py` - main module
* `maxp.py` - Updated version of the maxp.py module in pysal/region 
* `smoothing.py` - adventures in spatial autocorrelation

##### Secondary (Samples from Chicago)

* `blocks/` - Census block shapefiles (Chicago only)
* `tracts/` - Census tract shapefiles (Chicago only)
* `configure.py` - example code for building the shapefiles and associated files using Census data
* `Chicago Census.csv` - Census IDs and population data for blocks in Chicago
* `master311.csv` - sample Plenario data on Chicago 311 calls by Census tract

  
Example Usage
----------------------
######  Using Chicago data, included

```python
import blobs

############
## tract level
############

# download the data
d = blobs.Blobs_Data('Chicago Census.csv', 'tract', 'tracts/CensusTractsTIGER2010.shp', 
  'tractce10', ['crimes_2001_to_present', 
  '311_service_requests_vacant_and_abandoned_building', 
  '311_service_requests_rodent_baiting'])

# create blobs (minimum population of 10,000 in each blob)
b = blobs.Blobs(d, 'pop', 10000)

# cluster the blobs along similarities in the data
cl = blobs.Cluster_blobs(b.blobs_data, blobs_per_cluster=10)

# create blobs with a minimum of 30 tracts in each blob, and cluster
b = blobs.Blobs(d, 'areas', 30)
cl = blobs.Cluster_blobs(b, blobs_per_cluster=10)

# have around 3 blobs per cluster
cl.set_n_clusters(3)

# see blob assignments
b.assignments

# see cluster assignments
cl.clusters2blobs

# plot blobs along one of the variables
b.plot_blobs('crimes_2001_to_present')


############
## block level
############

d = blobs.Blobs_Data('Chicago Census.csv', 'block', 'blocks/CensusBlockTIGER2010.shp', 
  'geoid10', ['crimes_2001_to_present', 
  '311_service_requests_vacant_and_abandoned_building', 
  '311_service_requests_rodent_baiting'])
b = blobs.Blobs(d, 'pop', 10000)
cl = blobs.Cluster_blobs(b, blobs_per_cluster=10)



###########
## instant examples
###########

# if you don't want to download data now, you can test out blobs using some included data
# just run the following

import numpy as np
import pandas as pd
import pysal as ps
shp_link = 'tracts/CensusTractsTIGER2010.shp'
dbf = ps.open('tracts/CensusTractsTIGER2010.dbf')
cols = np.array([dbf.by_col(col) for col in dbf.header]).T
df = pd.DataFrame(cols)
df.columns = dbf.header
df.columns = df.columns.map(lambda x: x.lower())
df.commarea = df.commarea.astype('int')
df['order'] = df.index
w=ps.open('tracts/CensusTractsTIGER2010.gal').read()
init_calls = pd.read_csv('master311.csv', dtype=object)
for c in init_calls.columns[1:]:
    init_calls[c] = init_calls[c].astype('float')

# format data and merge on shapefile IDs
ordered_tracts = pd.DataFrame(df.loc[:,['tractce10', 'commarea', 'order']])
calls = pd.merge(init_calls, ordered_tracts, how='right', left_on='tractID', 
    right_on='tractce10', sort=False).fillna(0).sort(['order'])
calls = calls.drop(['order', 'commarea'],1)

class bd:
  data = calls
  w = w
  shp_link = shp_link
  id = 'tractce10'
  level = 'tract'

d = bd()
b = blobs.Blobs(d, 'pop', 10000, iterations=1)
cl = blobs.Cluster_blobs(b, blobs_per_cluster=10)


```



