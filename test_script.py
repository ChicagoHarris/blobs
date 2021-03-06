### script to test blobs
### contributors: jgiuffrida@uchicago.edu
### 6/18/15

import blobs

# tract level
d = blobs.Blobs_Data('Chicago Census.csv', 'tract', 'tracts/CensusTractsTIGER2010.shp', 
  'tractce10', ['crimes_2001_to_present', 
  '311_service_requests_vacant_and_abandoned_building', 
  '311_service_requests_rodent_baiting'])
b = blobs.Blobs(d, 'pop', 10000)

# block level
d = blobs.Blobs_Data('Chicago Census.csv', 'block', 'blocks/CensusBlockTIGER2010.shp', 
  'geoid10', ['crimes_2001_to_present', 
  '311_service_requests_vacant_and_abandoned_building', 
  '311_service_requests_rodent_baiting'])
b = blobs.Blobs(d, 'pop', 10000)
cl = blobs.Cluster_blobs(b, blobs_per_cluster=10)



# the following will let you get started immediately on pre-downloaded data
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

# example 1: blobs on three sanitation-related variables
b = blobs.Blobs(d, 'pop', 10000, vars_to_use=['sanitation_per1000', 
  'rodents_per1000', 'buildings_per1000'], iterations=1)

cl = blobs.Cluster_blobs(b, blobs_per_cluster=10)
# try clicking on the dots and stars
print cl.centers
# note the cluster that is off-the-charts high in rodents, but not in sanitation/buildings.
# this is ideal for a test-control situation by the sanitation department
print cl.assignments
# plot blobs on a map
cl.plot_map()  # by cluster number
cl.plot_map('sanitation_per1000_mean')  # by average sanitation calls per cluster

# can easily re-run with larger or smaller clusters:
cl = blobs.Cluster_blobs(b, blobs_per_cluster=5)
# can also set number of clusters directly:
cl.set_n_clusters(3)


# example 2: blobs on all variables
b = blobs.Blobs(d, 'pop', 10000, iterations=3)
cl = blobs.Cluster_blobs(b, blobs_per_cluster=20)
# note that the 3D graph cannot show all relevant data
# can request which variables to plot
cl.plot(['rodents_per1000_mean', 'garbage_per1000_mean', 'vehicles_per1000_mean'])



# # compare to an existing solution: community areas
# ca_regions = []
# for i in np.unique(calls['commarea']):
#     ca = []
#     for j in range(len(calls)):
#         if calls.ix[j, 'commarea'] == i:
#             ca.append(calls.ix[j, 'tractce10'])
#     ca_regions.append(ca)

# print('our solution: ' + str(b.r.objective_function()))
# print('existing solution: ' + str(b.r.objective_function(ca_regions)))




# # try using the interface too
# help(interface)

