# Import required libraries and assign root directory of blobs to root
root = "/home/jiajun/Desktop/blobs"
import os
os.chdir(root)
print os.getcwd()
import blobs
import numpy as np
import pandas as pd
import pysal as ps
import matplotlib as mpl

shp_link = root + '/blocks/CensusBlockTIGER2010.shp'
dbf = ps.open(root + '/blocks/CensusBlockTIGER2010.dbf')

from pysal.contrib.viz import mapping as maps

# Shapefile and data preparation.
cols = np.array([dbf.by_col(col) for col in dbf.header]).T
df = pd.DataFrame(cols)
df.columns = dbf.header
df.columns = df.columns.map(lambda x: x.lower())
df['order'] = df.index

df_pop = pd.read_csv(root + "/blocks/censusblockPop.csv",dtype = "object")
df_pop['Pop'] = df_pop['Pop'].astype('int')
df = pd.merge(df_pop, df, how='right', left_on='CENSUS BLOCK', 
    right_on='tract_bloc', sort=False).fillna(0).sort(['order'])
df.sort(columns = 'order')
# Assign spatial weight for census tracts.
w=ps.open(root + '/blocks/CensusBlockTIGER2010.gal').read()
crimes = pd.read_csv(root + '/usecase_temporalSlicing/intermediate_result.csv', dtype = object)
for c in crimes.columns[1:]:
    crimes[c] = crimes[c].astype('int')
crime_hour = []
for i in range(24):
    crime_hour.append(crimes[crimes['hour'] == 0])

# Format data and merge on shapefile IDs
ordered_blocks = pd.DataFrame(df.loc[:,['geoid10', 'tract_bloc','Pop', 'order']])
ordered_blocks.sort(['order'])
crime_hour_tables = []
for i in range(24):
    hourlyTable = pd.merge(crime_hour[i],ordered_blocks, how = 'right', left_on = 'geoid10',
                             right_on = 'geoid10', sort=False).fillna(0).sort(['order'])
    hourlyTable = hourlyTable.drop(['order', 'hour'],1)
    crime_hour_tables.append(hourlyTable)
i = 0
print i
class bd:
  data = crime_hour_tables[i]
  w = w
  shp_link = shp_link
  id = 'geoid10'
  level = 'block'
d = bd()
b = blobs.Blobs(d, 'Pop', 1000, plot=False, iterations=1)
b.generate_shapefile(filename='blob_hour_test_0_1000.shp')
