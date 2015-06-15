# This python file pulls the data for Chicago at a census tract level 
# from plenario's api for the listed 311 call datasets below
import urlparse
import requests
import urllib2
import numpy as np
import pandas as pd
import time
import sys

%cd 'Documents/tech/git/blobs'

datasets=['311_service_requests_abandoned_vehicles',
	'311_service_requests_alley_lights_out', '311_service_requests_tree_trims',
	'311_service_requests_street_lights_all_out', '311_service_requests_pot_holes_reported',
	'311_service_requests_garbage_carts', '311_service_requests_sanitation_code_complaints',
	'311_service_requests_graffiti_removal', '311_service_requests_vacant_and_abandoned_building',
	'311_service_requests_tree_debris', '311_service_requests_rodent_baiting',
	'311_service_requests_street_lights_one_out']
census = pd.read_csv('census tracts.csv', dtype=object)

# construct url
prefix_url = 'http://plenar.io/v1/api/timeseries/?obs_date__ge=2011-01-01&dataset_name='
mid_url = '&agg=month&census_block__ilike=17031'
suffix_url = '%&data_type=csv'

# plenario can technically handle querying all datasets at once using /v1/api/timeseries,
# but in practice it is too slow
for name in datasets:
	census[name] = 0
	t = 0
	for census_tract in census.Tract:
		url = prefix_url + name + mid_url + census_tract + suffix_url
		cr = pd.read_csv(url)
		try:
			census.ix[census.Tract == census_tract, name] = cr[name].sum()
		except:
			census.ix[census.Tract == census_tract, name] = 0
		t += 1
		sys.stdout.write('\r' + name + ': tract ' + str(t) + ' of ' + str(census.shape[0]))
		sys.stdout.flush()
		time.sleep(0.02) # to avoid overwhelming server
	print('\n')

census.to_csv('plenario data.csv', index=False)







# This python file pulls the data for Chicago at a census block level 
# from plenario's api for the listed 311 call datasets below
import urlparse
import requests
import urllib2
import numpy as np
import pandas as pd
import sys

%cd 'Documents/tech/git/blobs'

datasets=['311_service_requests_abandoned_vehicles',
	'311_service_requests_alley_lights_out', '311_service_requests_tree_trims',
	'311_service_requests_street_lights_all_out', '311_service_requests_pot_holes_reported',
	'311_service_requests_garbage_carts', '311_service_requests_sanitation_code_complaints',
	'311_service_requests_graffiti_removal', '311_service_requests_vacant_and_abandoned_building',
	'311_service_requests_tree_debris', '311_service_requests_rodent_baiting',
	'311_service_requests_street_lights_one_out']
census = pd.read_csv('census blocks.csv', dtype=object)

# construct url
prefix_url = 'http://plenar.io/v1/api/timeseries/?obs_date__ge=2011-01-01&dataset_name__in='
prefix_url += ','.join(datasets) + '&agg=month&data_type=csv&census_block='


census['block_id'] = ''
for name in datasets:
	census[name] = 0
t = 0
sys.stdout.write('\n')
for census_block in census.Id[:,]:
	census.ix[t, 'block_id'] = census_block[9:]
	url = prefix_url + census_block[9:]
	cr = pd.read_csv(url)
	for name in datasets:
		try:
			census.ix[t, name] = cr[name].sum()
		except:
			census.ix[t, name] = 0
	t += 1
	sys.stdout.write('\rtract ' + str(t) + ' of ' + str(census.shape[0]) + ' (' + 
		str(round(t * 100./census.shape[0], 2)) + '% done)')
	sys.stdout.flush()

census.to_csv('plenario data by block.csv', index=False)










import requests
import time
import sys
import csv
import os
 
datasets=['311_service_requests_abandoned_vehicles',
    '311_service_requests_alley_lights_out', '311_service_requests_tree_trims',
    '311_service_requests_street_lights_all_out', '311_service_requests_pot_holes_reported',
    '311_service_requests_garbage_carts', '311_service_requests_sanitation_code_complaints',
    '311_service_requests_graffiti_removal', '311_service_requests_vacant_and_abandoned_building',
    '311_service_requests_tree_debris', '311_service_requests_rodent_baiting',
    '311_service_requests_street_lights_one_out']
 
# plenario can technically handle querying all datasets at once using /v1/api/timeseries,
# but in practice it is too slow
tract_info = []
with open('tracts.csv', 'r') as f:
    tract_info = list(csv.DictReader(f))
 
url = 'http://plenar.io/v1/api/timeseries/'
 
totals = []
 
for census_tract in tract_info:
    outp_name = '%s.csv' % census_tract['Tract']
    print('working on %s' % census_tract['Tract'])
    if not os.path.exists(outp_name):
        params = {
            'agg': 'month',
            'obs_date__ge': '2011-01-01',
            'dataset_name__in': ','.join(datasets),
            'census_block__like': '17031%s%%' % census_tract['Tract'],
            'data_type': 'csv'
        }
        resp = requests.get(url, params=params)
        with open(outp_name, 'w') as outp:
            outp.write(resp.content.decode('utf-8'))
 
    with open(outp_name, 'r') as f:
        rows = list(csv.DictReader(f))
        tract_totals = {d: 0 for d in datasets}
        for dataset in datasets:
            try:
                tract_totals[dataset] = sum([int(r[dataset]) for r in rows])
            except KeyError:
                tract_totals[dataset] = 0
        outp_row = {
            'tract': census_tract['Tract'],
            'population': census_tract['Population'],
            'community_area': census_tract['CA']
        }
        outp_row.update(tract_totals)
        totals.append(outp_row)
 
with open('totals.csv', 'w') as f:
    writer = csv.DictWriter(f, fieldnames=totals[0].keys())
    writer.writeheader()
    writer.writerows(totals)









