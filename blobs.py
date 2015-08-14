### blobs code, v0.9
### contributors: 
### Jonathan Giuffrida: jgiuffrida@uchicago.edu
### Jiajun Shen: jiajun@uchicago.edu
### 8/14/15

import pysal as ps
import numpy as np
import pandas as pd
from pysal.contrib.viz import mapping as maps
import matplotlib.pyplot as plt
import time
import datetime
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
import sys
from shapely.geometry import mapping, Polygon
import Polygon as pl
import fiona

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

# extend Maxp with new method
ps.Maxp.sort_regions = sort_regions



# class for blobs data
class Blobs_Data:
    """ This pulls and preps the data for blobs from the http://plenar.io API
    for the given datasets and unit of analysis

    Parameters
    ----------
    census_data : string
                  relative file path and name for a CSV holding block-level
                    census data, with FIPS ID and population by block (for
                    an example, see the Chicago Census.csv file in this folder)

    level       : {'tract', 'block group', 'block'}
                  the desired unit of analysis (should match the census_data
                    and shp that you give)
    
    shp         : string
                  relative file path and name for the .shp file you want to use
                    (note: there should also be an associated .gal weights file
                    and .dbf data file with the same name and location; see
                    PySAL documentation for details on these formats)

    shp_id      : string
                  the name of the unique ID in the associated .dbf file

    datasets    : array
                  a list of the machine names of the datasets you want to use
                    from Plenario, e.g. 
                    ['crimes_2001_to_present', 'business_licenses']
                  if empty, will automatically include all available datasets

    temporal_agg: {'day', 'week', 'month', 'quarter', 'year', 'decade'}
                  the desired level of temporal aggregation for the Plenario
                    data (currently does not make a difference in the end)

    time_start  : 'yyyy-mm-dd'
                  the start date for Plenario data; default is Jan 1, 2000

    time_end    : 'yyyy-mm-dd'
                  the end date for Plenario data; default is today

    Attributes
    ----------
    data        : pandas DataFrame
                  contains Census IDs at the appropriate unit of analysis,
                    population, and a count of observations for each dataset 
                    by unit of analysis

    Sample usage
    ------------

    >>> d = Blobs_Data('Chicago Census.csv', 'block', 
          'blocks/CensusBlockTIGER2010.shp', 'geoid10', 
          ['crimes_2001_to_present', '311_service_requests_rodent_baiting'])

    """

    def __init__(self, census_data, level, shp, shp_id, datasets=[], temporal_agg='month', 
        time_start='2000-01-01', time_end=None):
        self.shp_link = shp
        if not time_end:
            time_end = time.strftime('%Y') + '-' + time.strftime('%m') + \
            '-' + time.strftime('%d')
        self.prefix_url = 'http://plenar.io/v1/api/timeseries/?obs_date__ge='+\
            time_start + '&obs_date__le=' + time_end + '&agg=' + \
            temporal_agg + '&data_type=csv'
        if len(datasets) > 0:
            self.prefix_url += '&dataset_name__in=' + ','.join(datasets)

        # data preparation
        census = pd.read_csv(census_data, dtype=object)
        final = []
        done = set([])

        sys.stdout.write('\n')
        sys.stdout.flush()

        for i in range(0, len(census)):
            # assemble the various IDs using the FIPS code
            block = census.id[i][9:]
            if level == 'tract':
                check = block[0:11]
            elif level == 'block group':
                check = block[0:12]
            elif level == 'block':
                check = block[0:]
            else:
                print "error: 'level' must be in {'tract', 'block group', 'block'}"
                return False
            if check not in done:
                done.add(check)
                row = [check, block[0:2], block[2:5], block[5:11]]
                try:
                    row.append(int(census['pop'][i]))
                except: #some pop data is screwed up
                    row.append(0)
                final.append(row)
                sys.stdout.write('\rpreparing line ' + str(i+1) + ' of ' + str(census.shape[0]) + 
                    ' (' + str(round(i * 100./census.shape[0], 2)) + '% done)')
                sys.stdout.flush()
            else:
                # just add the population
                try:
                    final[-1][4] += int(census['pop'][i])
                except:
                    pass

        sys.stdout.write('\rdata preparation complete' + (' ' * 30) + '\n')
        sys.stdout.flush()

        final = pd.DataFrame(np.array(final))

        for d in datasets:
            final[d] = 0

        final.columns = ['ID', 'stateID', 'countyID', 'tractID', 'pop'] + datasets

        times = []
        sys.stdout.write('\rdownloading data...')
        sys.stdout.flush()
        for t in range(0, len(final)):
            start = time.time()
            if level == 'tract' or level == 'block group':
                url = self.prefix_url + '&census_block__ilike=' + str(final.ix[t, 0]) + '%'
            elif level == 'block':
                url = self.prefix_url + '&census_block=' + str(final.ix[t, 0])
            cr = pd.read_csv(url)
            for name in datasets:
                try:
                    final.ix[t, name] = cr[name].sum()
                except:
                    final.ix[t, name] = 0
            end = time.time()
            times.append(end - start)
            sys.stdout.write('\rdownloading data for ' + level + ' ' + str(t+1) + ' of ' + 
                str(final.shape[0]) + ' (' + 
                str(round(t * 100./final.shape[0], 2)) + '% done, ' + 
                str(int(np.mean(times)*(len(final)-t-1)/60.))+' minutes remaining)')
            sys.stdout.flush()
        sys.stdout.write('\rdata download complete' + (' ' * 30) + '\n')
        sys.stdout.flush()

        final.to_csv('plenario data by block.csv', index=False)

        self.dbf = ps.open(self.shp_link[:-3] + 'dbf')
        self.w = ps.open(self.shp_link[:-3] + 'gal').read()
        self.id = shp_id
        self.level = level
        cols = np.array([dbf.by_col(col) for col in dbf.header]).T
        df = pd.DataFrame(cols)
        df.columns = dbf.header
        df.columns = df.columns.map(lambda x: x.lower())
        df['order'] = df.index # mark the order of the shapes
        for c in final.columns[4:]:
            final[c] = final[c].astype('float')

        # merge on shapefile IDs
        if level == 'tract':
            ordered = pd.DataFrame(df.loc[:,[shp_id, 'order']])
            self.data = pd.merge(final, ordered, how='right', left_on='tractID', 
                right_on=shp_id, sort=False).fillna(0).sort(['order'])
        elif level == 'block':
            ordered = pd.DataFrame(df.loc[:,[shp_id, 'order']])
            self.data = pd.merge(final, ordered, how='right', left_on='ID', 
                right_on=shp_id, sort=False).fillna(0).sort(['order'])
        self.data = self.data.drop(['order'])
        print('\rdata ready to use\n\n')
        return True


# main blobs class
class Blobs:
    """Create a max-p regions solution for a given shapefile and associated 
    dataset. Builds on pysal.Maxp with improvements to the user interface, 
    flexibility, and mapping. 

    Original solution from "The Max-p-Regions Problem," Duque, Anselin, and Rey, 
    JRS, October 2010, available at http://geography.sdsu.edu/Research/
    Projects/IPC/publication/MaxP_authored.pdf

    Parameters
    ----------
    bd          : Blobs_Data

    floor_var   : variable to use for the floor (minimum size of the blobs)
                  can be any variable in the dataset, or 'areas' to use the 
                  number of areas

    floor       : minimum size of each blob, as measured by floor_var
                  if floor_var is 'areas', this is the minimum number of areas
                    in each blob

    vars_to_use : variables on which to create blobs
                  default is all variables in the dataset, except for ID ones
                    and population

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
    verbose     : boolean
                  will print out comprehensive information about the progress
                    of the solution


    Attributes
    ----------

    regions     : numpy array
                  The assigned blob for each area, in the original order

    r           : PySAL MaxP instance
                  The best solution found.
    
    contours    : List of Shapely.Polygon instances
                  This is the contours of all the blobs. Notice that the number of contours might be larger
                  than the number of blobs since some blobs contain two contours ("Island situation")l



    plot_blobs(variable=None, k=None): method
                  Will plot the blobs. If a variable name is entered, will plot
                    blobs along that variable. If k is entered, it will be
                    the number of buckets used in the choropleth. Example:

                    b.plot_blobs('pop', 5)

                    Will plot blobs by population into five buckets. 

    build_blobs(): Will rebuild blobs. Currently there is no way to reassign
                     parameters here - to do that, create another Blobs object.
                     However, this method may return a different solution 
                     because of randomness. 

    generate_contours(): Generate contours for blobs after building blobs.

    generate_shpfile(filename): Generate shape file based on the blobs.

    Sample usage
    ------------

    >>> d = blobs.Blobs_Data('Chicago Census.csv', 'block', 'blocks/CensusBlockTIGER2010.shp', 
          'geoid10', ['crimes_2001_to_present', 
          '311_service_requests_vacant_and_abandoned_building', 
          '311_service_requests_rodent_baiting'])
    >>> b = blobs.Blobs(d, 'pop', 10000)

    """
    def __init__(self, bd, floor_var, floor, vars_to_use=[], iterations=10, 
    method='equal votes', weights=[], initial=10, plot=True, savedata=False, 
    plot_values=False, verbose=False):
        self.d = bd.data
        self.w = bd.w
        self.shp_link = bd.shp_link
        self.level = bd.level
        self.floor_var = floor_var
        self.floor = floor
        self.id_var = bd.id
        self.vars_to_use = vars_to_use
        if self.vars_to_use == []:
            self.vars_to_use = [v for v in self.d.columns if v not in \
            ['ID', 'stateID', 'countyID', 'tractID', 'pop', bd.id]]
        print self.vars_to_use
        self.iterations = iterations
        self.method = method
        self.weights = weights
        self.initial = initial
        self.plot = plot
        self.savedata = savedata
        self.plot_values = plot_values
        self.verbose = verbose
        self.r = None
        self.regions = None
        self.blobs_data = None
        self.build_blobs()
        self.generate_contours()

    def _get_floor_var(self):
        return self.d['pop']

    def build_blobs(self):
        """ Method to create a blobs solution.
        """
        solutions = []
        top_scores = []
        times = []
        num_blobs = []
        current_time = []
        iteration = []
        best_score = -1
        best_solution = None
        if self.floor_var == 'areas':
            floor_var_array = np.ones((self.d.shape[0], 1))
        else:
            floor_var_array = self.d[self.floor_var]
        blob_vars = np.array(self.d.loc[:, self.vars_to_use], np.float64)
        
        if len(self.vars_to_use) == 1:
            # add shape to the array
            blob_vars.shape = (blob_vars.shape[0], 1)
        print('\n### CREATING BLOBS FROM ' + str(len(self.vars_to_use)) + 
            ' VARIABLES ###\n    PARAMETERS:\n     # Minimum ' + self.floor_var + ' in each blob: ' + 
            str(int(self.floor)) + '\n     # Iterations: ' + str(self.iterations) +
            '\n     # Method: ' + self.method + '\n     # Plot blobs: ' + str(self.plot) + 
            '\n     # Save blobs data: ' + str(self.savedata) + '\n')

        for i in range(0,self.iterations):
            start = time.time()
            r=ps.Maxp(self.w, self._format_blobs(blob_vars),
                floor=self.floor, floor_variable=floor_var_array, 
                initial=self.initial, verbose=self.verbose)
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
            msg = '\n# ITERATION '+str(i+1)+'                 \n  Score: ' + \
                str(round(current_score,2)) + '\n  Created '+str(r.k)+' blobs (' + \
                str(int(self.d.shape[0]/r.k)) + ' tracts per blob)\n  Best solution so far: ' + \
                str(round(best_score,2))
            msg += '\n  Time taken: '+str(round(end-start,1))+' seconds ('+ \
                str(int(np.mean(times)*(self.iterations-i-1)))+' seconds remaining)\n'
            print msg
        
        r = best_solution
        print('\r# BEST SOLUTION:                      \n  Score: '+
            str(round(r.objective_function(),2)) + 
            '\n  '+str(r.k)+' blobs ('+str(int(self.d.shape[0]/r.k))+
            ' tracts per blob)')
        self.r = r
        # prep for plotting
        ids=np.array(self.d['tractce10']).astype(str)
        if self.plot_values:
            self.r.sort_regions(method='mean')  # sort regions by intensity of the variable
        regions=np.empty(self.d.shape[0])
        for j in range(0,self.d.shape[0]):
            reg=r.area2region[ids[j]]
            regions[j]=reg
        self.regions = regions
        if self.plot:
            self.plot_blobs()
        self.build_data_structure(self.savedata)

    # helper function to assign weights to variables
    def _format_blobs(self, data):
        if self.method == 'default':
            # use max p as originally designed
            # variables will be implicitly weighted in proportion to their means
            return data
        elif self.method == 'equal votes':
            # give equal weight to all variables by standardizing them
            x = np.zeros(data.shape)
            for v in range(data.shape[1]):
                x[:,v] = (data[:,v] - np.mean(data[:,v])) / np.std(data[:,v])
            return x
        elif self.method == 'weighted':
            # assign explicit weights to standardized variables
            x = np.zeros(data.shape)
            for v in range(data.shape[1]):
                x[:,v] = ((data[:,v] - np.mean(data[:,v])) / \
                    np.std(data[:,v])) * np.sqrt(self.weights[v])
            return x

    def plot_blobs(self, blob_shp=None, variable=None, k=None, mapType=None):
        # show blobs we created
        if blob_shp:
            k = len(self.contours)
            if not variable:
                variable = 'blob ID'
                data = np.arange(len(self.contours))
                mapType = 'unique_values'
            else:
                data = []
                try:
                    a = self.blobs_data[variable]
                except:
                    try:
                        a = self.blobs_data[variable + '_mean']
                        variable = variable + '_mean'
                    except:
                        print("invalid variable name; variables are the following: \n" + 
                            ', '.join(self.blobs_data.columns))
                        return False
                for i in range(len(self.contours)):
                    try:
                        data.append(self.blobs_data.ix[self.contours_to_blobs[i],variable])
                    except KeyError:
                        data.append(0)
                data = np.array(data)
                mapType = 'quantiles'
        elif not variable:
            data = self.regions
            variable = 'blob ID'
        else:
            data = []
            try:
                a = self.blobs_data[variable]
            except:
                try:
                    a = self.blobs_data[variable + '_mean']
                    variable = variable + '_mean'
                except:
                    print("invalid variable name; variables are the following: \n" + 
                        ', '.join(self.blobs_data.columns))
                    return False
            for i in self.d.tractID:
                try:
                    data.append(self.blobs_data.ix[self.r.area2region[str(i)],variable])
                except KeyError:
                    data.append(0)
            data = np.array(data)

        if not k:
            k = self.r.p
        if not mapType:
            mapType = 'quantiles'
        print('  Plotting...')
        map_shp = None
        
        if blob_shp:
            map_shp = blob_shp
        else:
            map_shp = self.shp_link
        maps.plot_choropleth(map_shp, data, type=mapType,
            title='Blobs from Census ' + self.level + 's\nby ' + variable + 
                ' (' + str(self.r.p)+' blobs)', k=k, figsize=(6,8))
        print('\r             \n')

    def build_data_structure(self, savedata=True):
        #build data structure
        sr = np.zeros([self.r.k, len(self.vars_to_use)*2+4])
        for region in range(0,self.r.k):
            # blob ID
            sr[region][0] = region
            selectionIDs = [self.r.w.id_order.index(i) for i in self.r.regions[region]]
            m = self.r.z[selectionIDs, :]
            # objective function
            var = m.var(axis=0)
            sr[region][1] = sum(np.transpose(var)) * len(self.r.regions[region])
            # blob size (number of places in blob)
            sr[region][2] = len(self.r.regions[region])
            # blob population
            if self.floor_var == 'areas':
                sr[region][3] = len(self.r.regions[region])
            else:    
                sr[region][3] = self.d.loc[selectionIDs, self.floor_var].sum()
            # variable means and standard deviations
            for j in range(0,len(self.vars_to_use)):
                sr[region][4+j*2] = m[:,j].mean()
                sr[region][5+j*2] = m[:,j].std()
        srdf = pd.DataFrame(sr)
        cols = ['Blob', 'Score', 'Number of Regions', self.floor_var]
        for j in range(0, len(self.vars_to_use)):
            cols.append(self.vars_to_use[j]+'_mean')
            cols.append(self.vars_to_use[j]+'_stdev')
        srdf.columns = cols
        if savedata:
            srdf.to_csv('Blobs data ' + datetime.datetime.now().strftime('%Y%m%d %H%M') + \
                '.csv', index=False)
        self.blobs_data = srdf

    # helper function to retrieve original, non-standardized data
    def retrieve_raw_data(self):
        pass  # todo

    def generate_contours(self):
        allPoly = ps.open(self.shp_link)
        blobPoly = [None for i in range(len(np.unique(self.regions)))]
        for i in range(len(allPoly)):
            if(blobPoly[int(self.regions[i])] == None):
                blobPoly[int(self.regions[i])] = pl.Polygon(allPoly[i].vertices)
            else:
                blobPoly[int(self.regions[i])] = blobPoly[int(self.regions[i])] + pl.Polygon(allPoly[i].vertices)
        outputPoly = []
        contours_to_blobs = []
        counter = 0
        for poly in blobPoly:
            for i in range(len(poly)):
                outputPoly.append(Polygon(poly.contour(i)))
                contours_to_blobs.append(counter)
            counter+=1
        self.contours = outputPoly
        self.contours_to_blobs = contours_to_blobs
    
    def generate_shapefile(self, filename='./blob_shapefile.shp'):
        self._generate_shapefile(self.contours, filename)

    def _generate_shapefile(self, polygons, filename):
        """Generate a shape file given a list of polygons.""" 
        # Define a polygon feature geometry with one attribute
        schema = {
            'geometry': 'Polygon',
            'properties': {'id': 'int'}
        }
        
        # write a new Shapefile
        with fiona.open(filename, 'w', 'ESRI Shapefile', schema) as c:
            for i in range(len(polygons)):
                c.write({
                    'geometry': mapping(polygons[i]),
                    'properties': {'id': i},
                }) 
        

# cluster the blobs data (k-means)
class Cluster_blobs:
    """Use k-means to cluster blobs along the explanatory variables.

    Parameters
    ----------

    b               : Blobs object
                      should be a Blobs object, accessible by calling Blobs()

    variables       : array
                      an array of variable names to use for clustering; by 
                        default, all variables ending in "_mean" will be used

    n_clusters      : int
                      (optional) number of clusters to form

    blobs_per_      : int
      cluster         (optional) average number of blobs per cluster. if both
                        n_clusters and blobs_per_cluster have values, the 
                        former will be ignored

    Attributes
    ----------

    assignments     : numpy array
                      an n*1 array of cluster labels, in order

    centers         : pandas DataFrame
                      the coordinates for the cluster centers, in order

    inertia         : float
                      the "inertia" for the final solution; lower is better

    set_n_clusters(n_clusters): Will use n_clusters clusters and re-cluster.

    set_blobs_per_cluster(blobs_per_cluster): Will use blobs_per_cluster and 
                      re-cluster.

    plot(variables=[]): Will plot the clusters. If a list of variables is 
                      supplied, it will plot along the first three supplied. 
                      The default list is all variables used to cluster. 

    plot_map(variable=None): Will plot the clusters on a map. If a variable
                      is supplied, the clusters will be colored by their 
                      values along that variable. 

    Sample usage
    ------------

    >>> solution = Blobs(dataset, min_pop=10000)
    >>> cl = Cluster_blobs(solution.data, blobs_per_cluster=10)
    >>> print(cl.centers)

    """

    def __init__(self, b, variables=[], n_clusters=0, blobs_per_cluster=0):
        """Initialize, run k-means, and plot."""
        # build list of variables on which to cluster (if not provided)
        self.b = b
        if variables == []:
            self.cluster_vars = []
            for c in self.b.blobs_data.columns:
                if c.find('_mean') > 0:
                    self.cluster_vars.append(c)
        else:
            self.cluster_vars = variables

        self.x = np.array(self.b.blobs_data[self.cluster_vars])
        self.n_clusters = n_clusters
        self.blobs_per_cluster = blobs_per_cluster
        self._set_clusters()

        self.assignments = None
        self.centers = None
        self.inertia = -1
        self.clusters2blobs = {}
        self.kmeans()
        self.plot()

    def _set_clusters(self):
        """Recalculate number of clusters."""
        # set n_clusters
        if self.n_clusters == 0 and self.blobs_per_cluster == 0:
            self.blobs_per_cluster = 10  # totally arbitrary
            self.n_clusters = int(np.round(self.x.shape[0] / 
                self.blobs_per_cluster))
        elif self.blobs_per_cluster:
            self.n_clusters = int(np.round(self.x.shape[0] / self.blobs_per_cluster))
        elif self.n_clusters:
            self.blobs_per_cluster = int(np.round(self.x.shape[0] / self.n_clusters))

    def kmeans(self):
        """Run k-means with current settings."""
        # run k-means with all available CPU cores
        self.est = KMeans(n_clusters=self.n_clusters, n_jobs=-1) 
        # the preceding can throw an error in multiprocessing.py in pycharm; change to n_jobs=1 to fix
        self.est.fit(self.x)
        # save a lot of data to work with, save, etc.
        self.assignments = self.est.labels_ 
        self.centers = pd.DataFrame(self.est.cluster_centers_, columns=self.cluster_vars)
        self.inertia = self.est.inertia_
        self.clusters2blobs = {}
        for c in np.unique(self.assignments):
            self.clusters2blobs[str(c)] = []
        for i in range(0, len(self.x)):
            self.clusters2blobs[str(self.assignments[i])].append(str(int(self.x[i][0])))


    def plot(self, variables=[]):
        """Plot the most recent k-means solution.

        Parameters
        ----------

        variables   : array
                      list of variables to plot (by default, uses variables 
                        used in kmeans). if only two are provided, will plot in 2D.

        """
        fig = plt.figure(figsize=(10, 9))
        plt.clf()

        if len(variables) == 0:
            vars_to_plot = self.cluster_vars
        else:
            vars_to_plot = variables

        pos = {}

        
        if len(vars_to_plot) == 1:
            print('No graph shown because only one variable was used')
            return True
        elif len(vars_to_plot) == 2:
            ax = fig.add_subplot(111)
        else:
            # will automatically plot in 3D if 3+ variables
            ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=30, azim=134)


        plt.cla()
        labels = self.est.labels_

        # click event to print data about each point on user interaction
        def onpick(event):
            ind = event.ind
            for i in ind:
                type = event.artist.get_label()
                msg = ''
                if type == 'Blobs':
                    msg = 'Blob ' + str(i)
                    for c in range(self.x.shape[1]):
                        msg += '\n  ' + self.cluster_vars[c][:self.cluster_vars[c].find('_mean')]+\
                                ': ' + str(round(np.take(self.x[:,c], i), 2))
                    msg += '\n  (All values are z-scores)'
                    neighbors = np.where(self.est.labels_ == self.est.labels_[i])[0]
                    if len(neighbors) > 1:
                        msg += '\n  Other blobs in cluster: ' + \
                                ', '.join([k for k in neighbors.astype('str') if not k==str(i)]) + \
                                '\n'
                elif type == 'Clusters':
                    msg = 'Cluster ' + str(i)
                    msg += '\n  Center of cluster (all values in z-scores):'
                    for c in range(self.est.cluster_centers_.shape[1]):
                        msg += '\n  ' + self.cluster_vars[c][:self.cluster_vars[c].find('_mean')]+\
                                ': ' + str(round(np.take(self.est.cluster_centers_[:,c], i), 2))
                    inhabitants = np.where(self.est.labels_ == i)[0]
                    msg += '\n  Blobs in cluster: ' + \
                            ', '.join([k for k in inhabitants.astype('str') if not k==str(i)])+'\n'
                print msg

        # prepare scatterplots and axes
        if len(vars_to_plot) == 2:
            ax.scatter(x[:, self.cluster_vars.index(vars_to_plot[0])], 
                x[:, self.cluster_vars.index(vars_to_plot[1])], s=30, 
                c=labels.astype(np.float),label="Blobs",picker=True)
            ax.scatter(self.est.cluster_centers_[:,self.cluster_vars.index(vars_to_plot[0])], 
                self.est.cluster_centers_[:,self.cluster_vars.index(vars_to_plot[1])], s=40, 
                marker='*', c=range(self.est.n_clusters), label="Clusters", picker=True)
            ax.xaxis.set_ticklabels(ax.xaxis.get_ticklocs())
            ax.yaxis.set_ticklabels(ax.yaxis.get_ticklocs())
            ax.set_xlabel(vars_to_plot[0])
            ax.set_ylabel(vars_to_plot[1])
        elif len(vars_to_plot) > 2:
            ax.scatter(self.x[:, self.cluster_vars.index(vars_to_plot[0])], 
                self.x[:, self.cluster_vars.index(vars_to_plot[1])], 
                self.x[:, self.cluster_vars.index(vars_to_plot[2])], s=30, 
                c=labels.astype(np.float), label="Blobs", picker=True)
            ax.scatter(self.est.cluster_centers_[:,self.cluster_vars.index(vars_to_plot[0])], 
                self.est.cluster_centers_[:,self.cluster_vars.index(vars_to_plot[1])], 
                self.est.cluster_centers_[:,self.cluster_vars.index(vars_to_plot[2])], 
                s=40, marker='*', c=range(self.est.n_clusters), label="Clusters", picker=True)
            ax.w_xaxis.set_ticklabels(ax.w_xaxis.get_ticklocs())
            ax.w_yaxis.set_ticklabels(ax.w_yaxis.get_ticklocs())
            ax.w_zaxis.set_ticklabels(ax.w_zaxis.get_ticklocs())
            ax.set_xlabel(vars_to_plot[0])
            ax.set_ylabel(vars_to_plot[1])
            ax.set_zlabel(vars_to_plot[2])
        ax.set_axisbelow(True)
        fig.canvas.mpl_connect('pick_event', onpick)

        plt.title(str(self.x.shape[0]) + ' Blobs in ' + str(self.est.n_clusters) + 
            ' Clusters (Based on ' + str(len(vars_to_plot)) + 
                ' variables)\nClick on values for more information')
        plt.show()


    def set_n_clusters(self, n_clusters):
        """Set the desired number of clusters."""
        # reset number of clusters
        if int(n_clusters) > 0:
            self.n_clusters = n_clusters
            self.blobs_per_cluster = 0
            self._set_clusters()
        else:
            print("Error: please provide n_clusters as an int")

    def set_blobs_per_cluster(self, blobs_per_cluster):
        """Set the desired number of blobs per cluster."""
        # reset number of blobs per cluster
        if int(blobs_per_cluster) > 0:
            self.blobs_per_cluster = blobs_per_cluster
            self.n_clusters = 0
            self._set_clusters()
        else:
            print("Error: please provide blobs_per_cluster as an int")

    def plot_map(self, variable=None, blob_shp=None):
        # plot clusters on the map
        plots = np.zeros(self.b.regions.shape)
        if blob_shp:
            cluster = np.zeros(len(self.b.contours))
            for i in range(len(self.b.contours)):
                cluster[i] = self.assignments[self.b.contours_to_blobs[i]]
            maps.plot_choropleth(blob_shp, cluster, type='unique_values',
                title=('Clustered blobs from Census ' + self.b.level + 's'), 
                k=30, figsize=(6,8))
        elif not variable:
            for i in range(len(self.b.regions)):
                plots[i] = self.assignments[self.b.regions[i]]
            maps.plot_choropleth(self.b.shp_link, plots, type='unique_values',
                title=('Clustered blobs from Census ' + self.b.level + 's'), 
                k=30, figsize=(6,8))
        else:
            for i in range(len(self.b.regions)):
                plots[i] = self.centers.ix[self.assignments[self.b.regions[i]], variable]
            maps.plot_choropleth(self.b.shp_link, plots, type='equal_interval',
                title=('Clustered blobs from Census ' + self.b.level + 's'), 
                k=30, figsize=(6,9))


