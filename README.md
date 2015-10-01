Blobs: Smart Clustering at the Urban Level
=====

*v0.9: Oct 2015*

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
* shapely
* polygon
* fiona

Prerequisite packages can be installed using `install.sh`. You can install the packages by:

    [sudo] bash install.sh

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
We present an example usage in `tutorials/Blob Example.ipynb`

