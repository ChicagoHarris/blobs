blobs
=====

usage:
```
interface() # will guide user through the creation of blobs
blobs(['all'], 10000) # will create blobs using all variables with pop >= 10000
```

dependencies:

* pysal
* numpy
* pandas
* matplotlib.pyplot with mpl_toolkits.mplot3d
* sklearn

file structure:

* `blocks/` - Census block shapefiles (Chicago only)
* `tracts/` - Census tract shapefiles (Chicago only)
* `blobs.py` - main blobs module; includes examples at end
* `configure.py` - example code for setting up blobs to run on Census blocks
* `data311.py` - script to pull data from Plenario
* `final block data.csv` - Plenario data by Census block
* `master311.csv` - Plenario data by Census tract
* `maxp.py` - Updated version of the maxp.py module in pysal/region
* `run blobs on census blocks.py` - example code that does what it says it does
* `smoothing.py` - adventures in spatial autocorrelation

suggested route:

* make a copy of your `maxp.py` file in `pysal/region` and copy over the
updated copy here (it gives much more output)
* read in `blobs.py` down to where the examples begin
* run `interface()` to be guided through the process of creating blobs
* try the examples at the bottom of `blobs.py` to learn more about the options
* run `run blobs on census blocks.py` in its entirety to try creating 
blobs on census blocks (will take a very long time)

