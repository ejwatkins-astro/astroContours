To use contour codes, download all of astroContours into a folder called astroContours

Then in the same directory as the folder location, you can import astroContours

Code requirements
Python 3

Package requirements likely already met:
numpy
scipy
matplotlib
astropy
skimage

Package requirements I expect you'll need to install:

geopandas: awesome package for comparing, joining, analysing and visulising large shape catalogues (including polygons created from a mask catalogue

shapely:   package to convert contours into polygon objects for geopandas

descartes: package for converting shapley objects into matplotib patches

colorcet:  A huge list of perceptually uniform matplotlib colourmaps

regions:   astropy affiliated package for creating, opening, analysing and visualising and ds9 regions

pickle:    a way to save data

json       a way to save database data (usually dicts) in human readable format. It normally prefered over pickle for security reasons.


