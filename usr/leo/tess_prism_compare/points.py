#!/usr/bin/env python
# Print the computation points in spherical coordinatesi
import sys

import numpy

# Read the tesseroid from file
w, e, s, n, t, b, dens = numpy.loadtxt(sys.argv[1])

dz = t - b

# Make the grid 10x the tesseroid
ratio = 10.

lons = numpy.arange(ratio*w, ratio*e, 0.05*ratio*e, 'f')
lats = numpy.arange(ratio*s, ratio*n, 0.05*ratio*n, 'f')
heights = numpy.arange(0.1*dz, ratio*2.*dz, dz, 'f')

for h in heights:
    for lat in lats:
        for lon in lons:
            print lon, lat, h

#ratio = 10.
#lons = numpy.arange(ratio*w, ratio*e, 0.02*ratio*e, 'f')
#lats = numpy.arange(ratio*s, ratio*n, 0.02*ratio*n, 'f')
#for lat in lats:
    #for lon in lons:
        #print lon, lat, 0.001