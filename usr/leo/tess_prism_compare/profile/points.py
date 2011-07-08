#!/usr/bin/env python
# Print the computation points in spherical coordinates
import sys

import numpy

# Read the tesseroid from file
w, e, s, n, t, b, dens = numpy.loadtxt(sys.argv[1])
dz = t - b

# Make the grid 10x the tesseroid
ratio = 10.
nlon = 100
nh = 100
dlon = ratio*(e - w)/(nlon - 1)
hmin, hmax = 0.1*dz, ratio*2.*dz
dh = (hmax - hmin)/(nh - 1)

lons = numpy.arange(ratio*w, ratio*e, dlon, 'f')
heights = numpy.arange(hmin, hmax, dh, 'f')
angle = 1
for h in heights:
    for lon in lons:
        # use lat = lon so that it is diagonal profile.
        # need this because xy and xz are zero on the north-south profile
        print lon, angle*lon, h
        #print lon, 0, h
    if len(lons) < nlon:
        print ratio*e, angle*ratio*e, h
        #print ratio*e, 0, h

if len(heights) < nh:
    for lon in lons:
        print lon, angle*lon, hmax
        #print lon, 0, hmax
    if len(lons) < nlon:
        print ratio*e, angle*ratio*e, hmax
        #print ratio*e, 0, hmax
    
