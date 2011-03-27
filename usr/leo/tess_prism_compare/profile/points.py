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
dh = (ratio*2.*dz - 0.1*dz)/(nh - 1)

lons = numpy.arange(ratio*w, ratio*e, dlon, 'f')
heights = numpy.arange(0.1*dz, ratio*2.*dz, dh, 'f')

for h in heights:
    for lon in lons:
        print lon, 0., h
    if len(lons) < nlon:
        print ratio*e, 0., h

if len(heights) < nh:
    for lon in lons:
        print lon, 0., ratio*2.*dz
    if len(lons) < nlon:
        print ratio*e, 0., ratio*2.*dz
    