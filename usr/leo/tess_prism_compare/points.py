#!/usr/bin/env python
# Print the computation points in spherical coordinates
import numpy

lons = numpy.arange(-1, 1, 0.05)
lats = numpy.arange(-1, 1, 0.05)
heights = numpy.arange(1, 250000, 5000)

#for h in heights:
    #for lat in lats:
        #for lon in lons:
            #print lon, lat, h

for lat in lats:
    for lon in lons:
        print lon, lat, 20000