"""
Assign density values for the DEM points. 2.67 for continent (h>=0) and 1.67
for ocean (h<0).
"""

import numpy

lons, lats, heights = numpy.loadtxt('parana-dem-10sec.xyz', unpack=True)

for i in xrange(len(heights)):
    if heights[i] >=0:
        print "%lf %lf %lf %lf" % (lons[i], lats[i], heights[i], 2670.0)
    else:
        print "%lf %lf %lf %lf" % (lons[i], lats[i], heights[i], 1670.0)
