#!/usr/bin/env python
import sys

import numpy
import pylab


R = 6378137.0
d2r = numpy.pi/180.


def gamma(lat):
    "International Gravity Formula"
    return 9.780490*(1. + 0.0052884*(numpy.sin(d2r*lat)**2) - 0.0000059*(numpy.sin(d2r*2.*lat)**2))


tessdata = numpy.loadtxt(sys.argv[1], unpack=True)
prismdata = numpy.loadtxt(sys.argv[2], unpack=True)

# Computation point coordinates
lonp, latp, hp = tessdata[0], tessdata[1], tessdata[2]
rp = hp + R
lonp = d2r*lonp
latp = d2r*latp

# Tesseroid and prism coordinates
rc = R
latc = 0
lonc = 0

distances = numpy.sqrt(rp**2 + rc**2 - 2*rp*rc*(numpy.sin(latp)*numpy.sin(latc) +
                        numpy.cos(latp)*numpy.cos(latc)*numpy.cos(lonp - lonc)))

distances = distances*0.001

diffs = tessdata[3:] - prismdata[3:]

# Convert potential to geoid height
#diffs[0] = diffs[0]/gamma(latp)

for d, diff in zip(distances, diffs.T):
    print d, ' '.join(map(str, diff))