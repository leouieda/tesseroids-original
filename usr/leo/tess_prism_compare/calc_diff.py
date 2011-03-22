#!/usr/bin/env python
import sys

import numpy
import pylab

R = 6378137.0
d2r = numpy.pi/180.

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

for d, diff in zip(distances, diffs.T):
    print d, diff[0], diff[1], diff[2], diff[3], diff[4], diff[5], diff[6]