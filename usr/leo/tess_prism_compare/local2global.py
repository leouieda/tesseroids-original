#!/usr/bin/env python
# Convert the g and GGT from the prisms local system to the global one
import sys

import numpy

import rotations as rot

# Get the spherical coordinates of the prism
lonO, latO, hO = 0, 0, 0

prismdata = numpy.loadtxt(sys.argv[1]).T
lons, lats, heights = numpy.loadtxt(sys.argv[2], unpack=True, usecols=(0,1,2))

prismpot = prismdata[3]
prismgx = prismdata[4]
prismgy = prismdata[5]
prismgz = prismdata[6]
prismgxx = prismdata[7]
prismgxy = prismdata[8]
prismgxz = prismdata[9]
prismgyy = prismdata[10]
prismgyz = prismdata[11]
prismgzz = prismdata[12]

for i, coords in enumerate(zip(lons, lats, heights)):
    lon, lat, h = coords
    pot = prismpot[i]
    g = numpy.matrix([prismgx[i], prismgy[i], prismgz[i]]).T
    ggt = numpy.matrix([[prismgxx[i], prismgxy[i], prismgxz[i]],
                        [prismgxy[i], prismgyy[i], prismgyz[i]],
                        [prismgxz[i], prismgyz[i], prismgzz[i]]])
    # Rotate g and ggt to the global system
    trans = rot.P2()*rot.R2(lat - 90.)*rot.R3(lon - lonO)*rot.R2(90. - latO)*rot.P2()

    g = trans*g
    ggt = trans*ggt*(trans.T)

    ggt = ggt.tolist()
    g = g.T.tolist()[0]

    # xz and yz should have their signs changed because prism has z->down and tess has z->up
    print lon, lat, h, pot, g[2], ggt[0][0], ggt[0][1], -1*ggt[0][2], ggt[1][1], \
          -1*ggt[1][2], ggt[2][2]