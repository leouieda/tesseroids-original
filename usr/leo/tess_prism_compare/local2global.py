#!/usr/bin/env python
# Convert the g and GGT from the prisms local system to the global one
import sys

import numpy

import rotations as rot

assert len(sys.argv) == 4, "Missing cmd args. Pass lon lat height of prism in this order."

# Get the spherical coordinates of the prism from argv
lonO, latO, hO = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])

prismdata = numpy.loadtxt("prismgrav_local.txt").T
lons, lats, heights = numpy.loadtxt("tessgrav.txt", unpack=True, usecols=(0,1,2))

prismgx = prismdata[3]
prismgy = prismdata[4]
prismgz = prismdata[5]
prismgxx = prismdata[6]
prismgxy = prismdata[7]
prismgxz = prismdata[8]
prismgyy = prismdata[9]
prismgyz = prismdata[10]
prismgzz = prismdata[11]

for i, coords in enumerate(zip(lons, lats, heights)):
    lon, lat, h = coords
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
    
    print lon, lat, h, g[2], ggt[0][0], ggt[0][1], ggt[0][2], ggt[1][1], \
          ggt[1][2], ggt[2][2]