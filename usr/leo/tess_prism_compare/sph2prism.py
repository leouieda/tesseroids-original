#!/usr/bin/env python
# Convert the spherical coordinates to local coordinates of the prism
import sys
import rotations

# Get the spherical coordinates of the prism 
lonO, latO, hO = 0, 0, 0

f = sys.stdin
line = f.readline()
while line:
    if line[0] == '#' or len(line.split()) < 3:
        line = f.readline()
        continue
    lon, lat, h = [float(n) for n in line.split()[:3]]
    x, y, z = rotations.G2L(lon, lat, h, lonO, latO, hO)
    print y, x, z
    line = f.readline()