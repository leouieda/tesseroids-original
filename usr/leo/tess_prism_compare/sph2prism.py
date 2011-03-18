#!/usr/bin/env python
# Convert the spherical coordinates to local coordinates of the prism
import sys
import rotations

assert len(sys.argv) == 4, "Missing cmd args. Pass lon lat height of prism in this order."

# Get the spherical coordinates of the prism from argv
lonO, latO, hO = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])

f = sys.stdin
line = f.readline()
while line:
    if line[0] == '#' or len(line.split()) < 3:
        line = f.readline()
        continue
    lon, lat, h = [float(n) for n in line.split()[:3]]
    x, y, z = rotations.G2L(lon, lat, h, lonO, latO, hO)
    print x, y, -1*z
    line = f.readline()