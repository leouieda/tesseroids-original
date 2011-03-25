# Make a plot of the data file passed

import sys
import pylab
import numpy

d2r = numpy.pi/180.


def gamma(lat):
    "International Gravity Formula"
    return 9.780490*(1. + 0.0052884*(numpy.sin(d2r*lat)**2) - 0.0000059*(numpy.sin(d2r*2.*lat)**2))


assert len(sys.argv) > 4, "Need 4 arguments: filename, column, nlon, nlat"

fname = sys.argv[1]
col = int(sys.argv[2])
nlons = int(sys.argv[3])
nlats = int(sys.argv[4])

lons, lats, vals = pylab.loadtxt(fname, usecols=[0,1,col-1], unpack=True)

vals = vals/gamma(lats)

X = pylab.reshape(lons, (nlons, nlats))
Y = pylab.reshape(lats, (nlons, nlats))
Z = pylab.reshape(vals, (nlons, nlats))

pylab.figure()
pylab.axis("scaled")
pylab.title("%s: column %d" % (fname, col))
#pylab.contourf(X, Y, Z, 12)
pylab.pcolor(X, Y, Z)
pylab.colorbar()
ct_data = pylab.contour(X, Y, Z, 12, colors='k')
#ct_data = pylab.contourf(X, Y, Z, 12)
#ct_data = pylab.contour(X, Y, Z, ct_data.levels, colors='k')
ct_data.clabel(fmt='%g', fontsize=12)
pylab.xlabel("Longitude")
pylab.ylabel("Latitude")
pylab.xlim(X.min(), X.max())
pylab.ylim(Y.min(), Y.max())
pylab.show()