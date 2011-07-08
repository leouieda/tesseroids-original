#!/usr/bin/env python
import sys

import matplotlib
import numpy
import pylab

tessdata = numpy.loadtxt(sys.argv[2], unpack=True)
prismdata = numpy.loadtxt(sys.argv[1], unpack=True)
w, e, s, n, t, b, dens = numpy.loadtxt(sys.argv[3])

R = 6378137.0
d2r = numpy.pi/180.

lonp, latp, hp = tessdata[0], tessdata[1], tessdata[2]
rp = hp + R
lonp = d2r*lonp
latp = d2r*latp

rc = t + R
lonc = d2r*0.5*(e + w)
latc = d2r*0.5*(n + s)

size = max([0.001*R*d2r*(e-w), 0.001*R*d2r*(n-s), 0.001*(t-b)])
distances = 0.001*numpy.sqrt(rp**2 + rc**2 - 2*rp*rc*(numpy.sin(latp)*numpy.sin(latc) +
                 numpy.cos(latp)*numpy.cos(latc)*numpy.cos(lonp - lonc)))
distances = distances/size

# Need to use base 10 log!
diffs = numpy.log10(abs(tessdata[4:] - prismdata[4:]))

shape = (100,100)
X = numpy.reshape(lonp, shape)
Y = numpy.reshape(hp, shape)
Zdist = numpy.reshape(distances, shape)

components = ['gz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
errors = [-2, -3, -3, -3, -3, -3, -3]

pylab.figure(figsize=(14,9))
pylab.subplots_adjust(hspace=0.5, wspace=0.2, top=0.95, bottom=0.08)
# Plot gz
comp, error = components[0], errors[0]
Zdiff = numpy.reshape(diffs[0], shape)
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
levels = [1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16]
pylab.subplot(4,2,1)
pylab.title(comp)
pylab.pcolor(X, Y, Zdiff)
pylab.colorbar()
ct_data = pylab.contour(X, Y, Zdiff, levels=[error], colors='w')
pylab.setp(ct_data.collections[0], linewidth=3)
ct_data = pylab.contour(X, Y, Zdist, levels=levels, colors='k')
ct_data.clabel(fmt='%d')
pylab.xlim(X.min(), X.max())
pylab.ylim(Y.min(), Y.max())
pylab.xlabel(r"Longitude = Latitude ($^\circ$)")
pylab.ylabel("Height (m)")
# Plot tensor
for i in xrange(1, len(components)):
    comp, error = components[i], errors[i]
    Zdiff = numpy.reshape(diffs[i], shape)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    levels = [1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16]
    pylab.subplot(4,2,i+2)
    pylab.title(comp)
    pylab.pcolor(X, Y, Zdiff)
    pylab.colorbar()
    ct_data = pylab.contour(X, Y, Zdiff, levels=[error], colors='w')
    pylab.setp(ct_data.collections[0], linewidth=3)
    ct_data = pylab.contour(X, Y, Zdist, levels=levels, colors='k')
    ct_data.clabel(fmt='%d')
    pylab.xlim(X.min(), X.max())
    pylab.ylim(Y.min(), Y.max())
    pylab.xlabel(r"Longitude = Latitude ($^\circ$)")
    pylab.ylabel("Height (m)")

pylab.savefig("comparison.png", dpi=200)
#pylab.savefig("comparison.pdf")

pylab.show()
