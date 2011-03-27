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

diffgz = abs(tessdata[5] - prismdata[5])
diffgz = numpy.log(diffgz)
diffgzz = abs(tessdata[-1] - prismdata[-1])
diffgzz = numpy.log(diffgzz)

shape = (100,100)
X = numpy.reshape(lonp, shape)

Y = numpy.reshape(hp, shape)

Zdiffgz = numpy.reshape(diffgz, shape)
Zdiffgzz = numpy.reshape(diffgzz, shape)
Zdist = numpy.reshape(distances, shape)

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

levels = [1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16]

pylab.figure(figsize=(12,6))
pylab.subplot(211)
pylab.title("gz")
pylab.pcolor(X, Y, Zdiffgz, vmax=1)
pylab.colorbar()
ct_data = pylab.contour(X, Y, Zdiffgz, levels=[-2], colors='w')
pylab.setp(ct_data.collections[0], linewidth=3)
#ct_data.clabel(fmt='%d')
ct_data = pylab.contour(X, Y, Zdist, levels=levels, colors='k')
ct_data.clabel(fmt='%d')
pylab.xlim(X.min(), X.max())
pylab.ylim(Y.min(), Y.max())

pylab.subplot(212)
pylab.title("gzz")
pylab.pcolor(X, Y, Zdiffgzz, vmax=1)
pylab.colorbar()
ct_data = pylab.contour(X, Y, Zdiffgzz, levels=[-3], colors='w')
pylab.setp(ct_data.collections[0], linewidth=3)
#ct_data.clabel(fmt='%d')
ct_data = pylab.contour(X, Y, Zdist, levels=levels, colors='k')
ct_data.clabel(fmt='%d')
pylab.xlim(X.min(), X.max())
pylab.ylim(Y.min(), Y.max())


pylab.savefig("comparison.png", dpi=200)
#pylab.savefig("comparison.pdf")

pylab.show()