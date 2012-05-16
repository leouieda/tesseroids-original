import sys
from math import pi
import numpy
from matplotlib import pyplot

lon, lat, height, gxx, gyy, gzz = numpy.loadtxt(sys.argv[1], unpack=True)
trace = gxx + gyy + gzz
inside = 4*pi*numpy.ones_like(trace)
xmin, xmax = trace.min(), trace.max()
pyplot.figure(figsize=(8,6))
pyplot.plot(trace, height, '-k', linewidth=2, label="Trace")
#pyplot.plot(inside, height, '--r')
pyplot.hlines([0, -5000], xmin, xmax, 'r', '--', linewidth=2,
    label="Top and bottom of tesseroid")
pyplot.legend(loc='upper left')
pyplot.xlabel("Trace = gxx + gyy + gzz")
pyplot.ylabel("Height [m]")
pyplot.xlim(xmin, xmax)
pyplot.ylim(height.min(), height.max())
pyplot.savefig('output.png')
#pyplot.show()
