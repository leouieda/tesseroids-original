#!/usr/bin/env python
import sys

import numpy
import pylab

data = numpy.loadtxt(sys.argv[1], unpack=True)
w, e, s, n, t, b, dens = numpy.loadtxt(sys.argv[2])
R = 6378137.0
d2r = numpy.pi/180.

size = max([0.001*R*d2r*(e-w), 0.001*R*d2r*(n-s), 0.001*(t-b)])

distances = data[0]/size
diffs = data[1:]

#labels = [r'$g_z$', r'$g_{xx}$', r'$g_{xy}$', r'$g_{xz}$', r'$g_{yy}$', r'$g_{yz}$', r'$g_{zz}$']
labels = ['gz', 'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']

noise_si = 10**(-1)
noise_mgal = 10**(-2)
noise_eotvos = 10**(-3)

pylab.figure(figsize=(10,9))
pylab.subplots_adjust(top=0.95, bottom=0.08, left=0.1, right=0.95,
                      hspace=0.5, wspace=0.25)

fontsize = 11

#pylab.subplot(4,2,1)
#pylab.title('Potential', fontsize=12)
#pylab.plot(distances, numpy.abs(diffs[0]), '.k')
#pylab.plot([0, 1.1*distances.max()], [noise_si, noise_si], '-r', linewidth=1.5)
##pylab.xscale('log')
#pylab.yscale('log')
#pylab.xlim(0, 1.1*distances.max())
#pylab.ylim(10**(-11), 10)
#pylab.grid(True)
##pylab.gca().xaxis.grid(True, which='minor')
#pylab.gca().xaxis.set_ticks(numpy.arange(0, 1.1*distances.max(), 2))
#pylab.gca().yaxis.set_ticks([10**(i) for i in range(-11, 2, 2)])
#pylab.xlabel("Size Ratio", fontsize=fontsize)
#pylab.ylabel(r"Difference [SI]", fontsize=fontsize)

pylab.subplot(4,1,1)
pylab.title('gz', fontsize=12)
pylab.plot(distances, numpy.abs(diffs[1]), '.k')
pylab.plot([0, 1.1*distances.max()], [noise_mgal, noise_mgal], '-r', linewidth=1.5)
#pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(0, 1.1*distances.max())
pylab.ylim(10**(-11), 10)
pylab.grid(True)
#pylab.gca().xaxis.grid(True, which='minor')
pylab.gca().xaxis.set_ticks(numpy.arange(0, 1.1*distances.max(), 1))
pylab.gca().yaxis.set_ticks([10**(i) for i in range(-11, 2, 2)])
pylab.xlabel("Size Ratio", fontsize=fontsize)
pylab.ylabel(r"Difference [mGal]", fontsize=fontsize)

for i, tmp in enumerate(zip(diffs[2:], labels[1:])):
    diff, label = tmp
    pylab.subplot(4,2,i+3)
    pylab.title(label, fontsize=12)
    pylab.plot(distances, numpy.abs(diff), '.k')
    pylab.plot([0, 1.1*distances.max()], [noise_eotvos, noise_eotvos], '-r', linewidth=1.5)
    #pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlim(0, 1.1*distances.max())
    pylab.ylim(10**(-11), 10)
    pylab.grid(True)
    #pylab.gca().xaxis.grid(True, which='minor')
    pylab.gca().xaxis.set_ticks(numpy.arange(0, 1.1*distances.max(), 2))
    pylab.gca().yaxis.set_ticks([10**(i) for i in range(-11, 2, 2)])
    pylab.xlabel("Size Ratio", fontsize=fontsize)
    pylab.ylabel(r"Difference [Eotvos]", fontsize=fontsize)
    
pylab.savefig("comparison.png", dpi=200)
#pylab.savefig("comparison.pdf")

#pylab.show()