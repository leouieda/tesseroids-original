# Make a plot of prism_output.txt

import pylab
import matplotlib

output = pylab.loadtxt('prism_output.txt').T

shape = (100, 100)
X = pylab.reshape(output[0], shape)*0.001
Y = pylab.reshape(output[1], shape)*0.001

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

pylab.figure(figsize=(9,10))
pylab.subplots_adjust(hspace=0.35, wspace=0.01, bottom=0.06, top=0.95, left=0.05,
                      right=1)
sp = [7, 1, 2, 3, 5, 6, 9]
for i in range(3, len(output)):
    Z = pylab.reshape(output[i], shape)
    pylab.subplot(3, 3, sp[i - 3])
    pylab.axis("scaled")
    pylab.title("prism_output.txt: column %d" % (i + 1), fontsize=10)
    #pylab.contourf(X, Y, Z, 8)
    pylab.pcolor(X, Y, Z)
    ct = pylab.contour(X, Y, Z, 8, colors='k')
    ct.clabel(fmt='%g', fontsize=11)
    pylab.xlabel("Easting (km)", fontsize=10)
    pylab.ylabel("Northing (km)", fontsize=10)
    pylab.xlim(X.min(), X.max())
    pylab.ylim(Y.min(), Y.max())

pylab.show()