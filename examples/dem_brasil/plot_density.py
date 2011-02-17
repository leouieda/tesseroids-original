"""
Make plots of the DEM density using Matplotlib with Basemap toolkit
"""

import sys
import numpy
import pylab
from mpl_toolkits.basemap import Basemap

lons, lats, dens = pylab.loadtxt('parana-dem-10sec-density.txt', unpack=True, usecols=(0,1,3))
nlons = 151  # Number of points in the longitude direction
nlats = len(lats)/nlons

# Convert from SI to g cm^-3
dens = 0.001*dens

# Convert the lists to 2D grids
glons = numpy.reshape(lons, (nlats, nlons))
glats = numpy.reshape(lats, (nlats, nlons))
gdens = numpy.reshape(dens, (nlats, nlons))

pylab.figure()

# Plot with Mercator projection
bm = Basemap(projection='merc',
             llcrnrlon=lons[0], llcrnrlat=lats[-1],
             urcrnrlon=lons[-1], urcrnrlat=lats[0],
             lon_0=lons[nlons//2], lat_0=lats[len(lats)//2],
             resolution='l',
             area_thresh=10000)

bm.drawmeridians(numpy.arange(lons[0]+5., lons[-1], 5.),
                 labels=[0,0,0,1], fontsize=12, linewidth=0.5)#,
                 #labelstyle='+/-')
bm.drawparallels(numpy.arange(lats[-1]+5., lats[0], 5.),
                 labels=[1,0,0,0], fontsize=12, linewidth=0.5)#,
                 #labelstyle='+/-')
bm.drawcoastlines(linewidth=1)
#bm.fillcontinents(color='coral')
bm.drawmapboundary()
#bm.bluemarble()
bm.drawcountries(linewidth=0.8)
#bm.drawstates(linewidth=0.6)

# Do the pseudo-color plot
glons, glats = bm(glons, glats)
cf = bm.pcolor(glons, glats, gdens, cmap=pylab.cm.jet)
cb = pylab.colorbar()
cb.set_label(r"Density [$g.cm^{-3}$]")

pylab.savefig('parana-dem-10sec-density.png')