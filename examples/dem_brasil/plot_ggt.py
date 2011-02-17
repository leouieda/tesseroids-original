"""
Make plots of the DEM  calculated gravity effect using Matplotlib with Basemap
toolkit
"""

import numpy
import pylab
from mpl_toolkits.basemap import Basemap

data = pylab.loadtxt('parana-ggt-10sec.txt')
lons, lats, heights, gxx, gxy, gxz, gyy, gyz, gzz = data.T
nlons = 100  # Number of points in the longitude direction
nlats = len(lats)/nlons

# Convert the lists to 2D grids
glons = numpy.reshape(lons, (nlats, nlons))
glats = numpy.reshape(lats, (nlats, nlons))

fig = pylab.figure(figsize=(14,9))
pylab.subplots_adjust(wspace=0.35)#hspace=0.15,

bm = Basemap(projection='merc', \
             llcrnrlon=lons[0], llcrnrlat=lats[0], \
             urcrnrlon=lons[-1], urcrnrlat=lats[-1], \
             lon_0=lons[nlons//2], lat_0=lats[len(lats)//2],
             resolution='l', area_thresh=10000)
glons, glats = bm(glons, glats)
             
titles = [r"$g_{xx}$", r"$g_{xy}$", r"$g_{xz}$", r"$g_{yy}$", r"$g_{yz}$",
          r"$g_{zz}$"]
fields = [gxx, gxy, gxz, gyy, gyz, gzz]

for i, args in enumerate(zip(fields, titles)):

    field, title = args

    ax = pylab.subplot(2, 3, i + 1, aspect='equal')
    pylab.title(title, fontsize=18)

    # Make it a 2D grid
    gfield = numpy.reshape(field, (nlats, nlons))

    # Plot the coastlines and etc
    bm.drawmeridians(numpy.arange(lons[0]+2.5, lons[-1]-2.5, 2.5), 
                     labels=[0,0,0,1], fontsize=9, linewidth=0.5)#, 
                    #labelstyle='+/-')
    bm.drawparallels(numpy.arange(lats[0]+2.5, lats[-1]-2.5, 2.5),
                     labels=[1,0,0,0], fontsize=9, linewidth=0.5)#, 
                     #labelstyle='+/-')
    bm.drawcoastlines(linewidth=1)
    #bm.fillcontinents(color='coral')
    bm.drawmapboundary()
    #bm.bluemarble()
    bm.drawcountries(linewidth=1)
    bm.drawstates(linewidth=0.2)

    # Plot it
    cf = bm.pcolor(glons, glats, gfield, cmap=pylab.cm.jet)
    cb = pylab.colorbar(orientation='vertical', format='%.2f', shrink=0.73)
    #cb.set_label("Eotvos")
    
pylab.savefig('parana-ggt-10sec.png', fmt='png')
